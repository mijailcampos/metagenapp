import os
import subprocess
import shutil
import pickle
import time
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd
from metagenapp.metagen_config import NAIVE_MODEL_PATH, SILVA_TRAINSET_TAX


# ============================================================
# 2. NAIVE BAYES (Python) — Block Classification
# ============================================================
def clasificar_una_secuencia(args):
    sid, seq, cond_probs, taxa, k = args

    best_tax = None
    best_score = -1e99

    for taxon in taxa:
        score = 0
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            score += cond_probs[taxon].get(kmer, 0)
        if score > best_score:
            best_score = score
            best_tax = taxon

    return sid, best_tax


def classify_naive_por_bloques(
    fasta_path,
    output_path,
    modelo_path=NAIVE_MODEL_PATH,
    block_size=10000,
    n_threads=None   # 👈 se acepta pero NO se usa
):
    import os, time, pickle, gc
    from Bio import SeqIO

    if n_threads not in (None, 1):
        print("🧠 Naive Bayes classification (Linux-safe mode)")

    # 🔒 Cargar modelo UNA sola vez
    with open(modelo_path, "rb") as f:
        payload = pickle.load(f)

    # 🔬 Modelo antiguo (tuple)
    if isinstance(payload, tuple):
        cond_probs, taxa, k = payload

    # 🔬 Modelo PR2 nuevo (dict)
    elif isinstance(payload, dict):
        class_kmer_counts = payload.get("class_kmer_counts") or payload.get("kmer_counts")
        #class_total_kmers = payload["class_total_kmers"]
        k = payload["k"]

        # reconstruimos estructura compatible
        taxa = list(class_kmer_counts.keys())
        cond_probs = class_kmer_counts

    else:
        raise ValueError("Unknown model format")


    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    open(output_path, "w").close()

    secuencias = SeqIO.parse(fasta_path, "fasta")
    bloque = []
    total = 0
    bloque_id = 0

    start = time.time()

    for record in secuencias:
        bloque.append((record.id, str(record.seq).upper()))

        if len(bloque) == block_size:
            bloque_id += 1
            print(f"🟦 Processing block {bloque_id} ({len(bloque)} sequences)")

            resultados = [
                clasificar_una_secuencia((sid, seq, cond_probs, taxa, k))
                for sid, seq in bloque
            ]

            with open(output_path, "a") as fout:
                for sid, taxon in resultados:
                    fout.write(f"{sid}\t{taxon}\n")

            total += len(bloque)
            bloque.clear()
            del resultados
            gc.collect()  # 🔑 clave en Linux

    if bloque:
        bloque_id += 1
        print(f"🟩 Processing final block ({len(bloque)} sequences)")

        resultados = [
            clasificar_una_secuencia((sid, seq, cond_probs, taxa, k))
            for sid, seq in bloque
        ]

        with open(output_path, "a") as fout:
            for sid, taxon in resultados:
                fout.write(f"{sid}\t{taxon}\n")

        total += len(bloque)

    end = time.time()
    print(f"🏁 Classification complete: {total} sequences in {end - start:.2f}s")

    return output_path, None

# ============================================================
# 2B. NAIVE BAYES (Sklearn Incremental Model) — PRO Mode
# ============================================================

def classify_naive_sklearn(
    fasta_path,
    output_path,
    modelo_path,
    block_size=5000
):
    import joblib
    import gc
    from Bio import SeqIO
    import os
    import time

    print("🧠 Naive Bayes PRO (sklearn incremental model)")

    payload = joblib.load(modelo_path)
    clf = payload["model"]
    vectorizer = payload["vectorizer"]

    dirpath = os.path.dirname(output_path)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)

    open(output_path, "w").close()

    bloque = []
    total = 0
    bloque_id = 0

    start = time.time()

    for record in SeqIO.parse(fasta_path, "fasta"):
        bloque.append((record.id, str(record.seq).upper()))

        if len(bloque) == block_size:
            bloque_id += 1
            print(f"🟦 Processing block {bloque_id} ({len(bloque)} sequences)")

            ids = [x[0] for x in bloque]
            seqs = [x[1] for x in bloque]

            X = vectorizer.transform(seqs)
            preds = clf.predict(X)

            with open(output_path, "a") as fout:
                for sid, taxon in zip(ids, preds):
                    fout.write(f"{sid}\t{taxon}\n")

            total += len(bloque)
            bloque.clear()
            del X, preds
            gc.collect()

    if bloque:
        bloque_id += 1
        print(f"🟩 Processing final block ({len(bloque)} sequences)")

        ids = [x[0] for x in bloque]
        seqs = [x[1] for x in bloque]

        X = vectorizer.transform(seqs)
        preds = clf.predict(X)

        with open(output_path, "a") as fout:
            for sid, taxon in zip(ids, preds):
                fout.write(f"{sid}\t{taxon}\n")

        total += len(bloque)

    end = time.time()
    print(f"🏁 PRO Classification complete: {total} sequences in {end - start:.2f}s")

    return output_path, None


# ============================================================
# 3. VSEARCH SINTAX — Fast Classifier
# ============================================================
def classify_vsearch_sintax(fasta_path, reference_fasta, output_path, n_threads=8):
    import subprocess
    import re

    fasta_path = str(fasta_path)
    reference_fasta = str(reference_fasta)
    output_path = str(output_path)

    raw_output = output_path + ".raw"

    cmd = [
        "vsearch",
        "--sintax", fasta_path,
        "--db", reference_fasta,
        "--tabbedout", raw_output,
        "--sintax_cutoff", "0.8",
        "--strand", "both",
        "--threads", str(n_threads),
    ]

    subprocess.run(cmd, check=True)

    with open(raw_output, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = [p.strip() for p in line.split("\t")]

            # Primer campo = seq_id
            seq_id = parts[0]

            # Buscar campo de taxonomía útil
            taxonomy = ""
            for p in parts[1:]:
                if "tax=" in p:
                    # Si alguna vez aparece formato con tax=
                    taxonomy = p.split("tax=", 1)[1]
                    break
                elif any(prefix in p for prefix in ["d:", "k:", "p:", "c:", "o:", "f:", "g:", "s:"]):
                    taxonomy = p

            if not taxonomy or taxonomy == "+":
                taxonomy = "unknown;"

            # Limpiar scores/confidencias tipo (0.98)
            taxonomy = re.sub(r"\([0-9.]+\)", "", taxonomy)

            # Quitar prefijos de nivel
            taxonomy = (
                taxonomy.replace("d:", "")
                        .replace("k:", "")
                        .replace("p:", "")
                        .replace("c:", "")
                        .replace("o:", "")
                        .replace("f:", "")
                        .replace("g:", "")
                        .replace("s:", "")
            )

            # Convertir comas en ;
            taxonomy = taxonomy.replace(",", ";")

            # Limpiar separadores repetidos
            taxonomy = re.sub(r";+", ";", taxonomy).strip(";")

            if not taxonomy:
                taxonomy = "unknown"

            taxonomy = taxonomy + ";"

            # ESCRIBIR SOLO 2 COLUMNAS
            fout.write(f"{seq_id}\t{taxonomy}\n")

    return output_path

# ============================================================
# 4. Convert SINTAX → Mothur .taxonomy format
# ============================================================
def convertir_sintax_a_mothur(tabbed_file, salida_final):

    df = pd.read_csv(
        tabbed_file,
        sep="\t",
        header=None,
        names=["seq_id", "sintax"]
    )

    def extract_taxonomy(s):
        if "tax=" not in s:
            return "Unclassified"

        taxpart = s.split("tax=")[1].split(";")[0]
        levels = taxpart.split(",")

        names = []
        for n in levels:
            if ":" in n:
                names.append(n.split(":")[1])
        return ";".join(names)

    df["taxonomy"] = df["sintax"].apply(extract_taxonomy)
    df = df[["seq_id", "taxonomy"]]

    df.to_csv(salida_final, sep="\t", header=False, index=False)


# ============================================================
# 5. ALIGNMENT FALLBACK
# Upgrades features that didn't reach genus level by using
# the best edlib alignment hit against the SILVA reference.
# Only called in ref mode after naive-v2 classification.
# ============================================================

def _taxonomy_depth(tax_str: str) -> int:
    """Count non-empty levels in a semicolon-separated taxonomy string."""
    return sum(1 for p in tax_str.split(";") if p.strip())


def apply_alignment_fallback(
    taxonomy_file: str,
    edlib_tsv: str,
    silva_tax_file: str = None,
    min_identity: float = 0.94,
    genus_level: int = 6,
) -> dict:
    """
    For features classified below genus level, attempt to assign taxonomy
    from the best edlib alignment hit if identity >= min_identity.

    Works for any 16S sample type (pharyngeal, gut, skin, environmental)
    because SILVA NR99 covers all environments.

    min_identity=0.94 is appropriate for genus-level assignment in 16S
    (97% = species, 94-95% = genus, per Yarza et al. 2014).

    Modifies taxonomy_file in-place.
    Returns stats dict with keys: n_total, n_below_genus, n_upgraded.
    """
    if silva_tax_file is None:
        silva_tax_file = str(SILVA_TRAINSET_TAX)

    if not os.path.exists(silva_tax_file):
        print(f"   [fallback] WARNING: SILVA taxonomy not found at {silva_tax_file} — skipping fallback")
        return {"n_total": 0, "n_below_genus": 0, "n_upgraded": 0}

    if not os.path.exists(edlib_tsv):
        print(f"   [fallback] WARNING: edlib TSV not found — skipping fallback")
        return {"n_total": 0, "n_below_genus": 0, "n_upgraded": 0}

    # Load SILVA taxonomy indexed by bare accession (no version, no region coords).
    # SILVA keys look like "AY727530_1_1478"; edlib ref_ids look like "AY727530.1".
    # We strip to the raw accession for matching.
    print("   [fallback] Loading SILVA taxonomy index...")
    silva_tax = {}   # accession -> taxonomy string up to genus
    with open(silva_tax_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            raw_key = parts[0]                          # e.g. AY727530_1_1478
            accession = raw_key.split("_")[0]           # e.g. AY727530
            tax = parts[1].rstrip(";")
            levels = [p.strip() for p in tax.split(";") if p.strip()]
            if len(levels) >= genus_level and accession not in silva_tax:
                silva_tax[accession] = ";".join(levels[:genus_level])
    print(f"   [fallback] {len(silva_tax):,} accessions indexed")

    # Load edlib TSV: query -> (ref_id, identity)  — keep best hit per query
    print("   [fallback] Loading edlib alignment hits...")
    best_hit = {}
    with open(edlib_tsv) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            query, ref_id, identity = parts[0], parts[1], float(parts[2])
            if identity >= min_identity:
                if query not in best_hit or identity > best_hit[query][1]:
                    best_hit[query] = (ref_id, identity)
    print(f"   [fallback] {len(best_hit):,} hits above identity {min_identity}")

    # Read current taxonomy, upgrade where needed
    entries = []
    with open(taxonomy_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            seq_id = parts[0]
            tax    = parts[1] if len(parts) > 1 else "Unclassified"
            entries.append([seq_id, tax])

    n_total       = len(entries)
    n_below_genus = 0
    n_upgraded    = 0

    for entry in entries:
        seq_id, tax = entry
        if _taxonomy_depth(tax) < genus_level:
            n_below_genus += 1
            hit = best_hit.get(seq_id)
            if hit:
                ref_id, identity = hit
                # edlib ref_id is like "AY727530.1" — strip version suffix to get accession
                accession = ref_id.split(".")[0]
                fallback_tax = silva_tax.get(accession)
                if fallback_tax:
                    entry[1] = fallback_tax
                    n_upgraded += 1

    # Write back
    with open(taxonomy_file, "w") as fh:
        for seq_id, tax in entries:
            fh.write(f"{seq_id}\t{tax}\n")

    print(f"   [fallback] {n_below_genus} features below genus | "
          f"{n_upgraded} upgraded ({100*n_upgraded/n_total:.1f}% of total)")

    return {"n_total": n_total, "n_below_genus": n_below_genus, "n_upgraded": n_upgraded}
