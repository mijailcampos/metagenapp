import os
import subprocess
import shutil
import pickle
import time
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd


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
    modelo_path="modelos/naive_model.pkl",
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
        class_kmer_counts = payload["class_kmer_counts"]
        class_total_kmers = payload["class_total_kmers"]
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

    cmd = [
        "vsearch",
        "--sintax", fasta_path,
        "--db", reference_fasta,
        "--tabbedout", output_path,
        "--sintax_cutoff", "0.6",
        "--threads", str(n_threads)
    ]

    try:
        subprocess.run(cmd, check=True)
        return output_path, None

    except subprocess.CalledProcessError as e:
        return None, f"Error running VSEARCH:\n{e}"


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
