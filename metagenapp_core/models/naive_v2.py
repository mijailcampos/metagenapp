MODEL_NAME = "naive_v2"
MODEL_TYPE = "taxonomy_classifier"

import os
from metagenapp_core.utils.parallel import run_parallel

# ============================================================
# GLOBALS para Naive V2 Parallel
# ============================================================

_V2_MODEL = None

def _v2_init_worker(model):
    import random
    global _V2_MODEL
    _V2_MODEL = model
    random.seed()  # re-seed desde OS entropy — evita que workers fork compartan el mismo estado PRNG

def _v2_worker(pair):
    from metagenapp_core.models.naive_v2_engine import classify_seq_v2
    sid, seq = pair
    return classify_seq_v2(sid, seq, _V2_MODEL)

# ============================================================
# 2. NAIVE BAYES (Python) — Block Classification (FAST + PARALLEL)
# ============================================================

import multiprocessing as mp
from collections import Counter

# --- globals para workers (Linux fork: comparten memoria por COW) ---
_NAIVE = {"cond_probs": None, "taxa": None, "k": None}

def _naive_init_worker(cond_probs, taxa, k):
    _NAIVE["cond_probs"] = cond_probs
    _NAIVE["taxa"] = taxa
    _NAIVE["k"] = k

def _naive_worker(pair):
    sid, seq = pair
    return clasificar_una_secuencia((sid, seq, _NAIVE["cond_probs"], _NAIVE["taxa"], _NAIVE["k"]))


def clasificar_una_secuencia(args):
    """
    Clasificación naive Bayes "manual" con:
    - Precomputo de k-mers UNA sola vez por secuencia
    - Conteo de k-mers (Counter)
    - Suma ponderada por frecuencia de k-mers
    """
    sid, seq, cond_probs, taxa, k = args

    n = len(seq) - k + 1
    if n <= 0:
        return sid, None

    # k-mers UNA sola vez
    kmer_counts = Counter(seq[i:i+k] for i in range(n))
    items = list(kmer_counts.items())  # lista reutilizable

    best_tax = None
    best_score = float("-inf")

    for taxon in taxa:
        probs = cond_probs[taxon]
        getp = probs.get
        score = 0.0

        # iteras solo sobre k-mers únicos
        for kmer, c in items:
            w = getp(kmer, 0)
            if w:
                score += w * c

        if score > best_score:
            best_score = score
            best_tax = taxon

    return sid, best_tax


def classify_naive_por_bloques(
    fasta_path,
    output_path,
    modelo_path="modelos/naive_model.pkl",
    block_size=10000,
    n_threads=None,
):
    """
    Clasifica FASTA por bloques con multiprocessing (Linux fork).
    - Carga modelo 1 vez
    - Usa Pool con initializer para no pasar el modelo por tarea
    - Escribe resultados incrementalmente
    """
    import os, time, pickle, gc
    from Bio import SeqIO

    # defaults sensatos
    if n_threads is None:
        n_threads = min(12, os.cpu_count() or 8)

    print(f"🧠 Naive Bayes classification (fork pool) | threads={n_threads} | block={block_size}")

    # Cargar modelo UNA sola vez
    with open(modelo_path, "rb") as f:
        payload = pickle.load(f)

    # Modelo antiguo (tuple)
    if isinstance(payload, tuple):
        cond_probs, taxa, k = payload

    # Modelo PR2 nuevo (dict)
    elif isinstance(payload, dict):
        class_kmer_counts = payload.get("class_kmer_counts") or payload.get("kmer_counts")
        k = payload["k"]
        taxa = list(class_kmer_counts.keys())
        cond_probs = class_kmer_counts

    else:
        raise ValueError("Unknown model format")

    # preparar salida
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    open(output_path, "w").close()

    # Pool fork (CRÍTICO en Linux para no duplicar RAM)
    ctx = mp.get_context("fork")
    pool = ctx.Pool(
        processes=n_threads,
        initializer=_naive_init_worker,
        initargs=(cond_probs, taxa, k),
        maxtasksperchild=200,  # ayuda contra fragmentación/fugas
    )

    total = 0
    bloque_id = 0
    start = time.time()

    try:
        bloque = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            bloque.append((record.id, str(record.seq).upper()))

            if len(bloque) >= block_size:
                bloque_id += 1
                print(f"🟦 Processing block {bloque_id} ({len(bloque)} sequences)")

                resultados_iter = pool.imap_unordered(_naive_worker, bloque, chunksize=200)

                with open(output_path, "a") as fout:
                    for sid, taxon in resultados_iter:
                        fout.write(f"{sid}\t{taxon}\n")

                total += len(bloque)
                bloque.clear()
                gc.collect()

        # bloque final
        if bloque:
            bloque_id += 1
            print(f"🟩 Processing final block ({len(bloque)} sequences)")

            resultados_iter = pool.imap_unordered(_naive_worker, bloque, chunksize=200)

            with open(output_path, "a") as fout:
                for sid, taxon in resultados_iter:
                    fout.write(f"{sid}\t{taxon}\n")

            total += len(bloque)

    finally:
        pool.close()
        pool.join()

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



# ============================================================
# 5. NAIVE BAYES V2 — Inverted Index Engine
# ============================================================

def classify_naive_v2_por_bloques(
    fasta_path,
    output_path,
    modelo_path="/data/databases/metagenapp_refs/16S/naive_model_v2.pkl",
    block_size=10000,
):
    import os, time, pickle, gc
    from Bio import SeqIO
    from metagenapp_core.models.naive_v2_engine import classify_seq_v2

    print(f"🧠 Naive V2 (inverted index) | block={block_size}")

    with open(modelo_path, "rb") as f:
        model_v2 = pickle.load(f)

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    open(output_path, "w").close()

    total = 0
    bloque_id = 0
    start = time.time()

    bloque = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        bloque.append((record.id, str(record.seq).upper()))

        if len(bloque) >= block_size:
            bloque_id += 1
            print(f"🟦 Block {bloque_id} ({len(bloque)} seqs)")

            with open(output_path, "a") as fout:
                for sid, seq in bloque:
                    sid2, tax = classify_seq_v2(sid, seq, model_v2)
                    fout.write(f"{sid2}\t{tax}\n")

            total += len(bloque)
            bloque.clear()
            gc.collect()

    if bloque:
        bloque_id += 1
        print(f"🟩 Final block ({len(bloque)} seqs)")

        with open(output_path, "a") as fout:
            for sid, seq in bloque:
                sid2, tax = classify_seq_v2(sid, seq, model_v2)
                fout.write(f"{sid2}\t{tax}\n")

        total += len(bloque)

    end = time.time()
    print(f"🏁 V2 complete: {total} sequences in {end - start:.2f}s")

    return output_path, None

# ============================================================
# 6. NAIVE BAYES V2 — Parallel (fork-safe)
# ============================================================

def classify_naive_v2_parallel(
    fasta_path,
    output_path,
    modelo_path="/data/databases/metagenapp_refs/16S/naive_model_v2.pkl",
    block_size=5000,
    n_threads=None
):
    import os
    import time
    import pickle
    import gc
    import multiprocessing as mp
    from Bio import SeqIO

    # decidir número de hilos si no se especifica
    if n_threads is None:
        n_threads = min(12, os.cpu_count() or 8)

    print(f"🧠 Naive V2 PARALLEL | threads={n_threads} | block={block_size}")

    # cargar modelo UNA sola vez
    with open(modelo_path, "rb") as f:
        model_v2 = pickle.load(f)

    # preparar carpeta de salida
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # limpiar archivo salida
    open(output_path, "w").close()

    # 🔥 usar fork (clave para compartir memoria en Linux)
    ctx = mp.get_context("fork")

    pool = ctx.Pool(
        processes=n_threads,
        initializer=_v2_init_worker,
        initargs=(model_v2,),
        maxtasksperchild=200
    )

    total = 0
    bloque_id = 0
    start = time.time()

    try:
        bloque = []

        for record in SeqIO.parse(fasta_path, "fasta"):
            bloque.append((record.id, str(record.seq).upper()))

            if len(bloque) >= block_size:
                bloque_id += 1
                print(f"🟦 Block {bloque_id} ({len(bloque)} seqs)")

                resultados = pool.imap_unordered(_v2_worker, bloque, chunksize=200)

                with open(output_path, "a") as fout:
                    for sid, tax in resultados:
                        fout.write(f"{sid}\t{tax}\n")

                total += len(bloque)
                bloque.clear()
                gc.collect()

        # procesar último bloque
        if bloque:
            bloque_id += 1
            print(f"🟩 Final block ({len(bloque)} seqs)")

            resultados = pool.imap_unordered(_v2_worker, bloque, chunksize=200)

            with open(output_path, "a") as fout:
                for sid, tax in resultados:
                    fout.write(f"{sid}\t{tax}\n")

            total += len(bloque)

    finally:
        pool.close()
        pool.join()

    end = time.time()
    print(f"🏁 V2 PARALLEL complete: {total} sequences in {end - start:.2f}s")

    return output_path, None

