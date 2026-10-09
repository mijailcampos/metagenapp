"""
Clasificador bayesiano ingenuo de Wang et al. (2007), implementación propia de MetagenApp.

Wang Q, Garrity GM, Tiedje JM, Cole JR (2007). Naive Bayesian classifier for rapid
assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol 73:5261-5267.

Algoritmo:
- Palabras = 8-mers solapados (sin bases ambiguas), contados una vez por secuencia.
- Prior por palabra: P_i = (n(w_i) + 0.5) / (N + 1)   (N = secuencias de entrenamiento)
- Prob. condicional: P(w_i | G) = (m(w_i) + P_i) / (M + 1)   (M = secuencias del taxón G)
- Asignación: argmax_G sum_i log P(w_i | G) con todas las palabras de la consulta.
- Confianza: 100 bootstraps, cada uno con len(palabras) // 8 palabras tomadas con reemplazo;
  confianza de cada rango = fracción de bootstraps cuyo taxón coincide en ese rango.

Clases = línea taxonómica completa del archivo .tax (p. ej. hasta género o hasta especie).

Modelo disperso (formato "sparse-v1", 7 oct 2026). Como
  log P(w|G) = log P_w - log(M_G + 1) + log(1 + m(w,G) / P_w),
el puntaje de una clase es una constante común a todas las clases, menos |W| * log(M_G + 1),
más la suma de delta(w,G) = log(1 + m/P_w) sobre las palabras de la consulta que la clase sí
contiene. Solo se guardan los pares (palabra, clase) con m > 0, en formato CSR por palabra.
El argmax (y por lo tanto la asignación y las confianzas) es el mismo que con la matriz densa,
pero la memoria escala con los pares no nulos y no con clases x 4^8, lo que permite entrenar
hasta especie. Los modelos densos anteriores ("logp") se siguen pudiendo cargar.

Semilla: cada secuencia usa su propio generador, derivado de la semilla global y del contenido
de la secuencia, de modo que el resultado no depende del orden del archivo, del tamaño de lote
ni del número de procesos.

Uso:
  python wang_engine.py train    --fasta ref.fasta --tax ref.tax --output modelo.npz
  python wang_engine.py classify --model modelo.npz --fasta query.fasta --output salida.taxonomy
La salida de `classify` usa el formato de mothur: Nombre(confianza);Nombre(confianza);...
En el pipeline (--classifier wang) se usa classify_to_taxonomy().
"""
import argparse
import hashlib
import multiprocessing as mp
import os
import sys
import time
from collections import Counter

import numpy as np

K = 8
N_WORDS = 4 ** K
_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[ord(chr(_b).lower())] = _i
_CODE[ord("U")] = _CODE[ord("u")] = 3


def kmer_words(seq, unique=True):
    """Índices (0..4^8-1) de los 8-mers válidos de la secuencia."""
    a = _CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    if a.size < K:
        return np.empty(0, dtype=np.int64)
    win = np.lib.stride_tricks.sliding_window_view(a, K)
    valid = (win < 4).all(axis=1)
    codes = (win.astype(np.int64) * (4 ** np.arange(K - 1, -1, -1))).sum(axis=1)[valid]
    return np.unique(codes) if unique else codes


def read_fasta(path):
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line.replace("-", "").replace(".", ""))
    if name is not None:
        yield name, "".join(parts)


def read_tax(path):
    tax = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            sid, lineage = line.split("\t", 1)
            tax[sid] = ";".join(t.strip().replace(" ", "_") for t in lineage.strip().rstrip(";").split(";"))
    return tax


# ----------------------------------------------------------------------------- entrenamiento

def _count_chunk(args):
    """Pares (clase, palabra, m) con m = nº de secuencias de la clase que contienen la palabra."""
    class_ids, seqs = args
    cls_out, word_out, m_out = [], [], []
    i = 0
    while i < len(class_ids):
        j = i
        while j < len(class_ids) and class_ids[j] == class_ids[i]:
            j += 1
        w = np.concatenate([kmer_words(s) for s in seqs[i:j]]) if j > i else np.empty(0, np.int64)
        u, m = np.unique(w, return_counts=True)
        cls_out.append(np.full(u.size, class_ids[i], dtype=np.int32))
        word_out.append(u.astype(np.int32))
        m_out.append(m.astype(np.int32))
        i = j
    return np.concatenate(cls_out), np.concatenate(word_out), np.concatenate(m_out)


def train(fasta, taxfile, output, threads=16):
    t0 = time.time()
    tax = read_tax(taxfile)
    records = [(tax[n], s) for n, s in read_fasta(fasta) if n in tax]
    missing = sum(1 for _ in tax) - len(records)
    lineages = sorted({lin for lin, _ in records})
    idx = {lin: i for i, lin in enumerate(lineages)}
    records.sort(key=lambda r: idx[r[0]])
    class_of = np.array([idx[lin] for lin, _ in records], dtype=np.int64)
    G, N = len(lineages), len(records)
    print(f"[wang] {N} secuencias, {G} clases ({missing} IDs del .tax sin secuencia)", flush=True)

    # Trozos que no parten una clase entre dos procesos
    bounds = np.flatnonzero(np.diff(class_of)) + 1
    starts = np.concatenate([[0], bounds])
    target = max(1, N // (threads * 8))
    chunks, cur = [], 0
    for st in list(starts[1:]) + [N]:
        if st - cur >= target or st == N:
            chunks.append((cur, st))
            cur = st
    jobs = [(class_of[a:b].tolist(), [records[i][1] for i in range(a, b)]) for a, b in chunks]
    parts = []
    with mp.get_context("fork").Pool(threads) as pool:
        for res in pool.imap(_count_chunk, jobs):
            parts.append(res)
    del jobs, records
    cls = np.concatenate([p[0] for p in parts])
    word = np.concatenate([p[1] for p in parts])
    m = np.concatenate([p[2] for p in parts])
    del parts

    M = np.bincount(class_of, minlength=G).astype(np.float64)
    n_w = np.bincount(word, weights=m, minlength=N_WORDS)          # secuencias que contienen w
    prior = (n_w + 0.5) / (N + 1)
    order = np.argsort(word, kind="stable")                         # CSR por palabra
    word, cls, m = word[order], cls[order], m[order]
    del order
    indptr = np.zeros(N_WORDS + 1, dtype=np.int64)
    np.cumsum(np.bincount(word, minlength=N_WORDS), out=indptr[1:])
    delta = np.log1p(m / prior[word]).astype(np.float32)
    del word, m
    np.savez(output, format="sparse-v1", indptr=indptr, cls=cls, delta=delta,
             log_m1=np.log(M + 1.0), lineages=np.array(lineages, dtype=object),
             class_size=M.astype(np.int64), n_train=N, k=K)
    print(f"[wang] modelo disperso: {cls.size:,} pares palabra-clase "
          f"({(cls.nbytes + delta.nbytes) / 1e9:.2f} GB); guardado en {output} ({time.time() - t0:.0f} s)", flush=True)


# ----------------------------------------------------------------------------- clasificación

_MODEL = {}


def load_model(path):
    m = np.load(path, allow_pickle=True)
    lineages = [l.split(";") for l in m["lineages"]]
    depth = max(len(l) for l in lineages)
    # nombre de cada rango por clase; None si la línea es más corta
    ranks = [np.array([";".join(l[:r + 1]) if len(l) > r else None for l in lineages], dtype=object)
             for r in range(depth)]
    model = {"lineages": lineages, "ranks": ranks}
    if "logp" in m.files:                                           # modelo denso (v1.0)
        model["logp"] = m["logp"]
    else:
        model.update(indptr=m["indptr"], cls=m["cls"], delta=m["delta"], log_m1=m["log_m1"],
                     log_m1_f32=m["log_m1"].astype(np.float32))
    return model


def seq_rng(seq, seed):
    """Generador propio de cada secuencia: depende solo de la semilla y del contenido."""
    h = int.from_bytes(hashlib.blake2b(seq.upper().encode("ascii"), digest_size=8).digest(), "little")
    return np.random.default_rng([seed, h])


def _bootstrap_sel(words_n, iters, rng):
    n_pick = max(1, words_n // 8)
    pick = rng.integers(0, words_n, size=(iters, n_pick))
    sel = np.zeros((iters, words_n), dtype=np.float32)
    np.add.at(sel, (np.repeat(np.arange(iters), n_pick), pick.ravel()), 1.0)
    return sel, n_pick


def _scores_dense(words, model, sel):
    sub = model["logp"][words]                                      # (n_palabras x clases)
    return sub.sum(axis=0), sel @ sub


def _scores_sparse(words, model, sel):
    # Bloque denso solo con las palabras de la consulta (n_palabras x clases, float32 como el
    # motor denso) y multiplicación con BLAS; mismo argmax que la versión de matriz completa.
    ip, cls, delta = model["indptr"], model["cls"], model["delta"]
    a, b = ip[words], ip[words + 1]
    G = model["log_m1"].size
    block = np.zeros((words.size, G), dtype=np.float32)
    for i in range(words.size):
        block[i, cls[a[i]:b[i]]] = delta[a[i]:b[i]]
    full = block.sum(axis=0, dtype=np.float64) - words.size * model["log_m1"]
    boot = sel @ block
    boot -= sel.sum(axis=1, keepdims=True) * model["log_m1_f32"][None, :]
    return full, boot


def classify_one(name, seq, model, iters=100, rng=None):
    """Devuelve [(taxón, confianza 0-100), ...] para la línea asignada."""
    words = kmer_words(seq)
    if words.size == 0:
        return []
    sel, _ = _bootstrap_sel(words.size, iters, rng)
    full, boot = (_scores_dense if "logp" in model else _scores_sparse)(words, model, sel)
    best = int(np.argmax(full))
    winners = np.argmax(boot, axis=1)                               # clase ganadora de cada bootstrap
    out = []
    for r, taxon in enumerate(model["lineages"][best]):
        conf = int(round(100.0 * np.mean(model["ranks"][r][winners] == model["ranks"][r][best])))
        out.append((taxon, conf))
    return out


def _classify_batch(items):
    from threadpoolctl import threadpool_limits
    with threadpool_limits(1):                                      # un hilo de BLAS por proceso
        return [(n, classify_one(n, s, _MODEL, _MODEL["iters"], seq_rng(s, _MODEL["seed"]))) for n, s in items]


def _run(model_path, fasta, threads, iters, seed, batch=50):
    # El modelo se carga una vez en el proceso padre y los hijos lo comparten (fork, solo lectura).
    _MODEL.clear()
    _MODEL.update(load_model(model_path))
    _MODEL["iters"], _MODEL["seed"] = iters, seed
    seqs = list(read_fasta(fasta))
    batches = [seqs[i:i + batch] for i in range(0, len(seqs), batch)]
    with mp.get_context("fork").Pool(threads) as pool:
        for res in pool.imap(_classify_batch, batches):
            yield from res


def _mothur_line(lineage):
    return "".join(f"{t}({c});" for t, c in lineage) if lineage else "unknown(0);"


def classify(model_path, fasta, output, threads=16, iters=100, seed=42):
    """Salida con confianzas en todos los rangos (formato mothur, como classify.seqs cutoff=0)."""
    t0 = time.time()
    n = 0
    with open(output, "w") as out:
        for name, lineage in _run(model_path, fasta, threads, iters, seed):
            out.write(f"{name}\t{_mothur_line(lineage)}\n")
            n += 1
    print(f"[wang] {n} secuencias clasificadas en {time.time() - t0:.0f} s -> {output}", flush=True)


def classify_to_taxonomy(model_path, fasta, output, confidence=0.80, threads=16, iters=100, seed=42, max_rank=None):
    """
    Para el pipeline: escribe `output` con la línea truncada al último rango con
    confianza >= `confidence` (formato "Bacteria;Phylum;...;", como naive-v2), o
    "Unclassified" si ni el dominio la alcanza. Con `max_rank` (p. ej. 6 = género) la línea
    no baja de ese rango aunque el modelo tenga hojas más profundas (modelo hasta especie).
    Las confianzas completas, con todos los rangos, quedan en
    <output sin extensión>.wang_confidence.taxonomy.
    """
    t0 = time.time()
    cutoff = 100.0 * confidence
    conf_path = os.path.splitext(output)[0] + ".wang_confidence.taxonomy"
    depth = Counter()
    with open(output, "w") as out, open(conf_path, "w") as out_conf:
        for name, lineage in _run(model_path, fasta, threads, iters, seed):
            out_conf.write(f"{name}\t{_mothur_line(lineage)}\n")
            kept = []
            for taxon, c in lineage[:max_rank]:
                if c < cutoff:
                    break
                kept.append(taxon)
            depth[len(kept)] += 1
            out.write(f"{name}\t{''.join(t + ';' for t in kept) if kept else 'Unclassified'}\n")
    n = sum(depth.values())
    print(f"[wang] {n:,} sequences classified in {time.time() - t0:.0f} s (confidence {confidence:.2f})", flush=True)
    print("[wang] deepest rank reached: " + ", ".join(f"{d} ranks {depth[d]:,}" for d in sorted(depth)), flush=True)
    print(f"[wang] confidences -> {conf_path}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    t = sp.add_parser("train")
    t.add_argument("--fasta", required=True)
    t.add_argument("--tax", required=True)
    t.add_argument("--output", required=True)
    t.add_argument("--threads", type=int, default=16)
    c = sp.add_parser("classify")
    c.add_argument("--model", required=True)
    c.add_argument("--fasta", required=True)
    c.add_argument("--output", required=True)
    c.add_argument("--threads", type=int, default=16)
    c.add_argument("--iters", type=int, default=100)
    c.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    if a.cmd == "train":
        train(a.fasta, a.tax, a.output, a.threads)
    else:
        classify(a.model, a.fasta, a.output, a.threads, a.iters, a.seed)


if __name__ == "__main__":
    sys.exit(main())
