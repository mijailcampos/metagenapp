#!/usr/bin/env python3
import argparse
import re
import joblib
from collections import defaultdict

from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.naive_bayes import MultinomialNB


DNA_RE = re.compile(r"[^ACGTN]")


def read_fasta(path):
    """Yield (id, seq) from fasta, streaming."""
    seq_id = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if seq_id is not None:
            yield seq_id, "".join(chunks)


def load_taxonomy_map(tax_path):
    """
    Expected taxonomy.tsv format:
    <id>\t<k__...;p__...;...;g__...;s__...>
    """
    tax = {}
    with open(tax_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            tax_id = parts[0]
            lineage = parts[1]
            tax[tax_id] = lineage
    return tax


def pick_label(lineage, rank):
    """
    rank in: kingdom/phylum/class/order/family/genus/species/full
    lineage like: k__Bacteria;p__...;...;g__X;s__Y
    """
    if rank == "full":
        return lineage

    want = {
        "kingdom": "k__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }[rank]

    items = lineage.split(";")
    value = ""
    for it in items:
        it = it.strip()
        if it.startswith(want):
            value = it[len(want):].strip()
            break

    if not value or value in {"", "unclassified", "unknown", "uncultured"}:
        return f"{rank}__unclassified"
    return f"{rank}__{value}"


def clean_seq(seq):
    seq = DNA_RE.sub("", seq.upper())
    # opcional: reemplazar N por vacío para ngramas exactos
    # seq = seq.replace("N", "")
    return seq


def iter_batches(fasta_path, tax_map, rank, batch_size):
    X_seqs, y = [], []
    missing = 0
    kept = 0

    for sid, seq in read_fasta(fasta_path):
        lineage = tax_map.get(sid)
        if lineage is None:
            missing += 1
            continue
        label = pick_label(lineage, rank)
        seq = clean_seq(seq)
        if len(seq) < 10:  # filtro mínimo absurdo para evitar basura
            continue

        X_seqs.append(seq)
        y.append(label)
        kept += 1

        if len(X_seqs) >= batch_size:
            yield X_seqs, y, kept, missing
            X_seqs, y = [], []

    if X_seqs:
        yield X_seqs, y, kept, missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--out", required=True)

    ap.add_argument("--k", type=int, default=8, help="k-mer size (ngram length)")
    ap.add_argument("--rank", default="genus",
                    choices=["kingdom", "phylum", "class", "order", "family", "genus", "species", "full"])
    ap.add_argument("--batch-size", type=int, default=2000)
    ap.add_argument("--n-features", type=int, default=2**18, help="Hashing space size (power of 2 recomendado)")
    ap.add_argument("--alpha", type=float, default=1.0, help="MultinomialNB alpha")
    ap.add_argument("--verbose-every", type=int, default=5, help="print cada N batches")
    args = ap.parse_args()

    print("Cargando taxonomía…")
    tax_map = load_taxonomy_map(args.taxonomy)
    print(f"Tax entries: {len(tax_map):,}")

    # 1) Pre-scan clases (sin cargar secuencias completas)
    print("Pre-scan de clases…")
    class_set = set()
    for _, lineage in tax_map.items():
        class_set.add(pick_label(lineage, args.rank))
    classes = sorted(class_set)
    print(f"Clases únicas ({args.rank}): {len(classes):,}")

    # 2) Vectorizador hashing (stateless, no crece)
    vectorizer = HashingVectorizer(
        analyzer="char",
        ngram_range=(args.k, args.k),
        n_features=args.n_features,
        alternate_sign=False,  # importante para NB
        lowercase=False
    )

    clf = MultinomialNB(alpha=args.alpha)

    # 3) Entrenamiento incremental
    print("Entrenando incremental…")
    batch_i = 0
    total_seen = 0
    last_missing = 0

    for seqs, y, kept, missing in iter_batches(args.fasta, tax_map, args.rank, args.batch_size):
        X = vectorizer.transform(seqs)

        if batch_i == 0:
            clf.partial_fit(X, y, classes=classes)
        else:
            clf.partial_fit(X, y)

        batch_i += 1
        total_seen = kept
        last_missing = missing

        if batch_i % args.verbose_every == 0:
            print(f"  batches={batch_i}  seen={total_seen:,}  missing_tax={last_missing:,}")

        # liberar referencias
        del X, seqs, y

    print(f"Listo. Total entrenado: {total_seen:,} | sin taxonomía: {last_missing:,}")

    payload = {
        "model": clf,
        "vectorizer": vectorizer,
        "rank": args.rank,
        "k": args.k,
        "n_features": args.n_features,
        "classes": classes,
    }

    joblib.dump(payload, args.out, compress=3)
    print(f"Modelo guardado en: {args.out}")


if __name__ == "__main__":
    main()
