import argparse
import joblib
import re
from Bio import SeqIO
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.linear_model import SGDClassifier
import numpy as np

DNA_RE = re.compile(r"[^ACGTN]")

def clean_seq(seq):
    return DNA_RE.sub("", seq.upper())

def read_fasta(path):
    for record in SeqIO.parse(path, "fasta"):
        yield record.id, clean_seq(str(record.seq))

def load_taxonomy_map(tax_path):
    tax = {}
    with open(tax_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                tax[parts[0]] = parts[1]
    return tax

def pick_genus(lineage):
    for part in lineage.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            val = part.split("g__")[1].strip()
            return f"genus__{val}" if val else "genus__unclassified"
    return "genus__unclassified"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--k", type=int, default=7)
    parser.add_argument("--n-features", type=int, default=131072)
    parser.add_argument("--batch-size", type=int, default=2000)
    args = parser.parse_args()

    print("Loading taxonomy...")
    tax_map = load_taxonomy_map(args.taxonomy)

    print("Collecting classes...")
    classes = set()
    for lineage in tax_map.values():
        classes.add(pick_genus(lineage))
    classes = np.array(sorted(classes))

    print(f"Total classes: {len(classes)}")

    vectorizer = HashingVectorizer(
        analyzer="char",
        ngram_range=(args.k, args.k),
        n_features=args.n_features,
        alternate_sign=False,
        lowercase=False
    )

    clf = SGDClassifier(
        loss="hinge",
        max_iter=5,
        tol=None
    )

    print("Training SGDClassifier incrementally...")
    first_batch = True
    batch_sequences = []
    batch_labels = []
    total = 0

    for sid, seq in read_fasta(args.fasta):
        if sid in tax_map:
            batch_sequences.append(seq)
            batch_labels.append(pick_genus(tax_map[sid]))

        if len(batch_sequences) == args.batch_size:
            X = vectorizer.transform(batch_sequences)
            y = np.array(batch_labels)

            if first_batch:
                clf.partial_fit(X, y, classes=classes)
                first_batch = False
            else:
                clf.partial_fit(X, y)

            total += len(batch_sequences)
            print(f"Processed {total}")

            batch_sequences = []
            batch_labels = []

    if batch_sequences:
        X = vectorizer.transform(batch_sequences)
        y = np.array(batch_labels)

        if first_batch:
            clf.partial_fit(X, y, classes=classes)
        else:
            clf.partial_fit(X, y)

        total += len(batch_sequences)

    payload = {
        "model": clf,
        "vectorizer": vectorizer,
        "k": args.k,
        "n_features": args.n_features,
        "classes": classes
    }

    joblib.dump(payload, args.out, compress=3)
    print("Model saved to:", args.out)

if __name__ == "__main__":
    main()
