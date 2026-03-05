import argparse
import joblib
import re
from Bio import SeqIO
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.svm import LinearSVC

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
        if part.strip().startswith("g__"):
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
    args = parser.parse_args()

    print("Loading taxonomy...")
    tax_map = load_taxonomy_map(args.taxonomy)

    print("Reading sequences...")
    sequences = []
    labels = []

    for sid, seq in read_fasta(args.fasta):
        if sid in tax_map:
            sequences.append(seq)
            labels.append(pick_genus(tax_map[sid]))

    print("Vectorizing...")
    vectorizer = HashingVectorizer(
        analyzer="char",
        ngram_range=(args.k, args.k),
        n_features=args.n_features,
        alternate_sign=False,
        lowercase=False
    )

    X = vectorizer.transform(sequences)

    print("Training LinearSVC...")
    clf = LinearSVC()
    clf.fit(X, labels)

    payload = {
        "model": clf,
        "vectorizer": vectorizer,
        "k": args.k,
        "n_features": args.n_features
    }

    joblib.dump(payload, args.out, compress=3)
    print("Model saved to:", args.out)

if __name__ == "__main__":
    main()
