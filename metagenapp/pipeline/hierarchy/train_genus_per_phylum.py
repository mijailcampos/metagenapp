import argparse
import os
import re
import joblib
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.linear_model import SGDClassifier

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


def extract_phylum(lineage):
    for part in lineage.split(";"):
        part = part.strip()
        if part.startswith("p__"):
            return part
    return None


def extract_genus(lineage):
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
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--k", type=int, default=7)
    parser.add_argument("--n-features", type=int, default=131072)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--only", type=str, default=None)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading taxonomy...")
    tax_map = load_taxonomy_map(args.taxonomy)

    print("Grouping sequences by phylum...")
    phylum_groups = defaultdict(list)

    for sid, lineage in tax_map.items():
        phylum = extract_phylum(lineage)
        genus = extract_genus(lineage)

        if phylum:
            phylum_groups[phylum].append((sid, genus))

    print(f"Total phyla found: {len(phylum_groups)}")

    print("Indexing FASTA...")
    fasta_index = {sid: seq for sid, seq in read_fasta(args.fasta)}

    for phylum, entries in phylum_groups.items():

        phylum_clean = phylum.replace("p__", "")

        if args.only and phylum_clean != args.only:
            continue

        print(f"\n🧬 Training genus model for {phylum_clean}")
        print(f"Sequences: {len(entries)}")

        classes = sorted(set(genus for _, genus in entries))
        classes = np.array(classes)
        print(f"Genus classes: {len(classes)}")

        # Caso especial: solo una clase
        if len(classes) == 1:
            print("⚠ Only one genus. Saving constant classifier.")

            payload = {
                "constant_genus": classes[0],
                "phylum": phylum_clean
            }

            model_path = os.path.join(
                args.outdir,
                f"model_genus_{phylum_clean}.joblib"
            )

            joblib.dump(payload, model_path, compress=3)
            print(f"Saved constant model: {model_path}")
            continue

        # Modelo normal
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

        first_batch = True
        batch_sequences = []
        batch_labels = []
        total = 0

        for sid, genus in entries:
            if sid in fasta_index:
                batch_sequences.append(fasta_index[sid])
                batch_labels.append(genus)

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

        payload = {
            "model": clf,
            "vectorizer": vectorizer,
            "k": args.k,
            "n_features": args.n_features,
            "classes": classes,
            "phylum": phylum_clean
        }

        model_path = os.path.join(
            args.outdir,
            f"model_genus_{phylum_clean}.joblib"
        )

        joblib.dump(payload, model_path, compress=3)
        print(f"Saved: {model_path}")

    print("\n🏁 All genus models trained.")


if __name__ == "__main__":
    main()
