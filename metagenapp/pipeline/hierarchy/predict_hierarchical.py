import os
import joblib
from Bio import SeqIO
import time
from collections import defaultdict


def classify_hierarchical(
    fasta_path,
    output_path,
    phylum_model_path,
    genus_models_dir,
    block_size=2000
):

    print("🔬 Loading phylum model...")
    phylum_payload = joblib.load(phylum_model_path)
    phylum_clf = phylum_payload["model"]
    phylum_vectorizer = phylum_payload["vectorizer"]

    dirpath = os.path.dirname(output_path)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)

    open(output_path, "w").close()

    genus_cache = {}

    total = 0
    start = time.time()

    block_ids = []
    block_seqs = []

    def process_block(ids, seqs):

        results = {}

        # ---- PHYLUM prediction (batch)
        Xp = phylum_vectorizer.transform(seqs)
        phylum_preds = phylum_clf.predict(Xp)

        # ---- Agrupar por phylum
        phylum_groups = defaultdict(list)

        for i, phylum_pred in enumerate(phylum_preds):
            phylum_clean = phylum_pred.replace("phylum__", "")
            phylum_groups[phylum_clean].append(i)

        # ---- GENUS prediction por grupo
        for phylum_clean, indices in phylum_groups.items():

            model_path = os.path.join(
                genus_models_dir,
                f"model_genus_{phylum_clean}.joblib"
            )

            if phylum_clean not in genus_cache:
                if os.path.exists(model_path):
                    genus_cache[phylum_clean] = joblib.load(model_path)
                else:
                    genus_cache[phylum_clean] = None

            genus_payload = genus_cache[phylum_clean]

            if genus_payload is None:
                for idx in indices:
                    results[idx] = (
                        phylum_preds[idx],
                        "genus__unclassified"
                    )
                continue

            if "constant_genus" in genus_payload:
                for idx in indices:
                    results[idx] = (
                        phylum_preds[idx],
                        genus_payload["constant_genus"]
                    )
                continue

            genus_clf = genus_payload["model"]
            genus_vectorizer = genus_payload["vectorizer"]

            seq_subset = [seqs[i] for i in indices]
            Xg = genus_vectorizer.transform(seq_subset)
            genus_preds = genus_clf.predict(Xg)

            for j, idx in enumerate(indices):
                results[idx] = (
                    phylum_preds[idx],
                    genus_preds[j]
                )

        return results

    # ---- Leer FASTA en bloques
    for record in SeqIO.parse(fasta_path, "fasta"):

        block_ids.append(record.id)
        block_seqs.append(str(record.seq).upper())

        if len(block_ids) == block_size:

            results = process_block(block_ids, block_seqs)

            with open(output_path, "a") as fout:
                for i in range(len(block_ids)):
                    phylum_pred, genus_pred = results[i]
                    fout.write(
                        f"{block_ids[i]}\t{phylum_pred};{genus_pred}\n"
                    )

            total += len(block_ids)
            block_ids.clear()
            block_seqs.clear()

    # Procesar último bloque
    if block_ids:
        results = process_block(block_ids, block_seqs)

        with open(output_path, "a") as fout:
            for i in range(len(block_ids)):
                phylum_pred, genus_pred = results[i]
                fout.write(
                    f"{block_ids[i]}\t{phylum_pred};{genus_pred}\n"
                )

        total += len(block_ids)

    end = time.time()
    print(f"🏁 Hierarchical classification complete: {total} sequences in {end - start:.2f}s")

    return output_path
