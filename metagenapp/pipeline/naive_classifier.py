import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import pickle
import os
import gc


def entrenar_y_guardar_modelo(reference_fasta, taxonomy_path, k=10, modelo_path="modelos/naive_model.pkl"):

    # Read taxonomy file
    tax_df = pd.read_csv(
        taxonomy_path,
        sep="\t",
        header=None,
        names=["ID", "Taxonomy"]
    )

    tax_dict = dict(zip(tax_df["ID"], tax_df["Taxonomy"]))

    taxon_kmers = defaultdict(lambda: defaultdict(int))
    total_kmers_per_taxon = defaultdict(int)

    counter = 0

    # STREAMING FASTA (important)
    for record in SeqIO.parse(reference_fasta, "fasta"):

        taxon = tax_dict.get(record.id, "Unclassified")
        seq = str(record.seq).upper()

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            taxon_kmers[taxon][kmer] += 1
            total_kmers_per_taxon[taxon] += 1

        counter += 1

        if counter % 1000 == 0:
            gc.collect()
            print(f"Processed {counter} sequences")

    all_kmers = set(k for d in taxon_kmers.values() for k in d)

    cond_probs = defaultdict(dict)

    for taxon in taxon_kmers:
        total = total_kmers_per_taxon[taxon]

        for kmer in all_kmers:
            cond_probs[taxon][kmer] = (taxon_kmers[taxon].get(kmer, 0) + 1) / (total + len(all_kmers))

    os.makedirs(os.path.dirname(modelo_path), exist_ok=True)

    with open(modelo_path, "wb") as f:
        pickle.dump((cond_probs, list(taxon_kmers.keys()), k), f)

    print(f"Model saved at: {modelo_path}")