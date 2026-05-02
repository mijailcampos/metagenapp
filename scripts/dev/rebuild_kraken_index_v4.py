import pickle

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import (
    build_kraken_index,
    invert_kmer_counts_by_taxon
)

FASTA = "/data/projects/TrainsetBuilder/TrainsetBuilder_streamlit/outputs/trainsets/TSB_SILVA_16S_V4_515F_806R.fasta"
TAX   = "/data/projects/TrainsetBuilder/TrainsetBuilder_streamlit/outputs/trainsets/TSB_SILVA_16S_V4_515F_806R.tax"

OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index_v4.pkl"

print("🧱 Training raw kmer postings (V4)...")

raw_model = train_raw_kmer_postings(
    FASTA,
    TAX,
    k=13
)

print("🔄 Inverting postings...")

raw_model["kmer_counts"] = invert_kmer_counts_by_taxon(
    raw_model["kmer_counts"],
    raw_model["taxonomy_labels"]
)

print("🧠 Building kraken-lite index...")

index = build_kraken_index(
    raw_model,
    max_taxa_per_kmer=100000
)

print("💾 Saving index...")

with open(OUTPUT, "wb") as f:
    pickle.dump(index, f)

print("✅ Done.")
