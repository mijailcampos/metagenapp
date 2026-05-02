import pickle

from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings
from metagenapp_core.models.build_kraken_index import (
    build_kraken_index,
    invert_kmer_counts_by_taxon
)

FASTA = "/data/databases/metagenapp_refs/trainset9/trainset9_032012.pds.fasta"
TAX   = "/data/databases/metagenapp_refs/trainset9/trainset9_032012.pds.tax"

OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index.pkl"

print("🧱 Training raw kmer postings...")

raw_model = train_raw_kmer_postings(
    FASTA,
    TAX,
    k=15
)

print("🔄 Inverting kmer postings (kmer → taxon → kmer)...")

# convertir kmer → {taxid: count}
# a taxon → {kmer: weight}
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