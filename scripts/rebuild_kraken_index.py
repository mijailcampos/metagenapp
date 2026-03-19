import pickle
from metagenapp_core.models.build_kraken_index import build_kraken_index

RAW_MODEL = "/data/databases/metagenapp_refs/16S/naive_model_v4.pkl"
OUTPUT = "/data/databases/metagenapp_refs/16S/kraken_index.pkl"

print("Loading raw model...")
with open(RAW_MODEL, "rb") as f:
    raw_model = pickle.load(f)

print("Building kraken-lite index...")
index = build_kraken_index(raw_model, max_taxa_per_kmer=100000)

print("Saving index...")
with open(OUTPUT, "wb") as f:
    pickle.dump(index, f)

print("Done.")
