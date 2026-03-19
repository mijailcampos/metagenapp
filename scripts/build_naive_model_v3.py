import pickle
from metagenapp_core.models.train_kmer_postings import train_raw_kmer_postings, build_kid_to_lca
from metagenapp_core.models.model_adapters import adapt_naive_model

FASTA = "/data/databases/metagenapp_refs/trainset9/trainset9_032012.pds.fasta"
TAX   = "/data/databases/metagenapp_refs/trainset9/trainset9_032012.pds.tax"
OUT   = "/data/databases/metagenapp_refs/16S/naive_model_v4.pkl"

K = 15
MAX_TAXA_PER_KMER = 100000

print("🧱 Training raw postings...")
raw = train_raw_kmer_postings(
    fasta_path=FASTA,
    tax_path=TAX,
    k=K,
    max_taxa_per_kmer=MAX_TAXA_PER_KMER,
    verbose_every=5000
)

print("🧩 Adapting to v2 postings format...")
v3 = adapt_naive_model(raw)

print("🧠 Precomputing LCA per kmer...")
kmer_to_lca = build_kid_to_lca(raw)

kid_to_lca = [None] * len(v3["kmer_to_id"])
for kmer, kid in v3["kmer_to_id"].items():
    kid_to_lca[kid] = kmer_to_lca.get(kmer)

v3["kid_to_lca"] = kid_to_lca
v3["max_taxa_per_kmer"] = MAX_TAXA_PER_KMER

print(f"💾 Saving: {OUT}")
with open(OUT, "wb") as f:
    pickle.dump(v3, f, protocol=pickle.HIGHEST_PROTOCOL)

print("✅ Done")
print("k:", v3["k"], "kmers:", len(v3["kmer_to_id"]), "taxa:", len(v3["taxa"]))