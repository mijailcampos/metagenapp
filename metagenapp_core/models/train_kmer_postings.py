import re
from collections import defaultdict, Counter

DNA_RE = re.compile(r"[^ACGT]")  # quitamos N y todo lo raro para kmers exactos

def read_fasta(path):
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
    tax = {}
    with open(tax_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            tax[parts[0]] = parts[1]
    return tax

def split_lineage(lineage: str):
    parts = [p.strip().strip('"') for p in lineage.split(";") if p.strip()]
    return parts

def lca_lineage(lineages):
    if not lineages:
        return []
    m = min(len(x) for x in lineages)
    out = []
    for i in range(m):
        v = lineages[0][i]
        if all(x[i] == v for x in lineages[1:]):
            out.append(v)
        else:
            break
    return out

def train_raw_kmer_postings(
    fasta_path: str,
    tax_path: str,
    k: int = 12,
    max_taxa_per_kmer: int = 50,
    verbose_every: int = 5000,
):
    """
    Construye raw_model invertido:
      kmer_counts: dict[kmer] -> dict[taxid] -> count
      taxonomy_labels: list[taxid] -> lineage str
      kmer_size: k

    Control de RAM:
      - por cada kmer, limitamos a max_taxa_per_kmer taxids manteniendo los top counts.
    """
    tax_map = load_taxonomy_map(tax_path)

    lineage_to_taxid = {}
    taxonomy_labels = []  # taxid -> lineage

    kmer_counts = {}  # kmer -> dict[taxid] -> count

    seen = 0
    missing = 0

    for sid, seq in read_fasta(fasta_path):
        lineage = tax_map.get(sid)
        if lineage is None:
            missing += 1
            continue

        # map lineage -> taxid
        taxid = lineage_to_taxid.get(lineage)
        if taxid is None:
            taxid = len(taxonomy_labels)
            lineage_to_taxid[lineage] = taxid
            taxonomy_labels.append(lineage)

        # limpiar secuencia (quita N y raros)
        seq = DNA_RE.sub("", seq.upper())
        n = len(seq) - k + 1
        if n <= 0:
            continue

        local = Counter(seq[i:i+k] for i in range(n))

        for kmer, c in local.items():
            d = kmer_counts.get(kmer)
            if d is None:
                d = {taxid: c}
                kmer_counts[kmer] = d
                continue

            d[taxid] = d.get(taxid, 0) + c

            # cap por kmer para limitar RAM
            if len(d) > max_taxa_per_kmer:
                # deja solo los top max_taxa_per_kmer por count
                # (simple y efectivo)
                top = sorted(d.items(), key=lambda x: x[1], reverse=True)[:max_taxa_per_kmer]
                d.clear()
                d.update(top)

        seen += 1
        if verbose_every and seen % verbose_every == 0:
            print(f"seen={seen:,}  missing_tax={missing:,}  unique_kmers={len(kmer_counts):,}  taxa={len(taxonomy_labels):,}")

    raw_model = {
        "kmer_size": k,
        "taxonomy_labels": taxonomy_labels,
        "kmer_counts": kmer_counts,  # invertido
        "max_taxa_per_kmer": max_taxa_per_kmer,
        "seen": seen,
        "missing_tax": missing,
    }
    return raw_model

def build_kid_to_lca(raw_model):
    """
    Precalcula LCA (por lineage) para cada kmer del raw_model invertido.
    Retorna dict[kmer] -> lca_lineage_str
    """
    taxa = raw_model["taxonomy_labels"]
    kmer_counts = raw_model["kmer_counts"]

    kmer_to_lca = {}

    for kmer, counts in kmer_counts.items():
        taxids = list(counts.keys())
        if not taxids:
            continue
        lineages = [split_lineage(taxa[t]) for t in taxids]
        lca = lca_lineage(lineages)
        if not lca:
            continue
        kmer_to_lca[kmer] = ";".join(lca) + ";"

    return kmer_to_lca
