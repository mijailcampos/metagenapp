import re
import os
from collections import Counter
from multiprocessing import Pool, cpu_count

DNA_RE = re.compile(r"[^ACGT]")  # quitamos N y todo lo raro para kmers exactos

# ── I/O ───────────────────────────────────────────────────────────────────────

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


# ── Helpers de linaje ─────────────────────────────────────────────────────────

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


# ── Worker paralelo ───────────────────────────────────────────────────────────

_ACGT = frozenset("ACGT")
_RNA_TO_DNA = str.maketrans("Uu", "Tt")


def _worker(args):
    """
    Procesa un chunk de (taxid, seq) y devuelve kmer_counts local.
    No aplica cap aquí — se aplica una sola vez tras el merge.

    NOTAS CRÍTICAS:
    - Las secuencias SILVA están en formato RNA (U en lugar de T).
      Se convierten U→T antes de la extracción de k-mers para que los
      k-mers del índice coincidan con las queries de secuenciación (DNA, T).
    - Usa sliding-window saltando k-mers con bases no-ACGT (IUPAC ambiguity codes:
      R, Y, S, W, M, K, etc.) en lugar de DNA_RE.sub + Counter. El motivo:
      DNA_RE.sub eliminaba U y otros chars desplazando posiciones → k-mers
      artificiales que no coinciden con ninguna posición real de la secuencia.
    """
    records, k = args
    kmer_counts = {}

    for taxid, seq in records:
        seq = seq.upper().translate(_RNA_TO_DNA)  # U → T (RNA a DNA)
        n = len(seq) - k + 1
        if n <= 0:
            continue

        local = Counter()
        for i in range(n):
            kmer = seq[i:i + k]
            if not _ACGT.issuperset(kmer):  # salta k-mers con bases ambiguas
                continue
            local[kmer] += 1

        for kmer, c in local.items():
            d = kmer_counts.get(kmer)
            if d is None:
                kmer_counts[kmer] = {taxid: c}
            else:
                d[taxid] = d.get(taxid, 0) + c

    return kmer_counts


def _merge_results(results, max_taxa_per_kmer):
    """
    Fusiona N dicts kmer_counts sumando counts.
    Aplica el cap una sola vez al final (más eficiente que aplicarlo N veces).
    """
    merged = {}
    for kmer_counts in results:
        for kmer, d in kmer_counts.items():
            m = merged.get(kmer)
            if m is None:
                merged[kmer] = dict(d)
            else:
                for taxid, c in d.items():
                    m[taxid] = m.get(taxid, 0) + c

    # Cap final: solo los top-N taxids por kmer
    if max_taxa_per_kmer:
        for kmer, d in merged.items():
            if len(d) > max_taxa_per_kmer:
                top = sorted(d.items(), key=lambda x: x[1], reverse=True)[:max_taxa_per_kmer]
                merged[kmer] = dict(top)

    return merged


# ── API pública ───────────────────────────────────────────────────────────────

def _extract_phylum(lineage: str) -> str:
    """Extrae el phylum de un linaje SILVA (partes[1] en Domain;Phylum;Class;Order;Family;Genus;...)."""
    parts = [p.strip().strip('"') for p in lineage.split(";") if p.strip()]
    return parts[1] if len(parts) >= 2 else parts[0] if parts else ""


def _extract_genus(lineage: str) -> str:
    """Extrae el genus de un linaje SILVA (partes[5] en Domain;Phylum;Class;Order;Family;Genus;...)."""
    parts = [p.strip().strip('"') for p in lineage.split(";") if p.strip()]
    if len(parts) >= 6:
        return parts[5]
    if len(parts) >= 5:
        return parts[4]   # familia como fallback si no hay genus
    return parts[-1] if parts else ""


def train_raw_kmer_postings(
    fasta_path: str,
    tax_path: str,
    k: int = 12,
    max_taxa_per_kmer: int = 50,
    max_seqs_per_genus: int = None,
    max_seqs_per_phylum: int = None,
    verbose_every: int = 5000,
    n_jobs: int = None,
):
    """
    Construye raw_model invertido en paralelo:
      kmer_counts: dict[kmer] -> dict[taxid] -> count
      taxonomy_labels: list[taxid] -> lineage str
      kmer_size: k

    max_seqs_per_genus: limita a N secuencias por genus.
      Previene que genera sobre-representados (ej. Streptococcus con 500+ entradas)
      dominen el índice frente a genera menos abundantes.

    max_seqs_per_phylum: limita a N secuencias por phylum (aplicado después del cap
      por genus). Previene que phyla con mayor diversidad de géneros (ej. Firmicutes
      con ~1000 géneros × 150 seqs >> Bacteroidetes con ~200 géneros × 150 seqs)
      dominen el trainset. Se aplica respetando primero el cap por genus.

    n_jobs: workers a usar (default: cpu_count() - 4, mínimo 1).
    """
    if n_jobs is None:
        n_jobs = max(1, cpu_count() - 4)

    # 1. Cargar taxonomía
    print(f"  Cargando taxonomía...")
    tax_map = load_taxonomy_map(tax_path)

    lineage_to_taxid = {}
    taxonomy_labels = []

    # 2. Leer FASTA completo y mapear a taxids
    print(f"  Leyendo FASTA y mapeando taxids...")
    records = []
    missing = 0
    genus_counts: Counter = Counter()
    phylum_counts: Counter = Counter()
    skipped_by_genus_cap = 0
    skipped_by_phylum_cap = 0

    for sid, seq in read_fasta(fasta_path):
        lineage = tax_map.get(sid)
        if lineage is None:
            missing += 1
            continue

        # Cap por genus: si ya tenemos max_seqs_per_genus secuencias de este
        # genus, descartar la actual para evitar que genera muy abundantes en
        # SILVA (Streptococcus, Lachnospiraceae, etc.) sesguen el índice.
        if max_seqs_per_genus is not None:
            genus = _extract_genus(lineage)
            if genus_counts[genus] >= max_seqs_per_genus:
                skipped_by_genus_cap += 1
                continue
            genus_counts[genus] += 1

        # Cap por phylum: aplicado después del genus cap para nivelar phyla con
        # distinta diversidad de géneros (Firmicutes >> Bacteroidetes en SILVA).
        if max_seqs_per_phylum is not None:
            phylum = _extract_phylum(lineage)
            if phylum_counts[phylum] >= max_seqs_per_phylum:
                skipped_by_phylum_cap += 1
                # Revertir conteo de genus si se aplicó
                if max_seqs_per_genus is not None:
                    genus_counts[_extract_genus(lineage)] -= 1
                continue
            phylum_counts[phylum] += 1

        taxid = lineage_to_taxid.get(lineage)
        if taxid is None:
            taxid = len(taxonomy_labels)
            lineage_to_taxid[lineage] = taxid
            taxonomy_labels.append(lineage)

        records.append((taxid, seq))

    seen = len(records)
    print(f"  Secuencias válidas : {seen:,}  missing_tax={missing:,}")
    if max_seqs_per_genus is not None:
        print(f"  Descartadas (cap genus={max_seqs_per_genus}): {skipped_by_genus_cap:,}")
        print(f"  Genera únicas      : {len(genus_counts):,}")
    if max_seqs_per_phylum is not None:
        print(f"  Descartadas (cap phylum={max_seqs_per_phylum}): {skipped_by_phylum_cap:,}")
        top5 = phylum_counts.most_common(5)
        print(f"  Top phyla en trainset: {top5}")
    print(f"  Taxa únicas        : {len(taxonomy_labels):,}")
    print(f"  Lanzando {n_jobs} workers en paralelo...")

    # 3. Dividir en chunks y procesar en paralelo
    chunk_size = max(1, (seen + n_jobs - 1) // n_jobs)
    chunks = [records[i:i+chunk_size] for i in range(0, seen, chunk_size)]
    args = [(chunk, k) for chunk in chunks]

    with Pool(n_jobs) as pool:
        results = pool.map(_worker, args)

    print(f"  Mergeando {len(results)} resultados parciales...")
    kmer_counts = _merge_results(results, max_taxa_per_kmer)

    raw_model = {
        "kmer_size": k,
        "taxonomy_labels": taxonomy_labels,
        "kmer_counts": kmer_counts,
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
