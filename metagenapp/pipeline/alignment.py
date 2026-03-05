# metagenapp/pipeline/alignment.py
from __future__ import annotations

import os
import math
import time
import tempfile
import subprocess
import multiprocessing as mp
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

def _worker_block_C(args):
    block_records, ref_items, index, k, max_candidates, mode = args

    import edlib
    import math

    out_rows = []

    for q_id, q_seq in block_records:
        qlen = len(q_seq)

        # candidate selection
        counts = {}
        for i in range(0, len(q_seq) - k + 1):
            kmer = q_seq[i:i+k]
            for ridx in index.get(kmer, []):
                counts[ridx] = counts.get(ridx, 0) + 1

        if not counts:
            out_rows.append((q_id, "", 0.0, -1, -1, -1, qlen))
            continue

        ranked = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        candidates = [ridx for ridx, _ in ranked[:max_candidates]]

        best_ref = ""
        best_ident = -1.0
        best_ed = 10**9
        best_start, best_end = -1, -1

        for ridx in candidates:
            r_id, r_seq = ref_items[ridx]
            res = edlib.align(q_seq, r_seq, mode=mode, task="locations")
            ed = res["editDistance"]
            locs = res.get("locations") or []

            if locs and locs[0][0] is not None:
                rs, re = locs[0]
                align_len = len(q_seq)
            else:
                rs, re = -1, -1
                align_len = len(q_seq)

            ident = max(0.0, 1.0 - (ed / align_len))

            if ident > best_ident or (math.isclose(ident, best_ident) and ed < best_ed):
                best_ident = ident
                best_ref = r_id
                best_ed = ed
                best_start, best_end = rs, re

        out_rows.append((q_id, best_ref, best_ident, best_ed, best_start, best_end, qlen))

    return out_rows

# -----------------------------
# FASTA utils (sin biopython)
# -----------------------------
def iter_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (id, sequence) from a FASTA. id is first token after '>'."""
    path = str(path)
    seq_id = None
    chunks: List[str] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks).upper()
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if seq_id is not None:
            yield seq_id, "".join(chunks).upper()


def read_fasta_dict(path: str | Path, max_seqs: Optional[int] = None) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for i, (sid, seq) in enumerate(iter_fasta(path)):
        d[sid] = seq
        if max_seqs is not None and (i + 1) >= max_seqs:
            break
    return d


def count_fasta_seqs(path: str | Path) -> int:
    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                c += 1
    return c


def write_fasta(records: Iterable[Tuple[str, str]], out_path: str | Path) -> None:
    out_path = str(out_path)
    with open(out_path, "w") as out:
        for sid, seq in records:
            out.write(f">{sid}\n")
            # wrap 80
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


# -----------------------------
# VSEARCH alignment (tu versión)
# -----------------------------
def align_with_vsearch(
    fasta_path: str | Path,
    reference_path: str | Path,
    output_path: str | Path,
    threads: Optional[int] = None,
    min_identity: float = 0.97,
) -> str:
    if threads is None:
        threads = os.cpu_count() or 1

    cmd = [
        "vsearch",
        "--usearch_global", str(fasta_path),
        "--db", str(reference_path),
        "--strand", "both",
        "--id", str(min_identity),
        "--threads", str(threads),
        "--alnout", str(output_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"VSEARCH failed:\n{result.stderr}")

    return str(output_path)


# ==========================================================
#  A) EDLib: mejor hit directo (simple, correcto, lento)
# ==========================================================
def _compute_identity_from_edlib(edit_distance: int, align_len: int) -> float:
    if align_len <= 0:
        return 0.0
    return max(0.0, 1.0 - (edit_distance / align_len))


def _edlib_best_hit_one_query(
    q_id: str,
    q_seq: str,
    ref_items: List[Tuple[str, str]],
    mode: str = "HW",          # semi-global: query dentro de ref
    task: str = "locations",   # necesitamos coords
) -> Tuple[str, Optional[str], float, int, int, int]:
    """
    Return: (q_id, best_ref_id, identity, editDistance, start, end)
    start/end are positions on reference (0-based) if available else -1.
    """
    import edlib

    best_ref = None
    best_ident = -1.0
    best_ed = 10**9
    best_start, best_end = -1, -1

    # Recorre todas las refs (A = brute force)
    for r_id, r_seq in ref_items:
        res = edlib.align(q_seq, r_seq, mode=mode, task=task)
        ed = res["editDistance"]
        locs = res.get("locations") or []
        # edlib devuelve (start, end) en el target (ref) para modo HW
        if locs and locs[0][0] is not None:
            start, end = locs[0]
            align_len = len(q_seq)  # en HW suele cubrir query
        else:
            start, end = -1, -1
            align_len = len(q_seq)

        ident = _compute_identity_from_edlib(ed, align_len)

        # criterio: mayor identidad, si empata menor edit distance
        if ident > best_ident or (math.isclose(ident, best_ident) and ed < best_ed):
            best_ident = ident
            best_ref = r_id
            best_ed = ed
            best_start, best_end = start, end

    return (q_id, best_ref, best_ident, best_ed, best_start, best_end)


def edlib_map_queries_to_reference_A(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    out_tsv: str | Path,
    max_refs: Optional[int] = None,
    mode: str = "HW",
) -> str:
    """
    A) Brute-force: para cada query, prueba contra TODAS las refs.
    Guarda el mejor hit en TSV.
    """
    ref_dict = read_fasta_dict(reference_fasta, max_seqs=max_refs)
    ref_items = list(ref_dict.items())

    t0 = time.time()
    n = 0
    with open(out_tsv, "w") as out:
        out.write("query_id\tref_id\tidentity\tedit_distance\tref_start\tref_end\tquery_len\n")
        for q_id, q_seq in iter_fasta(query_fasta):
            n += 1
            qlen = len(q_seq)
            q_id, best_ref, ident, ed, rs, re = _edlib_best_hit_one_query(q_id, q_seq, ref_items, mode=mode)
            out.write(f"{q_id}\t{best_ref or ''}\t{ident:.6f}\t{ed}\t{rs}\t{re}\t{qlen}\n")

    dt = time.time() - t0
    print(f"✅ EDLib A done: {n} queries mapped in {dt/60:.2f} min -> {out_tsv}")
    return str(out_tsv)


# ==========================================================
#  B) EDLib paralelo por bloques (misma lógica, más rápido)
# ==========================================================
def _worker_block(args):
    """
    Worker for multiprocessing.
    """
    block_records, ref_items, mode = args
    results = []
    for q_id, q_seq in block_records:
        qlen = len(q_seq)
        q_id, best_ref, ident, ed, rs, re = _edlib_best_hit_one_query(q_id, q_seq, ref_items, mode=mode)
        results.append((q_id, best_ref or "", ident, ed, rs, re, qlen))
    return results


def edlib_map_queries_to_reference_B(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    out_tsv: str | Path,
    threads: int = 4,
    block_size: int = 2000,
    max_refs: Optional[int] = None,
    mode: str = "HW",
) -> str:
    """
    B) Igual que A, pero paralelizado por bloques.
    """
    ref_dict = read_fasta_dict(reference_fasta, max_seqs=max_refs)
    ref_items = list(ref_dict.items())

    # Prepara bloques
    blocks = []
    current = []
    for q_id, q_seq in iter_fasta(query_fasta):
        current.append((q_id, q_seq))
        if len(current) >= block_size:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)

    t0 = time.time()
    print(f"🔧 EDLib B: {len(blocks)} blocks, threads={threads}, block_size={block_size}")

    with mp.Pool(processes=threads) as pool:
        it = pool.imap_unordered(_worker_block, [(blk, ref_items, mode) for blk in blocks])

        with open(out_tsv, "w") as out:
            out.write("query_id\tref_id\tidentity\tedit_distance\tref_start\tref_end\tquery_len\n")
            done = 0
            for chunk in it:
                for row in chunk:
                    q_id, ref_id, ident, ed, rs, re, qlen = row
                    out.write(f"{q_id}\t{ref_id}\t{ident:.6f}\t{ed}\t{rs}\t{re}\t{qlen}\n")
                done += 1
                if done % 10 == 0:
                    elapsed = (time.time() - t0) / 60
                    print(f"… blocks done {done}/{len(blocks)} ({elapsed:.1f} min)")

    dt = time.time() - t0
    print(f"✅ EDLib B done: mapped queries in {dt/60:.2f} min -> {out_tsv}")
    return str(out_tsv)


# ==========================================================
#  C) Índice k-mer para prefiltrar refs (plataforma)
# ==========================================================
def _kmers(seq: str, k: int) -> Iterator[str]:
    for i in range(0, len(seq) - k + 1):
        yield seq[i:i+k]


def build_kmer_index(
    reference_fasta: str | Path,
    k: int = 9,
    stride: int = 1,
    max_refs: Optional[int] = None,
) -> Tuple[List[Tuple[str, str]], Dict[str, List[int]]]:
    """
    index: kmer -> list of ref_indices
    """
    ref_dict = read_fasta_dict(reference_fasta, max_seqs=max_refs)
    ref_items = list(ref_dict.items())

    index: Dict[str, List[int]] = {}
    for ridx, (_, rseq) in enumerate(ref_items):
        # stride: para hacerlo más ligero
        for i in range(0, len(rseq) - k + 1, stride):
            kmer = rseq[i:i+k]
            index.setdefault(kmer, []).append(ridx)

    return ref_items, index


def _candidate_refs_for_query(
    q_seq: str,
    index: Dict[str, List[int]],
    k: int = 9,
    max_candidates: int = 50,
) -> List[int]:
    counts: Dict[int, int] = {}

    for kmer in _kmers(q_seq, k):
        for ridx in index.get(kmer, []):
            counts[ridx] = counts.get(ridx, 0) + 1

    if not counts:
        return []

    # top refs por número de k-mers compartidos
    ranked = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    return [ridx for ridx, _ in ranked[:max_candidates]]


def edlib_map_queries_to_reference_C(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    out_tsv: str | Path,
    threads: int = 4,
    k: int = 9,
    stride: int = 3,
    max_candidates: int = 50,
    mode: str = "HW",
    max_refs: Optional[int] = None,
) -> str:
    """
    C) Índice k-mer + edlib sólo contra candidatos (rápido y escalable).
    """
    ref_items, index = build_kmer_index(reference_fasta, k=k, stride=stride, max_refs=max_refs)
    print(f"🧠 K-mer index built: k={k}, stride={stride}, refs={len(ref_items)}, keys={len(index)}")


    # bloques
    blocks = []
    current = []
    for q_id, q_seq in iter_fasta(query_fasta):
        current.append((q_id, q_seq))
        if len(current) >= 2000:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)

    t0 = time.time()
    print(f"🚀 EDLib C: blocks={len(blocks)}, threads={threads}, max_candidates={max_candidates}")

    with mp.Pool(processes=threads) as pool:
        args_iter = [
            (blk, ref_items, index, k, max_candidates, mode)
            for blk in blocks
        ]

        it = pool.imap_unordered(_worker_block_C, args_iter)

        with open(out_tsv, "w") as out:
            out.write("query_id\tref_id\tidentity\tedit_distance\tref_start\tref_end\tquery_len\tk\tstride\tmax_candidates\n")
            done = 0
            for chunk in it:
                for (q_id, ref_id, ident, ed, rs, re, qlen) in chunk:
                    out.write(f"{q_id}\t{ref_id}\t{ident:.6f}\t{ed}\t{rs}\t{re}\t{qlen}\t{k}\t{stride}\t{max_candidates}\n")
                done += 1
                if done % 10 == 0:
                    elapsed = (time.time() - t0) / 60
                    print(f"… blocks done {done}/{len(blocks)} ({elapsed:.1f} min)")

    dt = time.time() - t0
    print(f"✅ EDLib C done in {dt/60:.2f} min -> {out_tsv}")
    return str(out_tsv)


# ----------------------------------------------------------
# Punto de integración en pipeline (lo que llamará MetagenApp)
# ----------------------------------------------------------
def align_centroids_refmode_professional(
    centroid_fasta: str | Path,
    reference_fasta: str | Path,
    out_tsv: str | Path,
    threads: int = 4,
    auto: bool = True,
) -> str:
    """
    Decide estrategia:
      - pequeño: B brute-force paralelo (más simple)
      - grande: C k-mer index + edlib
    """
    n = count_fasta_seqs(centroid_fasta)
    print(f"🔢 centroids={n}")

    if not auto:
        # fuerza C por default
        return edlib_map_queries_to_reference_C(
            centroid_fasta, reference_fasta, out_tsv, threads=threads
        )

    if n <= 5000:
        # brute-force pero paralelo
        return edlib_map_queries_to_reference_B(
            centroid_fasta, reference_fasta, out_tsv, threads=threads, block_size=500
        )
    else:
        # plataforma
        return edlib_map_queries_to_reference_C(
            centroid_fasta, reference_fasta, out_tsv,
            threads=threads,
            k=9,
            stride=3,
            max_candidates=50,
        )