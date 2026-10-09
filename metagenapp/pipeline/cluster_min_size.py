"""
OTU clustering with a minimum unique-sequence abundance, in the UPARSE order
(Edgar 2013; VSEARCH --minsize + --usearch_global, as in --otutab).

1. Only unique sequences with at least `min_size` reads in total (all samples)
   form OTUs: they are clustered with cluster_seqs (vsearch --cluster_fast).
2. Every remaining unique sequence (e.g. singletons) is then mapped to those
   OTU centroids with vsearch --usearch_global at the same identity. If it
   matches, its reads are added to that OTU; if not, its reads are discarded.

So a rare sequence cannot create an OTU of its own, but its reads still count
when it belongs to an existing OTU. The UC file written to `output_uc` holds the
cluster records plus the H records of the mapping, so aggregate_cluster_counts()
works unchanged.
"""

import os
import subprocess
from pathlib import Path

import pandas as pd

from metagenapp.pipeline.cluster import cluster_seqs


def _read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def cluster_with_min_size(fasta_path, count_table_path, output_fasta, output_uc,
                          identity=0.97, threads=None, min_size=2):
    """
    Returns (output_fasta, stats) on success, (None, error_message) on failure.
    stats: uniques/reads that formed OTUs, mapped back, and discarded.
    """
    fasta_path, output_uc = Path(fasta_path), Path(output_uc)
    if threads is None:
        threads = os.cpu_count()

    totals = pd.read_csv(count_table_path, sep="\t", index_col=0).sum(axis=1)
    abundant_fa = fasta_path.with_name(fasta_path.stem + f"_minsize{min_size}" + fasta_path.suffix)
    rare_fa = fasta_path.with_name(fasta_path.stem + f"_below_minsize{min_size}" + fasta_path.suffix)
    n_ab, rare_names = 0, []
    with open(abundant_fa, "w") as fa, open(rare_fa, "w") as fr:
        for name, seq in _read_fasta(fasta_path):
            if totals.get(name, 0) >= min_size:
                fa.write(f">{name}\n{seq}\n"); n_ab += 1
            else:
                fr.write(f">{name}\n{seq}\n"); rare_names.append(name)
    if n_ab == 0:
        return None, f"No unique sequence has >= {min_size} reads; nothing to cluster."

    # 1. OTUs from the abundant uniques only
    cluster_uc = output_uc.with_name(output_uc.stem + "_minsize" + output_uc.suffix)
    out, msg = cluster_seqs(fasta_path=abundant_fa, output_fasta=output_fasta,
                            output_uc=cluster_uc, identity=identity, threads=threads)
    if not out:
        return None, msg

    # 2. Map the rare uniques back to the OTU centroids
    map_uc = output_uc.with_name(output_uc.stem + "_below_minsize_map" + output_uc.suffix)
    hits = {}
    if rare_names:
        cmd = ["vsearch", "--usearch_global", str(rare_fa), "--db", str(output_fasta),
               "--id", str(identity), "--uc", str(map_uc), "--threads", str(threads)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            return None, f"VSEARCH usearch_global error:\n{res.stderr}"
        with open(map_uc) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] == "H":
                    hits[p[8]] = line

    with open(output_uc, "w") as out_uc:
        with open(cluster_uc) as fh:
            out_uc.writelines(fh)
        out_uc.writelines(hits.values())

    mapped = [n for n in rare_names if n in hits]
    stats = {
        "min_size": min_size,
        "uniques_clustered": n_ab,
        "reads_clustered": int(totals[totals >= min_size].sum()),
        "uniques_below": len(rare_names),
        "uniques_mapped": len(mapped),
        "reads_mapped": int(totals.reindex(mapped).sum()),
        "reads_discarded": int(totals.reindex(rare_names).sum() - totals.reindex(mapped).sum()),
    }
    return output_fasta, stats
