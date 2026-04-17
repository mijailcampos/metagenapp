"""
MetaSpecies integration for MetagenApp pipeline.

Classifies representative sequences (centroids / ASVs) to species level
using the MetaSpecies SSI+ANI hybrid classifier, and writes a taxonomy
file compatible with MetagenApp's taxonomy pipeline.

Output format matches Mothur's classify.seqs output:
    seq_id<TAB>Domain;Phylum;Class;Order;Family;Genus;Species;
"""
from __future__ import annotations

import os
import sys
from typing import Optional

_METASPECIES_CODE = "/data/projects/metaspecies/code"
if _METASPECIES_CODE not in sys.path:
    sys.path.insert(0, _METASPECIES_CODE)


def classify_fasta_to_taxonomy(
    fasta_path: str,
    ref_dir: str,
    output_path: str,
    genome_map_path: Optional[str] = None,
    batch_mode: bool = False,
    min_vote_fraction: float = 0.5,
    top_k: int = 5,
    progress_callback=None,
) -> tuple[str, Optional[str]]:
    """
    Classifies all sequences in a FASTA file to species level and writes
    a taxonomy file compatible with MetagenApp's classification pipeline.

    Returns (output_path, None) on success, (None, error_message) on failure.
    """
    try:
        from metaspecies_core.classifier import MetaSpeciesClassifier
        from metaspecies_core.io import read_fasta
    except ImportError as e:
        return None, (
            f"MetaSpecies no encontrado en {_METASPECIES_CODE}. Error: {e}"
        )

    if not os.path.isdir(ref_dir):
        return None, f"Directorio de índice no encontrado: {ref_dir}"

    clf = MetaSpeciesClassifier.load(ref_dir)

    genome_map = None
    if genome_map_path and os.path.exists(genome_map_path):
        genome_map = _load_genome_map(genome_map_path)

    seqs = list(read_fasta(fasta_path))
    if not seqs:
        return None, f"No se encontraron secuencias en {fasta_path}"

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as out_f:
        if batch_mode:
            if progress_callback:
                progress_callback(0, 1, f"Clasificando {len(seqs)} secuencias (modo batch)...")

            result = clf.classify_batch(
                sequences=seqs,
                genome_map=genome_map,
                top_k=top_k,
                min_vote_fraction=min_vote_fraction,
            )
            lineage = _resolve_lineage(clf, result.consensus_species, result.vote_fraction, min_vote_fraction)

            for seq_id, _ in seqs:
                clean_id = seq_id.lstrip(">").split()[0]
                out_f.write(f"{clean_id}\t{lineage}\n")

            if progress_callback:
                progress_callback(1, 1, f"Consenso: {result.consensus_species} ({result.vote_fraction:.0%})")

        else:
            for i, (seq_id, seq) in enumerate(seqs):
                clean_id = seq_id.lstrip(">").split()[0]
                result = clf.classify(seq, query_id=clean_id, genome_map=genome_map, top_k=top_k)
                lineage = _resolve_lineage(clf, result.species, 1.0, 0.0)
                out_f.write(f"{clean_id}\t{lineage}\n")

                if progress_callback and i % 100 == 0:
                    progress_callback(i + 1, len(seqs), f"Clasificando... {i+1}/{len(seqs)}")

            if progress_callback:
                progress_callback(len(seqs), len(seqs), "Clasificación completada")

    return output_path, None


def _resolve_lineage(clf, species: str, vote_fraction: float, min_vote_fraction: float) -> str:
    if species in ("UNKNOWN", "AMBIGUOUS") or (min_vote_fraction > 0 and vote_fraction < min_vote_fraction):
        return "unclassified;"

    lineage = clf.get_lineage(species)
    if lineage:
        return lineage

    parts = species.split("_", 1)
    genus = parts[0] if parts else species
    return f"unknown_domain;unknown_phylum;unknown_class;unknown_order;unknown_family;{genus};{species};"


def _load_genome_map(path: str) -> dict:
    gmap = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                gmap[parts[0]] = parts[1]
    return gmap
