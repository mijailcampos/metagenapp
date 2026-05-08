# MetagenApp

**MetagenApp** is a 16S/18S metabarcoding pipeline that classifies microbial
communities from paired-end Illumina reads, producing OTU tables with taxonomic
annotation. It is designed to be fast, resource-efficient, and reproducible.

Developed by José Mijail Campos Compeán.

---

## Benchmark: MetagenApp vs Mothur vs QIIME 2

Dataset: pharyngeal microbiome, 36 samples, 1,519,852 input reads.

### Table 1 — Pipeline comparison

| Parameter | MetagenApp (SILVA v1) | Mothur | QIIME 2 |
|---|---|---|---|
| Clustering algorithm | VSEARCH 97% | OptiClust 97% | DADA2 (ASV) |
| Reads retained | 1,149,967 (75.7%) | 1,145,172 (75.4%) | 191,568 (12.6%) |
| OTUs / ASVs | 14,876 | 10,473 | 2,098 |
| Taxonomic DB | SILVA 138 NR99 (k-mer, 83K taxa) | SILVA 138 (Wang) | SILVA 138 (sklearn NB) |
| Genus resolution | **91.1%** | ~100% | ~98% |

### Table 2 — Computational resources

| Resource | MetagenApp (SILVA v1) | Mothur | QIIME 2 |
|---|---|---|---|
| Total time | **19 min 32s** | 7h 00m 05s | 28m 34s |
| Peak RAM | **5.3 GB** | 53.4 GB | 14.0 GB |
| CPUs | 16 | multicore | multicore |

In this benchmark, MetagenApp completed the analysis substantially faster (22x faster) and with lower memory consumption (10x less RAM) than mothur.

Benchmark data and figures are in [`results_v1/`](results_v1/).

---

## Features

- Custom classifier (`naive-v2`) — Wang-style bootstrap k-mer confidence scoring
- SILVA 138 NR99 model — 83,000 taxa, 91.1% genus resolution
- Alignment fallback — EDLib pairwise alignment for borderline sequences (94% threshold, Yarza et al. 2014)
- 26-step pipeline — merge → filter → cluster → classify → OTU table
- Multiple execution modes: `student`, `premium`, `turbo`, `ref`, `qa`
- Modern CLI with configurable profiles

---

## Installation

### Requirements

- Python >= 3.9
- The following external tools:

```
vsearch >= 2.22    https://github.com/torognes/vsearch
edlib              https://github.com/Martinsos/edlib  (ref mode)
mafft >= 7.5       https://mafft.cbrc.jp/alignment/software/  (optional)
```

On Linux/Mac, install via conda:
```bash
conda install -c bioconda vsearch mafft
pip install edlib
```

### Install

```bash
git clone https://github.com/mijailcampos/metagenapp
cd metagenapp
pip install -r requirements.txt
pip install -e .
```

### Download models

The SILVA v1 classifier (~2.3 GB) is hosted on Zenodo (DOI: [10.5281/zenodo.20076968](https://doi.org/10.5281/zenodo.20076968)):

```bash
python download_model.py
```

Models are saved to `~/.metagenapp/refs/` by default. To use a custom path:

```bash
python download_model.py --dest /path/to/refs
export METAGENAPP_REFS=/path/to/refs  # add to ~/.bashrc
```

### Platform support

| Platform | Support |
|---|---|
| Linux | Full |
| macOS | Full (install tools via conda/brew) |
| Windows | Requires WSL2 (vsearch has no native Windows binary) |

---

## Quick start

```bash
# Run on example data (3 samples, ~20k reads each)
metagenapp -i example_data/ -o results/ --threads 4 --mode ref

# Full run with a profile (auto-generates timestamped output directory)
metagenapp --profile ref --input /path/to/fastq/

# Resume from a specific step
metagenapp -i data/ -o results/ --from-step clasificacion_tax
```

---

## Pipeline overview

26 sequential steps, from paired FASTQ to OTU table with taxonomy:

```
FASTQ R1 + R2
    ↓ [01] Merge paired reads (VSEARCH mergepairs)
    ↓ [02] Quality filter (length, ambiguities)
    ↓ [03] Dereplicate + count per sample
    ↓ [04] Reference alignment (VSEARCH, 70% identity)
    ↓ [09] Trim region of interest
    ↓ [10] Cluster at 97% → centroids + UC file
    ↓ [11] Extract centroids (mode-dependent)
    ↓ [12–18] MAFFT alignment + chimera removal (student/premium/turbo)
    ↓ [20] Taxonomic classification  ←  core step
    ↓ [21] Remove non-target lineages (chloroplasts, mitochondria)
    ↓ [22–26] OTU table (per sample, merged)
```

### Execution modes

| Mode | Description |
|---|---|
| `ref` | Publication grade — EDLib alignment, no MAFFT |
| `premium` | Balanced — MAFFT alignment |
| `turbo` | Max speed — no alignment |
| `student` | Low resource — 10,000 centroid cap |
| `qa` | Quality control only |

---

## Classifier: naive-v2

Wang-style bootstrap confidence scoring over a k-mer inverted index.

1. Extract k-mers from query sequence
2. Discard k-mers with `psize >= 20` (too generic)
3. Sum k-mer weights per taxon at each taxonomic level
4. Bootstrap (100 iterations): resample k-mer pool → vote
5. Assign deepest level with confidence >= 0.60

**Model (v1.0):** SILVA 138 NR99, 83,000 taxa, confidence threshold 0.60.

---

## Repository structure

```
metagenapp/
├── metagenapp/           # CLI and pipeline orchestration
│   ├── cli/main.py       # Entry point (Typer)
│   ├── pipeline/         # 26-step pipeline modules
│   └── metagen_config.py # Global config and paths
├── metagenapp_core/      # Classification engines
│   ├── models/           # naive-v2, kraken-lite
│   └── utils/            # k-mer generation, LCA
├── scripts/              # Paper figures and model training
│   ├── Figure*.R         # Benchmark figures
│   ├── train_naive_v2.py # Train a new naive-v2 model
│   ├── build_silva_trainset.py
│   └── dev/              # Internal / experimental scripts
├── results_v1/           # Benchmark results (v1.0)
│   ├── tables/           # TSV tables (OTU, taxonomy, resources)
│   ├── figures/          # PDF/PNG figures
│   ├── logs/             # /usr/bin/time -v timing logs
│   ├── raw_outputs/      # Pipeline outputs (MetagenApp, Mothur, QIIME 2)
│   └── COMMIT_HASH.txt   # Exact commit and command used
├── example_data/         # 3 public 16S samples for testing (~20k reads)
├── requirements.txt
├── pyproject.toml
├── VERSION.txt
└── LICENSE
```

---

## Reproducing the benchmark

```bash
git clone https://github.com/mijailcampos/metagenapp
git checkout v1.0

Frozen benchmark release:
https://github.com/mijailcampos/metagenapp/releases/tag/v1.0

# See exact command and commit used:
cat results_v1/COMMIT_HASH.txt
```

Benchmark run commit: `8228551d` — tag `v1.0`

---

## Citation

> Campos Compeán, J.M. (2026). MetagenApp: a fast and resource-efficient 16S
> metabarcoding pipeline with SILVA-scale taxonomic resolution. *Manuscript in preparation.*

---

## License

See [LICENSE](LICENSE).
