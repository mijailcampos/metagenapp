# Example Data

Three 16S rRNA paired-end samples (Illumina MiSeq, V3-V4 region) subsampled to
20,000 reads each. Suitable for testing MetagenApp on a standard laptop.

| Accession    | Title       | Source              | Reads  |
|-------------|-------------|---------------------|--------|
| SRR37813418 | THAU_16S_2  | Seawater, Thau Lagoon, France | 20,000 |
| SRR37813419 | THAU_16S_1  | Seawater, Thau Lagoon, France | 20,000 |
| SRR37813420 | SYLT_16S_2  | Seawater, Sylt, Germany       | 20,000 |

Study: SRP686975 — Oyster-associated metagenome (Magallana gigas, 2021)

Original reads downloaded from NCBI SRA. Subsampled with `head -80000` (20k reads × 4 lines).

## Quick start

```bash
metagenapp -i example_data/ -o results_example/ --threads 4 --mode ref
```
