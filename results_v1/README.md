# MetagenApp v1.0 — Benchmark Results

**Paper:** MetagenApp: a fast, low-memory 16S metabarcoding pipeline  
**Code:** github.com/mijailcampos/metagenapp — git tag `v1.0` (commit 8228551)  
**Run date:** 2026-05-02  
**Author:** José Mijail Campos Compeán

---

## Dataset

- **Samples:** 36 human pharyngeal microbiome samples (Bustos et al.)
- **Marker:** 16S rRNA V3-V4 region
- **Total input reads:** 1,519,852 (paired-end Illumina)
- **Raw data:** `/data/raw_data/collaborators/samuel/faringe_nariz/`

---

## Pipelines compared

| Pipeline | Version | Classifier | Database |
|---|---|---|---|
| **MetagenApp** | v1.0 (tag v1.0) | naive-v2 Wang bootstrap + alignment fallback | SILVA 138 NR99 BacArc (83,759 taxa) |
| Mothur | 1.48 | Wang bootstrap | SILVA 138 NR99 |
| QIIME2 | 2023.9 | sklearn Naive Bayes | SILVA 138 NR99 (99% classifier) |

---

## Key results

| Metric | MetagenApp | Mothur | QIIME2 |
|---|---|---|---|
| Genus resolution | **91.1%** | ~100% | ~98% |
| Wall-clock time | **19m 32s** | 7h 00m 05s | 28m 34s |
| Max RAM | **5.3 GB** | 53.4 GB | 14.0 GB |

---

## Folder structure

```
results_v1/
├── README.md                        ← this file
├── tables/
│   ├── Table1_pipeline_comparison.tsv     ← Table 1 from paper (pipeline comparison)
│   ├── Table2_computational_resources.tsv ← Table 2 from paper (time & RAM)
│   ├── final_asv_table.tsv                ← OTU table (14,876 OTUs × 36 samples)
│   ├── summary_tax_phylum.tsv             ← MetagenApp phylum-level abundances
│   ├── asv_resumen.tsv                    ← Full classification summary per OTU
│   ├── TableS1_genus_resolution_breakdown.tsv  ← Supplementary Table S1
│   └── TableS1_genus_resolution_combined.tsv   ← Supplementary Table S1 (combined)
├── figures/
│   ├── Figure1_phylum_comparison_naivev2.{png,pdf}  ← Fig 1: phylum composition (3 pipelines)
│   ├── Figure2_genus_bacteroidetes.{png,pdf}         ← Fig 2: top genus within Bacteroidetes
│   ├── Figure_benchmark_phylum_stacked.{png,pdf}     ← Benchmark phylum stacked bar
│   ├── Figure_benchmark_genus_resolution.png         ← Genus resolution comparison
│   ├── Figure_benchmark_alpha_diversity.{png,pdf}    ← Alpha diversity (Shannon, Chao1)
│   └── Figure_benchmark_beta_diversity.png           ← Beta diversity (Bray-Curtis PCoA)
├── raw_outputs/
│   ├── metagenapp/
│   │   ├── final_asv.taxonomy          ← Taxonomy assignments for all OTUs
│   │   ├── final_asv_table.tsv         ← OTU count table
│   │   ├── final_clean.count_table     ← Mothur-format count table
│   │   ├── summary_tax_phylum.tsv      ← Per-sample phylum counts
│   │   └── centroids_vs_reference.edlib.tsv  ← Alignment hits (used by fallback)
│   ├── mothur/
│   │   ├── mothur_taxonomy.taxonomy    ← Mothur taxonomy assignments
│   │   └── mothur_phylum.tsv           ← Mothur phylum-level summary
│   └── qiime2/
│       ├── taxonomy_export/taxonomy.tsv     ← QIIME2 taxonomy assignments
│       ├── table_export/feature-table.biom  ← QIIME2 feature table (BIOM format)
│       ├── table_export/feature-table.tsv   ← QIIME2 feature table (TSV)
│       └── table.tsv                        ← QIIME2 raw OTU table
└── logs/
    ├── metagenapp_run.log   ← Full MetagenApp run log (timing, steps, fallback stats)
    └── mothur_run.log       ← Mothur run log
```

---

## How to reproduce the MetagenApp run

```bash
# Requires: MetagenApp v1.0 (git tag v1.0), SILVA model at /data/databases/metagenapp_refs/
/usr/bin/time -v metagenapp \
  -i /data/raw_data/collaborators/samuel/faringe_nariz \
  -o /data/results/runs/ref_runs/faringe_nariz_20260502_v1 \
  --mode ref --marker 16S --classifier naive-v2 --model-type silva \
  --min-length 430 --max-length 500 \
  --max-ambigs 0 --max-poly 8 -t 16
```

---

## Notes

- **Genus resolution** = fraction of OTUs classified to genus level or finer
- **Alignment fallback:** MetagenApp applies an edlib-based fallback (min_identity=0.94) for OTUs that don't reach genus level with the k-mer classifier. This upgraded 150 OTUs (+0.7 pp genus resolution vs. without fallback).
- **QIIME2 low retention (12.6%):** Expected — DADA2 applies strict quality filtering. At phylum level, composition agrees with other pipelines when normalized.
- **Pseudomonadota vs. Proteobacteria:** SILVA 138 uses modern nomenclature (GTDB-based). Proteobacteria is split into Pseudomonadota, Campylobacterota, and Bdellovibrionota. Figure 1 maps SILVA names to classical nomenclature for comparability.
