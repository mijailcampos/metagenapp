# Manuscript figures

Scripts and source data of Figures 1-6 of the MetagenApp manuscript (ECCB 2026 Research Topic).

```bash
pip install matplotlib numpy pandas pillow
python figures_manuscript/make_figure1_frontiers.py   # Figure 1
python figures_manuscript/make_figures_frontiers.py   # Figures 2-6
```

TIFF and PNG files are written to `figures_manuscript/output/`.

| Source | Used for |
|---|---|
| `../results_v1.1/tables/`, `../results_v1.1/rarefaction/` | Figures 2-5 |
| `data/sample_site_carrier.csv` | Figure 4: site and carriage status of each library |
| `data/mock_profiles.json` | Figure 6: expected and observed genus profiles |
| `data/mock_accuracy_table.tsv` | Table 7 values, for reference |
