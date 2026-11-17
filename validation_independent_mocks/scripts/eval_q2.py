import pandas as pd, numpy as np
B="/data/validation_independent_mocks_20261111"; O=f"{B}/outputs/qiime2"
zymo=pd.read_csv(f"{B}/Zymo_D6306_16S_theoretical.tsv",sep="\t").groupby("genus_SILVA138.2")["16S_only_pct"].sum(); ref=100*zymo/zymo.sum()
t=pd.read_csv(f"{O}/table_export/feature-table.tsv",sep="\t",skiprows=1,index_col=0)
tax=pd.read_csv(f"{O}/taxonomy_export/taxonomy.tsv",sep="\t",index_col=0)["Taxon"]
d=t.copy(); d["FullTaxonomy"]=tax.reindex(d.index); d["Genus"]=d.FullTaxonomy.str.split(";").str[-1].str.replace(r"^[a-z]__","",regex=True)
d.reset_index().rename(columns={"#OTU ID":"#OTU ID"}).to_csv(f"{B}/eval_qiime2_resumen.tsv",sep="\t",index=False)
def profile(cols):
    tot=d[cols].sum(axis=1); g=tot.groupby(d.Genus.fillna("NA").astype(str).str.strip()).sum(); return 100*g/g.sum(), int(g.sum()), int((tot>0).sum())
def bc(pct):
    fp=pct[[x for x in pct.index if x not in ref.index]].sum(); return (sum(abs(ref[x]-pct.get(x,0)) for x in ref.index)+fp)/200, fp
def bc_sal(pct):
    p=pct.copy(); p["Salmonella"]=p.get("Salmonella",0)+p.get("Enterobacteriaceae",0); p=p.drop("Enterobacteriaceae",errors="ignore"); return bc(p)[0]
samples=list(t.columns)
for c in samples:
    p,r,f=profile([c]); b,fp=bc(p); print(f"  {c:9s} lecturas {r:7d} ASVs {f:4d} generos {sum(p.get(x,0)>0 for x in ref.index)}/8  inesperado {fp:5.2f}%  BC {b:.3f}  BC_Sal {bc_sal(p):.3f}")
for s in ["A","D","E"]:
    cols=[c for c in samples if c.startswith(s+"_")]; p,r,f=profile(cols); b,fp=bc(p); reps=[bc(profile([c])[0])[0] for c in cols]
    print(f"Conjunto {s}: n={len(cols)} lecturas {r} ASVs {f} BC_pool {b:.3f} (Sal {bc_sal(p):.3f}) BC_rep {np.mean(reps):.3f}+-{np.std(reps,ddof=1):.3f} inesperado {fp:.2f}% generos_todas {sum(all(profile([c])[0].get(x,0)>0 for c in cols) for x in ref.index)}/8")
    print("   ",p.sort_values(ascending=False).head(11).round(2).to_dict())
p,r,f=profile(samples); b,fp=bc(p); print(f"GLOBAL lecturas {r} ASVs {f} BC_pool {b:.3f} (Sal {bc_sal(p):.3f}) inesperado {fp:.2f}%")
print("esperado",ref.round(1).to_dict())
