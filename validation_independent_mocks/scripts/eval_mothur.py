import pandas as pd, numpy as np
B="/data/validation_independent_mocks_20261111"; O=f"{B}/outputs/mothur"
F=f"{O}/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc"
# make.file truncó 3 nombres de grupo; mapeo verificado contra stability.files
REN={"A":"A_ZY1","D":"D_mock4","E":"E_mock1"}
zymo=pd.read_csv(f"{B}/Zymo_D6306_16S_theoretical.tsv",sep="\t").groupby("genus_SILVA138.2")["16S_only_pct"].sum(); ref=100*zymo/zymo.sum()
sh=pd.read_csv(f"{F}.shared",sep="\t").drop(columns=["label","numOtus"]).set_index("Group").rename(index=REN).T
tax=pd.read_csv(f"{F}.0.03.cons.taxonomy",sep="\t",index_col=0)["Taxonomy"]
def genus(t):
    g=t.split(";")[5]; g=g.split("(")[0]; return g.replace("_unclassified","")
d=sh.copy(); d["FullTaxonomy"]=tax.reindex(d.index); d["Genus"]=d.FullTaxonomy.map(genus)
d.rename_axis("OTU").reset_index().to_csv(f"{B}/eval_mothur_resumen.tsv",sep="\t",index=False)
def profile(cols):
    tot=d[cols].sum(axis=1); g=tot.groupby(d.Genus.fillna("NA").astype(str).str.strip()).sum(); return 100*g/g.sum(), int(g.sum()), int((tot>0).sum())
def bc(pct):
    fp=pct[[x for x in pct.index if x not in ref.index]].sum(); return (sum(abs(ref[x]-pct.get(x,0)) for x in ref.index)+fp)/200, fp
def bc_sal(pct):
    p=pct.copy(); p["Salmonella"]=p.get("Salmonella",0)+p.get("Enterobacteriaceae",0); p=p.drop("Enterobacteriaceae",errors="ignore"); return bc(p)[0]
samples=sorted(sh.columns)
for c in samples:
    p,r,f=profile([c]); b,fp=bc(p); print(f"  {c:9s} lecturas {r:7d} OTUs {f:5d} generos {sum(p.get(x,0)>0 for x in ref.index)}/8  inesperado {fp:5.2f}%  BC {b:.3f}  BC_Sal {bc_sal(p):.3f}")
for s in ["A","D","E"]:
    cols=[c for c in samples if c.startswith(s+"_")]; p,r,f=profile(cols); b,fp=bc(p); reps=[bc(profile([c])[0])[0] for c in cols]
    print(f"Conjunto {s}: n={len(cols)} lecturas {r} OTUs {f} BC_pool {b:.3f} (Sal {bc_sal(p):.3f}) BC_rep {np.mean(reps):.3f}+-{np.std(reps,ddof=1):.3f} inesperado {fp:.2f}% generos_todas {sum(all(profile([c])[0].get(x,0)>0 for c in cols) for x in ref.index)}/8")
    print("   ",p.sort_values(ascending=False).head(11).round(2).to_dict())
p,r,f=profile(samples); b,fp=bc(p); print(f"GLOBAL lecturas {r} OTUs {f} BC_pool {b:.3f} (Sal {bc_sal(p):.3f}) inesperado {fp:.2f}%")
# OTUs con <=10 lecturas
tot=d[samples].sum(axis=1); print(f"OTUs<=10 lecturas: {(tot<=10).sum()} de {(tot>0).sum()}; OTUs>10 concentran {100*tot[tot>10].sum()/tot.sum():.1f}% lecturas")
