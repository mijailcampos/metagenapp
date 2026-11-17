import pandas as pd, numpy as np, sys
B="/data/validation_independent_mocks_20261111"
zymo=pd.read_csv(f"{B}/Zymo_D6306_16S_theoretical.tsv",sep="\t").groupby("genus_SILVA138.2")["16S_only_pct"].sum(); ref=100*zymo/zymo.sum()
d=pd.read_csv(f"{B}/outputs/metagenapp_v11/asv_resumen.tsv",sep="\t")
def profile(cols):
    tot=d[cols].sum(axis=1); g=tot.groupby(d.Genus.fillna("NA").astype(str).str.strip()).sum(); return 100*g/g.sum(), int(g.sum()), int((tot>0).sum())
def bc(pct):
    fp=pct[[x for x in pct.index if x not in ref.index]].sum(); return (sum(abs(ref[x]-pct.get(x,0)) for x in ref.index)+fp)/200, fp
samples=[c for c in d.columns if c.endswith("_L001")]
for c in samples:
    p,r,f=profile([c]); b,fp=bc(p); print(f"  {c:18s} lecturas {r:7d} OTUs {f:4d} generos {sum(p.get(x,0)>0 for x in ref.index)}/8  inesperado {fp:5.2f}%  BC {b:.3f}")
for s in ["A","D","E"]:
    cols=[c for c in samples if c.startswith(s+"_")]; p,r,f=profile(cols); b,fp=bc(p); reps=[bc(profile([c])[0])[0] for c in cols]
    print(f"Conjunto {s}: n={len(cols)} BC_pool {b:.3f}  BC_rep {np.mean(reps):.3f}+-{np.std(reps,ddof=1):.3f}  inesperado {fp:.2f}%  generos_todas {sum(all(profile([c])[0].get(x,0)>0 for c in cols) for x in ref.index)}/8")
p,r,f=profile(samples); b,fp=bc(p); print(f"GLOBAL BC_pool {b:.3f} inesperado {fp:.2f}%")
for s in ["A","D","E"]:
    p,_,_=profile([c for c in samples if c.startswith(s+"_")]); print(s, p.sort_values(ascending=False).head(10).round(2).to_dict())
