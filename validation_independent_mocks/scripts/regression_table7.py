import pandas as pd, numpy as np
M="/data/mock_community_20261001"; V="/data/validation_independent_mocks_20261111/posthoc_primerfix"
hmp=pd.read_csv(f"{M}/HMP_mockB_genus_16S_theoretical.tsv",sep="\t").set_index("genus_SILVA138.2")
z=pd.read_csv(f"{M}/Zymo_D6306_16S_theoretical.tsv",sep="\t").groupby("genus_SILVA138.2")["16S_only_pct"].sum()
refs={"HM-782D":100*hmp["HM-782D_16S_pct"]/hmp["HM-782D_16S_pct"].sum(),"HM-783D":100*hmp["HM-783D_16S_pct"]/hmp["HM-783D_16S_pct"].sum(),"D6306":100*z/z.sum()}
raw={"HM-782D":205194,"HM-783D":211769,"D6306":128657}
runs={"old":{"bar":f"{M}/outputs/metagenapp_barcelona_wangsp","aa":f"{M}/outputs/metagenapp_aachen_v11_mergetest"},
      "new":{"bar":f"{V}/regresion_barcelona","aa":f"{V}/regresion_aachen"}}
def prof(d,cs):
    t=d[cs].sum(axis=1); g=t.groupby(d.Genus.fillna("NA").astype(str).str.strip()).sum(); return 100*g/g.sum(), int(t.sum()), int((t>0).sum())
def bc(p,ref):
    fp=p[[x for x in p.index if x not in ref.index]].sum(); return (sum(abs(ref[x]-p.get(x,0)) for x in ref.index)+fp)/200, fp
P={}
for v,r in runs.items():
    for mock,key,pre in [("HM-782D","bar","HM782"),("HM-783D","bar","HM783"),("D6306","aa","")]:
        d=pd.read_csv(f"{r[key]}/asv_resumen.tsv",sep="\t"); cs=[c for c in d.columns if c.endswith("_L001") and c.startswith(pre)]
        ref=refs[mock]; ref=ref[ref>0]; p,n,f=prof(d,cs); b,fp=bc(p,ref); reps=[bc(prof(d,[c])[0],ref)[0] for c in cs]
        allg=sum(all(prof(d,[c])[0].get(x,0)>0 for c in cs) for x in ref.index); pool=sum(p.get(x,0)>0 for x in ref.index)
        P[(v,mock)]=p
        print(f"{v} {mock}: {n:,} ({100*n/raw[mock]:.1f}%) feat {f} {pool}/{len(ref)} ({allg}/{len(ref)}) inesp {fp:.2f} BC {b:.3f} rep {np.mean(reps):.3f} ± {np.std(reps,ddof=1):.3f}")
for mock in refs:
    a,b=P[("old",mock)],P[("new",mock)]; idx=a.index.union(b.index); dd=(a.reindex(idx,fill_value=0)-b.reindex(idx,fill_value=0)).abs()
    print(mock,"max dif pp por genero (pool):",round(dd.max(),4),dd.idxmax())
