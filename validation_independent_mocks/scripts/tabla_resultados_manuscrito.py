import pandas as pd, numpy as np
B="/data/validation_independent_mocks_20261111"
zymo=pd.read_csv(f"{B}/Zymo_D6306_16S_theoretical.tsv",sep="\t").groupby("genus_SILVA138.2")["16S_only_pct"].sum(); ref=100*zymo/zymo.sum()
man=pd.read_csv(f"{B}/manifest_download.tsv",sep="\t").drop_duplicates("run")
imap=pd.read_csv(f"{B}/input_map.tsv",sep="\t",header=None,names=["lib","run"]).merge(man[["run","read_count"]])
raw=pd.read_csv("/data/validation_independent_mocks_20261111/pares_reales_fastq.tsv",sep="\t",header=None,names=["lib","run","r1","r2"]).set_index("lib").r1
def load_mg(path):
    d=pd.read_csv(path,sep="\t"); cols=[c for c in d.columns if c.endswith("_L001")]
    d=d.rename(columns={c:"_".join(c.split("_")[:2]) for c in cols}); return d,["_".join(c.split("_")[:2]) for c in cols]
def load_q2():
    d=pd.read_csv(f"{B}/eval_qiime2_resumen.tsv",sep="\t"); return d,[c for c in d.columns if c[:2] in("A_","D_","E_")]
def load_mo():
    d=pd.read_csv(f"{B}/eval_mothur_resumen.tsv",sep="\t"); return d,[c for c in d.columns if c[:2] in("A_","D_","E_")]
def prof(d,cols,sal=False):
    tot=d[cols].sum(axis=1); g=tot.groupby(d.Genus.fillna("NA").astype(str).str.strip()).sum()
    if sal: g["Salmonella"]=g.get("Salmonella",0)+g.get("Enterobacteriaceae",0); g=g.drop("Enterobacteriaceae",errors="ignore")
    return (100*g/g.sum()) if g.sum()>0 else g, int(tot.sum())
def bc(p):
    fp=p[[x for x in p.index if x not in ref.index]].sum(); return (sum(abs(ref[x]-p.get(x,0)) for x in ref.index)+fp)/200, fp
tools=[("MetagenApp v1.1 (prereg.)",*load_mg(f"{B}/outputs/metagenapp_v11/asv_resumen.tsv"),False),
       ("mothur (prereg.)",*load_mo(),False),
       ("QIIME 2 principal (prereg.)",*load_q2(),False),
       ("QIIME 2 secundario Sal (prereg.)",*load_q2(),True),
       ("MetagenApp corregido (post hoc)",*load_mg(f"{B}/posthoc_primerfix/validacion_ADE/asv_resumen.tsv"),False)]
print("herramienta\tconj\tn\tpares_crudos\tlect_finales\tretencion_%\tret_min-max_%\tBC_pool\tBC_rep_media\tBC_rep_SD\tinesperado_%\tgeneros_todas")
for name,d,cols,sal in tools:
    for s in "ADE":
        cs=[c for c in cols if c.startswith(s+"_")]
        p,r=prof(d,cs,sal); b,fp=bc(p); reps=[bc(prof(d,[c],sal)[0])[0] for c in cs]
        rets=[100*prof(d,[c],sal)[1]/raw[c] for c in cs]; R=raw[cs].sum()
        allg=sum(all(prof(d,[c],sal)[0].get(x,0)>0 for c in cs) for x in ref.index)
        print(f"{name}\t{s}\t{len(cs)}\t{R}\t{r}\t{100*r/R:.1f}\t{min(rets):.1f}-{max(rets):.1f}\t{b:.3f}\t{np.mean(reps):.3f}\t{np.std(reps,ddof=1):.3f}\t{fp:.2f}\t{allg}/8")
    p,r=prof(d,cols,sal); b,fp=bc(p); print(f"{name}\tGLOBAL\t{len(cols)}\t{raw[cols].sum()}\t{r}\t{100*r/raw[cols].sum():.1f}\t\t{b:.3f}\t\t\t{fp:.2f}\t")
