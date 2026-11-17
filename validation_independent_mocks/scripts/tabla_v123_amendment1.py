exec(open("tabla_resultados_manuscrito.py").read().split("tools=[")[0])
d,cols=load_mg(f"{B}/reeval_v123/validacion_ADE/asv_resumen.tsv")
for s in "ADE":
    cs=[c for c in cols if c.startswith(s+"_")]; tot=d[cs].sum(axis=1); p,r=prof(d,cs); b,fp=bc(p); reps=[bc(prof(d,[c])[0])[0] for c in cs]
    allg=sum(all(prof(d,[c])[0].get(x,0)>0 for c in cs) for x in ref.index); pool=sum(p.get(x,0)>0 for x in ref.index)
    big=(tot>=0.001*tot.sum()); print(f"{s}\t{r:,} ({100*r/raw[cs].sum():.1f}%)\t{int((tot>0).sum())}\t{pool}/8 ({allg}/8)\t{fp:.2f}\t{b:.3f}\t{np.mean(reps):.3f} ± {np.std(reps,ddof=1):.3f}\tOTUs>=0.1%: {big.sum()} = {100*tot[big].sum()/tot.sum():.1f}%\tEnt x")
