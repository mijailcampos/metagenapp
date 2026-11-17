exec(open("tabla_resultados_manuscrito.py").read().split("print(\"herramienta")[0])
for name,d,cols,sal in tools:
    if name.startswith("QIIME 2 secundario"): continue
    out=[]
    for s in "ADE":
        cs=[c for c in cols if c.startswith(s+"_")]; p0,_=prof(d,cs,False); p1,_=prof(d,cs,True); e=p0.get("Enterobacteriaceae",0)
        out.append(f"{s}: {bc(p0)[0]:.6f} -> {bc(p1)[0]:.6f}  dif {bc(p1)[0]-bc(p0)[0]:+.2e} (Ent {e:.2f}%)")
    print(name[:32].ljust(32)," ; ".join(out))
