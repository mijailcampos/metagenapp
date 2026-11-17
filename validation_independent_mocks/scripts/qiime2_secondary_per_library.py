exec(open("tabla_resultados_manuscrito.py").read().split("print(\"herramienta")[0])
d,cols=load_q2()
for s in "ADE":
    cs=[c for c in cols if c.startswith(s+"_")]; r=[bc(prof(d,[c],True)[0])[0] for c in cs]
    print(s, " ".join(f"{c}:{v:.3f}" for c,v in zip(cs,r)), f"media {np.mean(r):.3f} ± {np.std(r,ddof=1):.3f}")
