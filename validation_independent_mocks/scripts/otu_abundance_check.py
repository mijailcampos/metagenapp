import pandas as pd
B="/data/validation_independent_mocks_20261111"
exp={"Pseudomonas","Escherichia-Shigella","Salmonella","Limosilactobacillus","Enterococcus","Staphylococcus","Listeria","Bacillus"}
for lab,path in [("corr",f"{B}/posthoc_primerfix/validacion_ADE/asv_resumen.tsv"),("v11_paper_D6306","/data/mock_community_20261001/outputs/metagenapp_aachen_v11_mergetest/asv_resumen.tsv")]:
    try: d=pd.read_csv(path,sep="\t")
    except Exception as e: print(lab,e); continue
    cols=[c for c in d.columns if c.endswith("_L001")]
    sets={"A":"A_","D":"D_","E":"E_"} if lab=="corr" else {"all":""}
    for s,pre in sets.items():
        cs=[c for c in cols if c.startswith(pre)]; t=d[cs].sum(axis=1); t=t[t>0]; tot=t.sum(); g=d.loc[t.index,"Genus"].fillna("NA").astype(str).str.strip()
        print(f"{lab} {s}: OTUs {len(t)}; <=10 lect {(t<=10).sum()} ({100*t[t<=10].sum()/tot:.2f}% lect); <0.01% {(t<0.0001*tot).sum()}; >=0.1% {(t>=0.001*tot).sum()} ({100*t[t>=0.001*tot].sum()/tot:.1f}% lect); top8 {100*t.nlargest(8).sum()/tot:.1f}%")
        e=g.isin(exp); print(f"   OTUs de generos esperados {e.sum()} ({100*t[e].sum()/tot:.1f}% lect); OTUs inesperados {(~e).sum()}; inesperados >=0.1%: {((~e)&(t>=0.001*tot)).sum()}")
        print("   top inesperados:", (t[~e].groupby(g[~e]).sum()/tot*100).nlargest(6).round(2).to_dict())
