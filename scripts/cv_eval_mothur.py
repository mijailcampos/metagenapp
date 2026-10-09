# Evalúa la salida de mothur classify.seqs (cutoff=0) sobre la muestra de CV con los mismos criterios que cv_naive_v2.py
import re, sys
from collections import Counter
CV = "/data/mock_community_20261001/cv"
AMBIG = ("uncultured", "unclassified", "unidentified", "incertae", "metagenome", "unknown", "environmental")
truth = {}
for l in list(open(f"{CV}/test_truth.tsv"))[1:]:
    sid, t, g = l.rstrip("\n").split("\t"); truth[sid] = t.replace(" ", "_").split(";")
tax = {}
for l in open(sys.argv[1]):
    sid, t = l.rstrip("\n").split("\t")
    tax[sid] = [(m.group(1), int(m.group(2))) for m in re.finditer(r"([^;]+?)\((\d+)\);", t)]
print("clasificadas:", len(tax), "de", len(truth))
for conf in (60, 70, 80, 90):
    c = Counter()
    for sid, t in truth.items():
        ranks = []
        for name, b in tax.get(sid, []):
            if b < conf or name.endswith("_unclassified") or name == "unknown": break
            ranks.append(name)
        amb = len(t) < 6 or any(a in t[5].lower() for a in AMBIG)
        if len(ranks) >= 6:
            c["ok" if ranks[:6] == t[:6] else ("mal_ambig" if amb else "mal_real")] += 1
        elif not ranks: c["sin"] += 1
        else: c["arriba"] += 1
    n = len(truth)
    print(f"conf {conf/100}: genero_ok {100*c['ok']/n:.1f} | mal con verdad 'uncultured' etc. {100*c['mal_ambig']/n:.1f} | mal REAL {100*c['mal_real']/n:.1f} | arriba de genero {100*c['arriba']/n:.1f} | sin clasificar {100*c['sin']/n:.1f}")
