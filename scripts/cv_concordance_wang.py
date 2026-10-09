# Concordancia por secuencia entre dos salidas formato mothur (motor Wang de MetagenApp vs mothur classify.seqs)
import re, sys
from collections import Counter

def load(path):
    tax = {}
    for l in open(path):
        sid, t = l.rstrip("\n").split("\t")
        tax[sid] = [(m.group(1), int(m.group(2))) for m in re.finditer(r"([^;]+?)\((\d+)\);", t)]
    return tax

def call(lin, conf):
    """Línea truncada al último rango con confianza >= conf (máx. 6 rangos = género)."""
    out = []
    for name, b in lin[:6]:
        if b < conf or name.endswith("_unclassified") or name == "unknown":
            break
        out.append(name)
    return tuple(out)

a, b = load(sys.argv[1]), load(sys.argv[2])
ids = sorted(set(a) & set(b))
print(f"secuencias comunes: {len(ids)} (A {len(a)}, B {len(b)})")
raw = sum(1 for i in ids if [n for n, _ in a[i][:6]] == [n for n, _ in b[i][:6]])
print(f"genero asignado SIN umbral igual: {100*raw/len(ids):.1f} %")
for conf in (60, 70, 80, 90):
    c = Counter()
    for i in ids:
        x, y = call(a[i], conf), call(b[i], conf)
        if x == y:
            c["igual"] += 1
        elif x[:min(len(x), len(y))] == y[:min(len(x), len(y))]:
            c["misma_rama_otra_profundidad"] += 1
        else:
            c["distinto"] += 1
    gen = sum(1 for i in ids if call(a[i], conf)[:6] == call(b[i], conf)[:6] and len(call(a[i], conf)) == 6)
    n = len(ids)
    print(f"conf {conf/100}: llamada identica {100*c['igual']/n:.1f} % | misma rama, otra profundidad {100*c['misma_rama_otra_profundidad']/n:.1f} % | "
          f"distinta {100*c['distinto']/n:.1f} % | ambos al mismo genero {100*gen/n:.1f} %")
# diferencias de confianza en el género asignado cuando ambos asignan el mismo género
d = [a[i][5][1] - b[i][5][1] for i in ids if len(a[i]) > 5 and len(b[i]) > 5 and a[i][5][0] == b[i][5][0]]
if d:
    d.sort()
    print(f"conf genero (A - B) cuando coinciden: n {len(d)}, mediana {d[len(d)//2]}, p5 {d[len(d)//20]}, p95 {d[len(d)*19//20]}")
