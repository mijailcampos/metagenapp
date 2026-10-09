# Analiza los errores de género de la CV: tipo de etiqueta verdadera/predicha y efecto de un bootstrap estilo Wang (submuestra 1/8)
import sys, pickle, random
from collections import Counter, defaultdict
from multiprocessing import Pool
sys.path.insert(0, "/data/projects/paper_metagenapp_benchmark/metagenapp_source_classfix/scripts")
import cv_naive_v2 as cv

AMBIG = ("uncultured", "unclassified", "unidentified", "incertae", "metagenome", "unknown", "environmental")

def wang_variant(seq, psize_max, frac, seed):
    """Igual que votes_and_wins pero cada iteración del bootstrap usa una submuestra de tamaño n*frac."""
    m = cv.MODEL; random.seed(seed)
    from metagenapp_core.utils.kmers import generate_kmers
    pool, kid_data = [], {}
    for kmer, c in Counter(generate_kmers(seq, m["k"])).items():
        kid = m["kmer_to_id"].get(kmer)
        if kid is None: continue
        s, e = int(m["pt_indptr"][kid]), int(m["pt_indptr"][kid + 1]); ps = e - s
        if ps == 0 or ps >= psize_max: continue
        pool.extend([kid] * c)
        if kid not in kid_data:
            kid_data[kid] = [(1.0 / ps, cv._prefixes(m["taxa"][int(t)], 6)) for t in m["pt_data"][s:e]]
    if not pool: return {}, {}
    def vote(idxs):
        v = defaultdict(float)
        for i in idxs:
            for w, prefs in kid_data[pool[i]]:
                for p in prefs: v[p] += w
        return v
    n = len(pool); tv = vote(range(n))
    top = {}
    for lvl in range(1, 7):
        lv = {t: x for t, x in tv.items() if t.count(";") == lvl}
        if lv: top[lvl] = max(lv, key=lv.get)
    wins = defaultdict(int); b = max(1, int(n * frac))
    for _ in range(cv.N_BOOT):
        bs = vote(random.choices(range(n), k=b))
        for lvl, t in top.items():
            lv = {x: v for x, v in bs.items() if x.count(";") == lvl}
            if lv and max(lv, key=lv.get) == t: wins[lvl] += 1
    return top, dict(wins)

def work(a):
    sid, seq, frac = a
    return sid, frac, wang_variant(seq, 20, frac, hash(sid) % (2**31))

if __name__ == "__main__":
    cv.MODEL = pickle.load(open("/data/mock_community_20261001/cv/model_full_genus.pkl", "rb"))
    seqs = dict(cv.read_fasta(cv.CV / "test_v3v4.fasta"))
    truth = {l.split("\t")[0]: l.split("\t")[1].split(";") for l in list(open(cv.CV / "test_truth.tsv"))[1:]}
    with Pool(16) as p:
        res = p.map(work, [(s, q, f) for f in (1.0, 0.125) for s, q in seqs.items()], chunksize=20)
    for frac in (1.0, 0.125):
        print(f"\n===== bootstrap con submuestra {frac} de los k-mers (1.0 = motor actual; 0.125 = Wang) =====")
        for conf in (0.6, 0.7, 0.8, 0.9):
            c = Counter(); conf_pairs = Counter()
            for sid, f, (top, wins) in res:
                if f != frac: continue
                got = cv.assign(top, defaultdict(int, wins), conf); g = got.rstrip(";").split(";") if got else []
                t = truth[sid]; amb = any(a in t[5].lower() for a in AMBIG) if len(t) >= 6 else True
                if len(g) >= 6:
                    if g[:6] == t[:6]: c["ok"] += 1
                    else:
                        c["mal_ambig" if amb else "mal_real"] += 1
                        if not amb and conf == 0.8: conf_pairs[(t[5], g[5])] += 1
                else: c["arriba"] += 1
            n = sum(c.values())
            print(f"conf {conf}: genero_ok {100*c['ok']/n:.1f} | mal con verdad 'uncultured' etc. {100*c['mal_ambig']/n:.1f} | mal REAL {100*c['mal_real']/n:.1f} | arriba de genero {100*c['arriba']/n:.1f}")
            if conf == 0.8 and conf_pairs: print("   confusiones reales más comunes (verdad → predicho):", conf_pairs.most_common(8))
    print("\nverdades con género ambiguo en la muestra:", sum(1 for t in truth.values() if len(t) < 6 or any(a in t[5].lower() for a in AMBIG)), "de", len(truth))
