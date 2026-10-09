#!/usr/bin/env python3
"""
Validación cruzada (hold-out 10 %) del clasificador naive-v2 de MetagenApp a nivel GÉNERO.

Objetivo: elegir psize_max y confidence SIN usar las comunidades mock (que quedan como prueba
independiente). Se aparta el 10 % de SILVA 138.2 NR99 BacArc, se entrena con el 90 % y se
clasifican fragmentos V3–V4 (341F/805R) de las secuencias apartadas.

Subcomandos:
  split   → cv/train_full.{fasta,tax}, cv/train_v3v4.{fasta,tax}, cv/test_v3v4.fasta, cv/test_truth.tsv
  eval    → clasifica la muestra de prueba con un modelo y un rango de psize_max/confidence
            y escribe cv/eval_<nombre>.tsv

La evaluación replica los pasos de metagenapp_core.models.naive_v2_engine.classify_seq_v2
(votos con k-mers de psize < psize_max, bootstrap de 100) pero guarda las victorias por nivel,
así cada umbral de confianza se evalúa sin repetir el bootstrap. Antes de evaluar se comprueba
que la réplica da el mismo resultado que el motor original.
"""
import argparse, pickle, random, sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/data/projects/paper_metagenapp_benchmark/metagenapp_source_classfix")
from metagenapp_core.models.naive_v2_engine import classify_seq_v2, _prefixes
from metagenapp_core.utils.kmers import generate_kmers

FULL = Path("/data/databases/metagenapp_refs/16S_silva1382_full")
V3V4 = Path("/data/databases/metagenapp_refs/16S_silva1382_v3v4")
CV = Path("/data/mock_community_20261001/cv")
SEED = 42
TEST_FRAC = 0.10
N_TEST = 2000
N_BOOT = 100


def read_fasta(path):
    sid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(chunks)
                sid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if sid is not None:
        yield sid, "".join(chunks)


def read_tax(path):
    return dict(l.rstrip("\n").split("\t") for l in open(path) if "\t" in l)


def cmd_split(_):
    CV.mkdir(parents=True, exist_ok=True)
    tax = read_tax(FULL / "silva1382_full_all.tax")
    ids = sorted(tax)
    rng = random.Random(SEED)
    test = set(rng.sample(ids, int(len(ids) * TEST_FRAC)))
    for src, name in [(FULL / "silva1382_full_all.fasta", "full"), (V3V4 / "silva1382_v3v4_all.fasta", "v3v4")]:
        n = 0
        with open(CV / f"train_{name}.fasta", "w") as fa, open(CV / f"train_{name}.tax", "w") as ft:
            for sid, seq in read_fasta(src):
                if sid in test:
                    continue
                fa.write(f">{sid}\n{seq}\n")
                ft.write(f"{sid}\t{tax[sid]}\n")
                n += 1
        print(f"train_{name}: {n:,} secuencias")
    test_frag = {sid: seq for sid, seq in read_fasta(V3V4 / "silva1382_v3v4_all.fasta") if sid in test}
    sample = rng.sample(sorted(test_frag), min(N_TEST, len(test_frag)))
    train_genera = {tuple(tax[s].rstrip(";").split(";")[:6]) for s in ids if s not in test and len(tax[s].rstrip(";").split(";")) >= 6}
    with open(CV / "test_v3v4.fasta", "w") as fa, open(CV / "test_truth.tsv", "w") as ft:
        ft.write("id\ttruth\tgenus_in_train\n")
        for sid in sample:
            lin = tax[sid].rstrip(";").split(";")
            fa.write(f">{sid}\n{test_frag[sid]}\n")
            ft.write(f"{sid}\t{';'.join(lin)}\t{int(len(lin) >= 6 and tuple(lin[:6]) in train_genera)}\n")
    print(f"apartadas: {len(test):,} | con fragmento V3–V4: {len(test_frag):,} | muestra de prueba: {len(sample):,}")


MODEL = None


def votes_and_wins(seq, psize_max, seed):
    """Réplica de classify_seq_v2 (rama Wang) que devuelve top por nivel y victorias por nivel."""
    m = MODEL
    k, kmer_to_id, taxa = m["k"], m["kmer_to_id"], m["taxa"]
    pt_data, pt_indptr = m["pt_data"], m["pt_indptr"]
    random.seed(seed)
    pool, kid_data = [], {}
    for kmer, c in Counter(generate_kmers(seq, k)).items():
        kid = kmer_to_id.get(kmer)
        if kid is None:
            continue
        s, e = int(pt_indptr[kid]), int(pt_indptr[kid + 1])
        psize = e - s
        if psize == 0 or psize >= psize_max:
            continue
        pool.extend([kid] * c)
        if kid not in kid_data:
            w = 1.0 / psize
            kid_data[kid] = [(w, _prefixes(taxa[int(t)], max_level=6)) for t in pt_data[s:e]]
    if not pool:
        return {}, {}, 0
    tv = defaultdict(float)
    for kid in pool:
        for w, prefs in kid_data[kid]:
            for p in prefs:
                tv[p] += w
    by_level = defaultdict(dict)
    for t, v in tv.items():
        lvl = t.count(";")
        if 1 <= lvl <= 6:
            by_level[lvl][t] = v
    top = {lvl: max(by_level[lvl], key=by_level[lvl].get) for lvl in range(1, 7) if by_level.get(lvl)}
    wins = defaultdict(int)
    n = len(pool)
    for _ in range(N_BOOT):
        bs = defaultdict(float)
        for idx in random.choices(range(n), k=n):
            for w, prefs in kid_data[pool[idx]]:
                for p in prefs:
                    bs[p] += w
        for lvl, t in top.items():
            lv = {x: v for x, v in bs.items() if x.count(";") == lvl}
            if lv and max(lv, key=lv.get) == t:
                wins[lvl] += 1
    return top, wins, len(set(pool))


def assign(top, wins, conf):
    chosen = None
    for lvl in range(1, 7):
        if lvl not in top or wins[lvl] / N_BOOT < conf:
            break
        chosen = top[lvl]
    return chosen


def _work(args):
    sid, seq, psize_max = args
    top, wins, ninf = votes_and_wins(seq, psize_max, seed=hash(sid) % (2**31))
    return sid, psize_max, top, dict(wins), ninf


def cmd_eval(a):
    global MODEL
    MODEL = pickle.load(open(a.model, "rb"))
    print(f"modelo {a.model}: {len(MODEL['taxa']):,} taxa, k={MODEL['k']}")
    seqs = dict(read_fasta(CV / "test_v3v4.fasta"))
    truth = {}
    for l in list(open(CV / "test_truth.tsv"))[1:]:
        sid, t, g = l.rstrip("\n").split("\t")
        truth[sid] = (t.split(";"), int(g))

    # Comprobación: la réplica coincide con el motor original
    for sid in list(seqs)[:20]:
        for pm in (20, 200):
            top, wins, _ = votes_and_wins(seqs[sid], pm, seed=7)
            mine = assign(top, wins, 0.8)
            random.seed(7)
            _, ref = classify_seq_v2(sid, seqs[sid], MODEL, psize_max=pm, confidence=0.8)
            assert mine == ref, (sid, pm, mine, ref)
    print("réplica == motor original (20 secuencias × 2 psize_max): OK")

    jobs = [(sid, seq, pm) for pm in a.psize for sid, seq in seqs.items()]
    with Pool(a.threads) as p:
        res = p.map(_work, jobs, chunksize=20)

    rows = []
    for pm in a.psize:
        for conf in a.conf:
            c = Counter()
            for sid, p_m, top, wins, ninf in res:
                if p_m != pm:
                    continue
                lin, in_train = truth[sid]
                got = assign(top, defaultdict(int, wins), conf)
                got_l = got.rstrip(";").split(";") if got else []
                grp = "conocido" if in_train else "novedoso"
                c[grp + "_n"] += 1
                if not got_l:
                    c[grp + "_sin_clasificar"] += 1
                elif len(got_l) >= 6:
                    c[grp + ("_genero_ok" if got_l[:6] == lin[:6] else "_genero_mal")] += 1
                else:
                    ok = got_l == lin[:len(got_l)]
                    c[grp + ("_arriba_ok" if ok else "_arriba_mal")] += 1
            nk, nn = c["conocido_n"], c["novedoso_n"]
            rows.append([pm, conf, nk,
                         100 * c["conocido_genero_ok"] / nk, 100 * c["conocido_genero_mal"] / nk,
                         100 * (c["conocido_arriba_ok"] + c["conocido_arriba_mal"]) / nk,
                         100 * c["conocido_sin_clasificar"] / nk, nn,
                         100 * c["novedoso_genero_mal"] / nn if nn else 0])
    hdr = ["psize_max", "conf", "n_conocidos", "genero_ok_%", "genero_mal_%", "arriba_de_genero_%",
           "sin_clasificar_%", "n_novedosos", "novedosos_con_genero_falso_%"]
    out = CV / f"eval_{a.name}.tsv"
    with open(out, "w") as fo:
        fo.write("\t".join(hdr) + "\n")
        for r in rows:
            fo.write("\t".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    print(open(out).read())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("split")
    e = sp.add_parser("eval")
    e.add_argument("--model", required=True)
    e.add_argument("--name", required=True)
    e.add_argument("--psize", type=int, nargs="+", default=[20, 50, 100, 200, 500])
    e.add_argument("--conf", type=float, nargs="+", default=[0.6, 0.7, 0.8])
    e.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    {"split": cmd_split, "eval": cmd_eval}[a.cmd](a)
