MODEL_NAME = "kraken_lite"
MODEL_TYPE = "taxonomy_classifier"

from multiprocessing import Pool
import re

from metagenapp_core.models.kraken_lite_engine import classify_seq_kraken
from metagenapp_core.models.kraken_minimizer_engine import classify_seq_minimizer

MODEL = None

DNA_RE = re.compile("[^ACGT]")


def init_worker(model):

    global MODEL
    MODEL = model


def _worker(args):

    sid, seq = args

    seq = DNA_RE.sub("", seq.upper())

    if "minimizer_index" in MODEL:
        return classify_seq_minimizer(sid, seq, MODEL)

    if "kmer_index" in MODEL:
        return classify_seq_kraken(sid, seq, MODEL)

    return sid, "Unclassified"


def read_fasta(path):

    sid = None
    seq = []
    seqs = []

    with open(path) as f:

        for line in f:

            line = line.strip()

            if line.startswith(">"):

                if sid:
                    seqs.append((sid, "".join(seq)))

                sid = line[1:]
                seq = []

            else:
                seq.append(line)

    if sid:
        seqs.append((sid, "".join(seq)))

    return seqs


def classify_kraken_parallel(
        fasta_path,
        output_path,
        model,
        threads=4):

    seqs = read_fasta(fasta_path)

    args = [(sid, seq) for sid, seq in seqs]

    with Pool(
        processes=threads,
        initializer=init_worker,
        initargs=(model,)
    ) as p:

        results = p.map(_worker, args, chunksize=200)

    with open(output_path, "w") as out:

        for sid, tax in results:
            out.write(f"{sid}\t{tax}\n")

    return output_path, None