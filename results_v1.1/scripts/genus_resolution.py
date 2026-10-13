#!/usr/bin/env python3
"""Table 6 - genus-level taxonomic resolution by reads and by features.

Reconstruction (1 Oct 2026) of the ephemeral /tmp/genus_resolution.py used on
17 Sep for the manuscript. Criterion recorded at the time: the genus label (6th
rank) is "ambiguous/propagated" if it contains uncultured, unclassified,
unidentified, Incertae Sedis (with space or underscore), environmental, metagenome or unknown; it is
"unclassified" if absent or empty; otherwise "classified". Unique genera =
distinct classified genus names.

Usage: genus_resolution.py <metagenapp_run_dir> [out.tsv]
"""
import re
import sys
import pandas as pd

RAW = "data/samples36"
MOTHUR_TAX = ("data/samples36/mothur/stability.trim.contigs.good.unique.good.filter.unique."
              "precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy")
QIIME_TAX = f"{RAW}/qiime2/taxonomy.tsv"
QIIME_TABLE = f"{RAW}/qiime2/feature-table.tsv"
AMBIGUOUS = re.compile(
    r"uncultured|unclassified|unidentified|incertae[ _]sedis|environmental|metagenome|unknown",
    re.IGNORECASE,
)


def genus_of(lineage):
    ranks = [re.sub(r"\(\d+\)$", "", r.strip()) for r in str(lineage).strip().rstrip(";").split(";")]
    if len(ranks) < 6:
        return ""
    return re.sub(r"^[a-z]__", "", ranks[5]).strip()


def category(genus):
    if not genus:
        return "Unclassified"
    if AMBIGUOUS.search(genus):
        return "Ambiguous/propagated"
    return "Classified"


def summarize(name, genera, reads):
    df = pd.DataFrame({"genus": genera, "reads": reads})
    df["cat"] = df["genus"].map(category)
    out = {"pipeline": name, "features": len(df),
           "unique_genera": df.loc[df.cat == "Classified", "genus"].nunique()}
    for cat in ["Classified", "Ambiguous/propagated", "Unclassified"]:
        out[f"reads_{cat}"] = 100 * df.loc[df.cat == cat, "reads"].sum() / df["reads"].sum()
        out[f"features_{cat}"] = 100 * (df.cat == cat).mean()
    out["total_reads"] = int(df["reads"].sum())
    return out


def metagenapp(run_dir):
    tax = pd.read_csv(f"{run_dir}/final_asv.taxonomy", sep="\t")
    tab = pd.read_csv(f"{run_dir}/final_asv_table.tsv", sep="\t", index_col=0)
    reads = tab.sum(axis=1)
    tax = tax.set_index(tax.columns[0]).loc[tab.index]
    return summarize("MetagenApp", [genus_of(t) for t in tax.iloc[:, 0]], reads.values)


def mothur():
    tax = pd.read_csv(MOTHUR_TAX, sep="\t")
    return summarize("mothur", [genus_of(t) for t in tax["Taxonomy"]], tax["Size"].values)


def qiime2():
    tax = pd.read_csv(QIIME_TAX, sep="\t").set_index("Feature ID")
    tab = pd.read_csv(QIIME_TABLE, sep="\t", skiprows=1, index_col=0)
    reads = tab.sum(axis=1)
    return summarize("QIIME 2", [genus_of(t) for t in tax.loc[tab.index, "Taxon"]], reads.values)


if __name__ == "__main__":
    rows = [metagenapp(sys.argv[1]), mothur(), qiime2()]
    res = pd.DataFrame(rows).set_index("pipeline").T
    pd.set_option("display.float_format", lambda v: f"{v:,.1f}")
    print(res.to_string())
    if len(sys.argv) > 2:
        res.to_csv(sys.argv[2], sep="\t")
