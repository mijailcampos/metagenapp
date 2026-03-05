import os
import pandas as pd

def export_to_phyloseq(
    count_table,
    taxonomy,
    outdir
):
    os.makedirs(outdir, exist_ok=True)

    # ---------- OTU TABLE ----------
    otu = pd.read_csv(count_table, sep="\t")
    otu = otu.rename(columns={otu.columns[0]: "ASV"})
    otu.to_csv(os.path.join(outdir, "phyloseq_otu.tsv"),
               sep="\t", index=False)

    # ---------- TAXONOMY ----------
    tax = pd.read_csv(taxonomy, sep="\t", header=None,
                      names=["ASV", "Taxonomy"])

    # limpiar comillas raras
    tax["Taxonomy"] = (
        tax["Taxonomy"]
        .str.replace('"', '', regex=False)
        .str.rstrip(";")
    )

    ranks = ["Kingdom", "Phylum", "Class",
             "Order", "Family", "Genus"]

    tax_split = tax["Taxonomy"].str.split(";", expand=True)
    tax_split.columns = ranks[:tax_split.shape[1]]

    tax_final = pd.concat([tax["ASV"], tax_split], axis=1)

    tax_final.to_csv(
        os.path.join(outdir, "phyloseq_tax.tsv"),
        sep="\t", index=False
    )

    # ---------- METADATA (mínima) ----------
    samples = otu.columns[1:]
    metadata = pd.DataFrame({
        "SampleID": samples,
        "Group": samples
    })

    metadata.to_csv(
        os.path.join(outdir, "phyloseq_metadata.tsv"),
        sep="\t", index=False
    )

    # ---------- R SCRIPT ----------
    r_script = f"""
library(phyloseq)

otu <- read.table("phyloseq_otu.tsv", header=TRUE, sep="\\t", row.names=1)
tax <- read.table("phyloseq_tax.tsv", header=TRUE, sep="\\t", row.names=1)
meta <- read.table("phyloseq_metadata.tsv", header=TRUE, sep="\\t", row.names=1)

OTU  <- otu_table(as.matrix(otu), taxa_are_rows=TRUE)
TAX  <- tax_table(as.matrix(tax))
META <- sample_data(meta)

ps <- phyloseq(OTU, TAX, META)

ps
"""

    with open(os.path.join(outdir, "import_phyloseq.R"), "w") as f:
        f.write(r_script)

    return outdir
