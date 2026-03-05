import pandas as pd
from Bio import SeqIO

def build_final_count_table(
    source_count_table,
    final_fasta,
    output_count_table
):
    # leer IDs finales
    final_ids = {rec.id for rec in SeqIO.parse(final_fasta, "fasta")}

    # leer count table fuente
    df = pd.read_csv(source_count_table, sep="\t")

    id_col = df.columns[0]

    df_final = df[df[id_col].isin(final_ids)]

    df_final.to_csv(output_count_table, sep="\t", index=False)
