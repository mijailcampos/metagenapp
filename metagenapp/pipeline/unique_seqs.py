# ==========================================
# pipeline/unique_seqs.py
# ==========================================

import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict


def unique_seqs(
    fasta_path,
    count_table_path=None,
    output_fasta="user_data/outputs/unique_post_mafft.fasta",
    output_count="user_data/outputs/unique_post_mafft.count_table"
):
    import os
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from collections import defaultdict

    seq_to_ids = defaultdict(list)

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        seq_to_ids[seq_str].append(record.id)

    original_count_df = None
    sample_names = []

    if count_table_path and os.path.exists(count_table_path):
        original_count_df = pd.read_csv(
            count_table_path, sep="\t", index_col=0
        )

        # 🔴 CLAVE: forzar a numérico
        original_count_df = original_count_df.apply(
            pd.to_numeric, errors="coerce"
        ).fillna(0).astype(int)

        sample_names = original_count_df.columns.tolist()

    unique_records = []
    new_count_data = {}

    # 🔒 orden estable y reproducible
    for i, (seq_str, original_ids) in enumerate(sorted(seq_to_ids.items())):
        new_id = f"unique_{i+1}"

        unique_records.append(
            SeqRecord(Seq(seq_str), id=new_id, description="")
        )

        if original_count_df is not None:
            summed = original_count_df.loc[
                original_count_df.index.intersection(original_ids)
            ].sum(axis=0)

            new_count_data[new_id] = summed.tolist()

    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    SeqIO.write(unique_records, output_fasta, "fasta")

    if original_count_df is not None:
        os.makedirs(os.path.dirname(output_count), exist_ok=True)

        df = pd.DataFrame.from_dict(
            new_count_data, orient="index", columns=sample_names
        )
        df.index.name = "Representative_Sequence"
        df.to_csv(output_count, sep="\t")

    return output_fasta, output_count, len(unique_records)
