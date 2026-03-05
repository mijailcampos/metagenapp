from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os


def filtrar_mafft_por_posicion(
    fasta_path,
    output_path,
    start,
    end,
    count_table_path=None,
    output_count_path=None,
    min_info_threshold=0.02,
    high_ambig_threshold=0.50
):
    """
    Advanced MAFFT alignment filtering.

    Steps:
        1. Trim alignment to [start:end] positions.
        2. Remove columns composed entirely of gaps.
        3. Remove low-information columns (< min_info_threshold real bases).
        4. Remove highly ambiguous columns (> high_ambig_threshold).
        5. Remove sequences that become fully gap after filtering.
        6. Optionally filter count_table.

    Returns
    -------
    tuple
        (output_path, retained_sequences, total_sequences, message)
    """

    # ------------------------------------------------------------
    # 0. Validaciones básicas
    # ------------------------------------------------------------
    if not os.path.exists(fasta_path):
        return None, 0, 0, f"Input FASTA not found: {fasta_path}"

    if start < 1 or end <= start:
        return None, 0, 0, "Invalid start/end positions."

    # ------------------------------------------------------------
    # 1. Leer FASTA alineado
    # ------------------------------------------------------------
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total = len(records)

    if total == 0:
        return None, 0, 0, "No sequences found in FASTA file."

    # ------------------------------------------------------------
    # 2. Construir matriz de alineamiento recortada
    # ------------------------------------------------------------
    try:
        matriz = np.array([
            list(str(rec.seq)[start - 1:end])
            for rec in records
            if len(rec.seq) >= end
        ])
    except Exception as e:
        return None, 0, total, f"Error building alignment matrix: {e}"

    if matriz.size == 0:
        return None, 0, total, "No sequences long enough after trimming."

    n_seqs, n_cols = matriz.shape
    cols_iniciales = n_cols

    # ------------------------------------------------------------
    # 3. Eliminar columnas 100% gap
    # ------------------------------------------------------------
    mask_no_full_gap = ~(np.all(matriz == "-", axis=0))
    matriz = matriz[:, mask_no_full_gap]

    # ------------------------------------------------------------
    # 4. Eliminar columnas con baja información
    # ------------------------------------------------------------
    bases_reales = np.isin(matriz, list("ACGTacgt"))
    fraccion_info = bases_reales.sum(axis=0) / n_seqs
    mask_info = fraccion_info >= min_info_threshold
    matriz = matriz[:, mask_info]

    # ------------------------------------------------------------
    # 5. Eliminar columnas altamente ambiguas
    # ------------------------------------------------------------
    ambiguos = np.isin(matriz, list("NRYWSKMBDHVnrywskmbdhv"))
    fraccion_ambig = ambiguos.sum(axis=0) / n_seqs

    # Aplicar solo si realmente hay problema de ambigüedad
    if (fraccion_ambig > high_ambig_threshold).any():
        mask_ambig = fraccion_ambig <= high_ambig_threshold
        matriz = matriz[:, mask_ambig]

    # ------------------------------------------------------------
    # 6. Reconstruir secuencias filtradas
    # ------------------------------------------------------------
    filtradas = []

    for i, rec in enumerate(records):
        if i >= matriz.shape[0]:
            continue

        new_seq = "".join(matriz[i])

        # Eliminar secuencias que quedaron solo en gaps
        if set(new_seq) == {"-"}:
            continue

        rec.seq = Seq(new_seq)
        filtradas.append(rec)

    retenidas = len(filtradas)

    if retenidas == 0:
        return None, 0, total, "All sequences removed after filtering."

    SeqIO.write(filtradas, output_path, "fasta")

    # ------------------------------------------------------------
    # 7. Filtrar count_table (opcional)
    # ------------------------------------------------------------
    if count_table_path and output_count_path:
        df = pd.read_csv(count_table_path, sep="\t")

        if "OTU" not in df.columns:
            return (
                output_path,
                retenidas,
                total,
                "Filtered FASTA created, but count_table lacks 'OTU' column."
            )

        ids_retenidos = {rec.id for rec in filtradas}
        df_filtrado = df[df["OTU"].isin(ids_retenidos)]
        df_filtrado.to_csv(output_count_path, sep="\t", index=False)

    # ------------------------------------------------------------
    # 8. Mensaje final
    # ------------------------------------------------------------
    msg = (
        f"MAFFT positional filtering completed. "
        f"Columns: {cols_iniciales} → {matriz.shape[1]}. "
        f"Sequences retained: {retenidas}/{total}."
    )

    return output_path, retenidas, total, msg
