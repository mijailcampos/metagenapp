# ============================================================
# summarize_contigs.py — Biotech Professional Version
# ============================================================

import os
import pandas as pd
import re
from Bio import SeqIO
# Optional Streamlit support
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False


# Tracking del pipeline
from metagenapp.pipeline.step_tracker import start_step, end_step

# ============================================================
# RUTAS ESTÁNDAR (OBLIGATORIAS)
# ============================================================
OUTPUT_FASTA = "user_data/outputs/assembled_contigs.fasta"
OUTPUT_COUNT = "user_data/outputs/assembled_contigs.count_table"

SUMMARY_TABLE = "user_data/outputs/contigs_summary.tsv"
SUMMARY_IMAGE = "user_data/outputs/contigs_summary.png"  # (para uso futuro)

# ============================================================
# FUNCIÓN PRINCIPAL (LÓGICA)
# ============================================================
def summarize_contigs(fasta_path, count_table_path=None):
    """
    Genera un resumen estadístico de los contigs ensamblados.
    """

    count_data = {}

    # Leer tabla de conteos si existe
    if count_table_path and os.path.exists(count_table_path):
        try:
            count_df = pd.read_csv(count_table_path, sep="\t", index_col=0)
            count_data = count_df.sum(axis=1).to_dict()
        except Exception:
            count_data = {}

    stats = []

    # Procesar secuencias
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        nbases = len(seq)

        # Bases ambiguas
        ambigs = sum(1 for base in seq.upper() if base not in "ACGT")


        # Homopolímeros máximos
        max_poly = 0
        for nucleotide in "ACGT":
            for match in re.finditer(rf"({nucleotide}+)", seq, re.IGNORECASE):
                max_poly = max(max_poly, len(match.group(0)))

        count_value = count_data.get(record.id, 1)

        stats.append([1, nbases, nbases, ambigs, max_poly, count_value])

    # Si no hubo secuencias
    if not stats:
        return {
            "Minimum": {"Start": 0, "End": 0, "NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "Median": {"Start": 0, "End": 0, "NBases": 0, "Ambigs": 0, "Polymer": 0, "NumSeqs": 0},
            "# of contigs": 0,
            "total sequences": 0,
        }

    df = pd.DataFrame(
        stats,
        columns=["Start", "End", "NBases", "Ambigs", "Polymer", "NumSeqs"]
    )

    summary = {
        "Minimum": df.min().to_dict(),
        "Median": df.median().round().astype(int).to_dict(),
        "Mean": df.mean().round(2).to_dict(),
        "# of contigs": len(df),
        "total sequences": df["NumSeqs"].sum()
    }

    return summary


# ============================================================
# INTERFAZ STREAMLIT (BIOTECH PRO)
# ============================================================
def run_summarize_contigs():
    """Interfaz Streamlit para resumen de contigs ensamblados."""

    if not STREAMLIT_AVAILABLE:
        raise RuntimeError("Streamlit is not available in CLI mode.")

    st.markdown(
        """
        <h1 style='font-size: 28px; font-weight: 700;'>
            Contigs Summary
        </h1>
        """,
        unsafe_allow_html=True
    )

    st.info(
        "This module computes descriptive statistics for assembled contigs, "
        "including sequence length, ambiguous bases, homopolymer distribution, "
        "and total per-sample counts when available."
    )


    # Validación del FASTA de entrada
    if not os.path.exists(OUTPUT_FASTA):
        st.error("The file 'assembled_contigs.fasta' was not found.")
        st.info("Run the assembly step before generating this summary.")
        return

    # 🔹 Inicia tracking del paso 03
    start_step("03_Summarize_contigs")

    # Ejecutar análisis
    with st.spinner("Computing contig statistics..."):
        try:
            summary = summarize_contigs(
                fasta_path=OUTPUT_FASTA,
                count_table_path=OUTPUT_COUNT if os.path.exists(OUTPUT_COUNT) else None
            )
            end_step(success=True)
        except Exception as e:
            end_step(success=False)
            st.error(f"Error computing summary: {e}")
            return

    st.success("Summary generated successfully.")

    # =========================================================
    # CONSTRUIR TABLA RESUMEN
    # =========================================================
    rows = []

    for metric, values in summary.items():
        if isinstance(values, dict):
            row = {"Metric": metric}
            row.update(values)
            rows.append(row)

    df_summary = pd.DataFrame(rows)

    st.subheader("Contig summary statistics")
    st.dataframe(df_summary, use_container_width=True)

    # Guardar tabla (CLAVE)
    df_summary.to_csv(SUMMARY_TABLE, sep="\t", index=False)

    # =========================================================
    # DESCARGA
    # =========================================================
    with open(SUMMARY_TABLE, "rb") as f:
        st.download_button(
            "Download contig summary table",
            f,
            file_name="contigs_summary.tsv",
            mime="text/tab-separated-values"
        )

    # =========================================================
    # ESTADÍSTICAS GLOBALES
    # =========================================================
    st.markdown(
        "<h3 style='font-size: 22px; margin-top: 20px;'>Global Statistics</h3>",
        unsafe_allow_html=True
    )

    st.write(f"Number of contigs: {summary['# of contigs']}")
    st.write(f"Total sequences: {summary['total sequences']}")
