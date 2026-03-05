import os
import re
from Bio import SeqIO
import pandas as pd

# --------------------------------------------------
# Optional Streamlit support (GUI-safe / CLI-safe)
# --------------------------------------------------
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False

# Tracking del pipeline
from metagenapp.pipeline.step_tracker import start_step, end_step


# --------------------------------------------------
# Default paths (GUI)
# --------------------------------------------------
INPUT_FASTA = "user_data/outputs/assembled_contigs.fasta"
INPUT_COUNT = "user_data/outputs/assembled_contigs.count_table"

OUTPUT_FASTA = "user_data/outputs/filtered_contigs.fasta"
OUTPUT_COUNT = "user_data/outputs/filtered_contigs.count_table"


# ============================================================
#   summary_seqs (CORE, CLI + GUI)
# ============================================================
def summary_seqs(fasta_path):
    starts, ends, lengths, ambigs, polymers = [], [], [], [], []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).replace("-", "")
        length = len(seq)

        starts.append(1)
        ends.append(length)
        lengths.append(length)

        ambigs.append(sum(1 for b in seq if b in "RYWSKMBDHVrywskmbdhv"))

        homo = re.findall(r"(A+|C+|G+|T+)", seq, re.IGNORECASE)
        polymers.append(max((len(h) for h in homo), default=0))

    df = pd.DataFrame({
        "Start": starts,
        "End": ends,
        "NBases": lengths,
        "Ambigs": ambigs,
        "Polymer": polymers,
        "NumSeqs": [1] * len(lengths)
    })

    return pd.DataFrame({
        "Minimum": df.min(),
        "Median": df.median(),
        "Mean": df.mean()
    }).T.round(2)


# ============================================================
#   CORE FUNCTION — FILTER CONTIGS (CLI + GUI)
# ============================================================
def filter_contigs(
    fasta_path=INPUT_FASTA,
    count_table_path=INPUT_COUNT,
    output_fasta=OUTPUT_FASTA,
    output_count=OUTPUT_COUNT,
    max_ambigs=0,
    max_poly=8,
    min_length=250,
    max_length=600
):
    """
    Filter contigs by length, ambiguities and homopolymers.
    CLI- and GUI-safe core function.
    """

    kept_ids = []
    retained = 0

    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

    try:
        with open(output_fasta, "w") as out_f:
            for record in SeqIO.parse(fasta_path, "fasta"):

                seq = str(record.seq).replace("-", "")
                seq_len = len(seq)

                num_ambigs = sum(1 for b in seq if b in "RYWSKMBDHVrywskmbdhv")

                homo = re.search(
                    rf"(A{{{max_poly+1},}}|C{{{max_poly+1},}}|G{{{max_poly+1},}}|T{{{max_poly+1},}})",
                    seq,
                    re.IGNORECASE
                )

                if (
                    num_ambigs <= max_ambigs
                    and min_length <= seq_len <= max_length
                    and homo is None
                ):
                    SeqIO.write(record, out_f, "fasta")
                    kept_ids.append(record.id)
                    retained += 1

    except Exception as e:
        return None, None, 0, f"Error filtering contigs: {e}"

    # Filter count table if present
    if count_table_path and os.path.exists(count_table_path):
        try:
            df = pd.read_csv(count_table_path, sep="\t")
            df[df.iloc[:, 0].isin(kept_ids)].to_csv(
                output_count, sep="\t", index=False
            )
        except Exception:
            output_count = None
    else:
        output_count = None

    return output_fasta, output_count, retained, None


# ============================================================
#   STREAMLIT WRAPPER (GUI ONLY)
# ============================================================
def run_filter_contigs():
    """
    Streamlit UI wrapper.
    Not used by CLI.
    """

    if not STREAMLIT_AVAILABLE:
        raise RuntimeError("run_filter_contigs() requires Streamlit")

    st.subheader("Contig filtering")

    if not os.path.exists(INPUT_FASTA):
        st.error("assembled_contigs.fasta not found. Run the assembly step first.")
        return

    st.info(
        "Adjust filters or use automatic presets according to 16S / 18S target regions."
    )

    # ---------------------------------------------------------
    # Marker presets
    # ---------------------------------------------------------
    tipo = st.radio("Select the marker type:", ["16S", "18S"], horizontal=True)

    if tipo == "16S":
        default_min, default_max, default_amb, default_poly = 250, 500, 0, 8
    else:
        default_min, default_max, default_amb, default_poly = 300, 650, 0, 10

    st.markdown(f"""
    ### Preset values for **{tipo}**:
    - Minimum length: `{default_min}`
    - Maximum length: `{default_max}`
    - Maximum ambiguities: `{default_amb}`
    - Maximum homopolymer length: `{default_poly}`
    """)

    st.markdown("---")

    # ---------------------------------------------------------
    # Filtering parameters
    # ---------------------------------------------------------
    col1, col2 = st.columns(2)

    with col1:
        min_length = st.number_input(
            "Minimum length",
            min_value=100,
            max_value=2000,
            value=default_min,
            step=1
        )

    with col2:
        max_length = st.number_input(
            "Maximum length",
            min_value=min_length + 1,
            max_value=3000,
            value=default_max,
            step=1
        )

    max_ambigs = st.number_input(
        "Maximum ambiguities allowed",
        min_value=0,
        max_value=10,
        value=default_amb,
        step=1
    )

    max_poly = st.number_input(
        "Maximum homopolymer size",
        min_value=1,
        max_value=30,
        value=default_poly,
        step=1
    )

    st.markdown("---")

    # ---------------------------------------------------------
    # Run filtering
    # ---------------------------------------------------------
    if st.button("Apply filtering"):

        start_step("04_Filter_contigs")

        out_fasta, out_count, retained, error = filter_contigs(
            fasta_path=INPUT_FASTA,
            count_table_path=INPUT_COUNT,
            output_fasta=OUTPUT_FASTA,
            output_count=OUTPUT_COUNT,
            min_length=min_length,
            max_length=max_length,
            max_ambigs=max_ambigs,
            max_poly=max_poly
        )

        if error:
            end_step(success=False)
            st.error(error)
            return

        end_step(success=True)

        st.success(f"{retained} contigs were retained after filtering.")
        st.code(f"Filtered FASTA: {out_fasta}\nCount Table: {out_count}")

        # Summary AFTER filtering
        st.markdown("### Contig summary (after filtering)")
        try:
            summary_after = summary_seqs(out_fasta)
            st.dataframe(summary_after, use_container_width=True)
        except Exception:
            st.warning("Post-filter summary could not be generated.")
