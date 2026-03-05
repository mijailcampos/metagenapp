import os
import subprocess
import pandas as pd


def generar_otus_desde_asvs(
    fasta_asv,
    count_table_asv,
    outdir,
    identity=0.97,
    threads=8
):
    """
    Genera OTUs al 97% a partir de ASVs finales usando VSEARCH.
    Produce:
      - otu_97.uc
      - otu_97_centroids.fasta
      - otu_97.shared (estilo mothur)
      - otu_table_0_03.tsv
    """

    # ========================
    # Validaciones iniciales
    # ========================
    if not os.path.exists(fasta_asv):
        return None, f"No existe el archivo FASTA: {fasta_asv}"

    if not os.path.exists(count_table_asv):
        return None, f"No existe el archivo count_table: {count_table_asv}"

    # ========================
    # Preparar directorio salida
    # ========================
    os.makedirs(outdir, exist_ok=True)

    output_uc = os.path.join(outdir, "otu_97.uc")
    output_centroids = os.path.join(outdir, "otu_97_centroids.fasta")
    output_shared = os.path.join(outdir, "otu_97.shared")
    output_otu_table = os.path.join(outdir, "otu_table_0_03.tsv")

    # ========================
    # Ejecutar VSEARCH
    # ========================
    try:
        cmd = [
            "vsearch",
            "--cluster_fast", fasta_asv,
            "--id", str(identity),
            "--uc", output_uc,
            "--centroids", output_centroids,
            "--threads", str(threads)
        ]
        subprocess.run(cmd, check=True)

    except subprocess.CalledProcessError as e:
        return None, f"Error ejecutando VSEARCH: {e}"

    # ========================
    # Leer count_table ASV
    # ========================
    try:
        df = pd.read_csv(count_table_asv, sep="\t")
        df = df.rename(columns={df.columns[0]: "ASV"})
    except Exception as e:
        return None, f"Error leyendo count_table: {e}"

    # ========================
    # Leer archivo .uc
    # ========================
    asv_to_otu = {}

    try:
        with open(output_uc) as f:
            for line in f:
                parts = line.strip().split("\t")

                if line.startswith("S"):   # centroid
                    asv = parts[8]
                    otu = parts[8]
                    asv_to_otu[asv] = otu

                elif line.startswith("H"):  # hit
                    asv = parts[8]
                    otu = parts[9]
                    asv_to_otu[asv] = otu

    except Exception as e:
        return None, f"Error leyendo archivo .uc: {e}"

    # ========================
    # Asignar OTUs a ASVs
    # ========================
    df["OTU"] = df["ASV"].map(asv_to_otu)

    if df["OTU"].isna().sum() > 0:
        return None, "Hay ASVs sin OTU asignada"

    # ========================
    # Sumar abundancias por OTU
    # ========================
    samples = df.columns[1:-1]
    df_otu = df.groupby("OTU")[samples].sum().reset_index()

    # ========================
    # Guardar tabla OTU
    # ========================
    try:
        df_otu.to_csv(output_otu_table, sep="\t", index=False)
    except Exception as e:
        return None, f"Error guardando tabla OTU: {e}"

    # ========================
    # Crear archivo .shared (estilo mothur)
    # ========================
    try:
        with open(output_shared, "w") as out:
            out.write(
                "label\tGroup\tnumOtus\t" +
                "\t".join(df_otu["OTU"]) + "\n"
            )

            for sample in samples:
                out.write(f"0.03\t{sample}\t{df_otu.shape[0]}")
                for otu in df_otu["OTU"]:
                    value = df_otu.loc[df_otu["OTU"] == otu, sample].values[0]
                    out.write("\t" + str(value))
                out.write("\n")

    except Exception as e:
        return None, f"Error generando archivo .shared: {e}"

    return output_otu_table, None
