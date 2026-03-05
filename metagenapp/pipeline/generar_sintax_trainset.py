import pandas as pd
from Bio import SeqIO


def convertir_trainset_a_sintax(fasta_in, tax_in, fasta_out):
    """
    Convierte el trainset de Mothur (fasta + tax) en un archivo SINTAX
    compatible con VSEARCH.

    Salida:
        >seqID;tax=k:...,p:...,c:...,o:...,f:...,g:...,s:...
        ATGCATGCATGC...
    """

    # 1. Cargar taxonomía ".tax" de Mothur
    tax_dict = {}
    with open(tax_in) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue

            seq_id, taxonomy = parts
            niveles = taxonomy.split(";")

            sintax_string = []
            etiquetas = ["k", "p", "c", "o", "f", "g", "s"]

            for i, nivel in enumerate(niveles):
                if nivel.strip() == "":
                    continue

                nombre = nivel.strip()
                if i < len(etiquetas):
                    sintax_string.append(f"{etiquetas[i]}:{nombre}")

            tax_dict[seq_id] = ",".join(sintax_string)

    # 2. Escribir FASTA de salida compatible con SINTAX
    with open(fasta_out, "w") as fout:
        for record in SeqIO.parse(fasta_in, "fasta"):
            seqid = record.id
            if seqid not in tax_dict:
                continue

            linea = f">{seqid};tax={tax_dict[seqid]};\n"
            fout.write(linea)
            fout.write(str(record.seq) + "\n")

    print("✔ Archivo SINTAX generado:", fasta_out)
