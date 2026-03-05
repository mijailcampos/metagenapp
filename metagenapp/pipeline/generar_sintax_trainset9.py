import csv

# === ENTRADAS ORIGINALES ===
fasta_in = "trainset9_032012.pds.fasta"
tax_in = "trainset9_032012.pds.tax"

# === SALIDA SINTAX ===
sintax_out = "trainset9_032012.sintax.fasta"


# -------------------------------------------------------------
# 1. Cargar taxonomía trainset9
# -------------------------------------------------------------
taxonomy = {}

with open(tax_in, "r", encoding="utf-8") as taxfile:
    reader = csv.reader(taxfile, delimiter="\t")
    for row in reader:
        seq_id = row[0]
        taxa_raw = row[1].strip(";")

        # Convertir taxonomía a formato SINTAX
        niveles = taxa_raw.split(";")

        # Colocar prefijos SINTAX
        prefix = ["k", "p", "c", "o", "f", "g"]
        sintax_tax = []

        for idx, taxon in enumerate(niveles):
            if idx >= len(prefix):
                break
            clean = taxon.replace('"', '').replace(" ", "_")
            sintax_tax.append(f"{prefix[idx]}:{clean}")

        taxonomy[seq_id] = ";".join(sintax_tax) + ";"


# -------------------------------------------------------------
# 2. Leer FASTA original y generar FASTA SINTAX
# -------------------------------------------------------------
def fasta_reader(path):
    with open(path, "r") as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    yield seq_id, "".join(seq)
                seq_id = line[1:].split()[0]   # Solo ID antes del espacio
                seq = []
            else:
                seq.append(line)
        if seq_id:
            yield seq_id, "".join(seq)


with open(sintax_out, "w") as out:
    for seq_id, sequence in fasta_reader(fasta_in):

        if seq_id not in taxonomy:
            print(f"⚠ Advertencia: {seq_id} no está en el archivo TAX")
            continue

        tax_string = taxonomy[seq_id]

        # Encabezado SINTAX correcto
        header = f">{seq_id};tax={tax_string}"

        out.write(header + "\n")
        out.write(sequence + "\n")

print("\n✅ Archivo generado exitosamente:")
print(f"   {sintax_out}")
