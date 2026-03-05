def convertir_trainset_a_sintax(trainset_fasta, trainset_tax, salida_sintax_fasta):

    tax = {}
    with open(trainset_tax) as f:
        for line in f:
            seq_id, taxonomy = line.strip().split("\t")
            parts = taxonomy.split(";")
            sintax = []
            ranks = ["k", "p", "c", "o", "f", "g", "s"]
            for i, p in enumerate(parts):
                if p:
                    sintax.append(f"{ranks[i]}:{p}")
            tax[seq_id] = ";".join(sintax)

    with open(salida_sintax_fasta, "w") as out:
        with open(trainset_fasta) as f:
            seq_id = None
            seq = []
            for line in f:
                if line.startswith(">"):
                    if seq_id:
                        out.write(f">{seq_id};tax={tax[seq_id]};\n")
                        out.write("".join(seq) + "\n")
                    seq_id = line[1:].strip()
                    seq = []
                else:
                    seq.append(line.strip())

            if seq_id:
                out.write(f">{seq_id};tax={tax[seq_id]};\n")
                out.write("".join(seq) + "\n")
