def parse_taxonomy(taxon_str):
    """
    Convierte taxonomías de distintos formatos (Mothur, SINTAX, PR2, SILVA)
    a los niveles estándar:
    Kingdom, Phylum, Class, Order, Family, Genus, Species
    """

    niveles = ["Kingdom", "Phylum", "Class",
               "Order", "Family", "Genus", "Species"]
    resultado = {n: "Unclassified" for n in niveles}

    if not isinstance(taxon_str, str) or len(taxon_str.strip()) == 0:
        return resultado

    # ======================================
    # 1. Normalizar string
    # ======================================
    taxon_str = taxon_str.replace('"', "").replace(" ", "").strip()

    # ======================================
    # 2. Detectar formato
    # ======================================

    # Formato SINTAX → separado por comas y con ":score"
    if "," in taxon_str and ":" in taxon_str:
        partes = [p.split(":")[0].strip("_") for p in taxon_str.split(",")]

    # Formato Mothur/SILVA/PR2 → separado por ;
    else:
        partes = [p for p in taxon_str.split(";") if p != ""]

    # ======================================
    # 3. Asignación directa por posición si existe
    # ======================================
    for i, nivel in enumerate(niveles):
        if i < len(partes):
            resultado[nivel] = partes[i]

    # ======================================
    # 4. Rangos 'Incertae_Sedis' (SILVA): no son un taxón; se etiquetan con el
    # rango superior con nombre, p. ej. 'Muribaculaceae (incertae sedis)', para
    # que no se junten familias distintas en un falso género (BIOS-478/479).
    # ======================================
    for i, nivel in enumerate(niveles[1:], start=1):
        if resultado[nivel].lower().startswith("incertae_sedis"):
            padre = next((resultado[n] for n in reversed(niveles[:i])
                          if resultado[n] != "Unclassified"
                          and not resultado[n].lower().startswith("incertae_sedis")), None)
            if padre:
                resultado[nivel] = f"{padre.split(' (')[0]} (incertae sedis)"

    return resultado
