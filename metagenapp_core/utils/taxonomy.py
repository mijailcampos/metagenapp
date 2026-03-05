def lca_taxa(taxa_list):

    if not taxa_list:
        return "Unclassified"

    # limpiar ; final
    split_taxa = [t.rstrip(";").split(";") for t in taxa_list]

    lca = []

    for levels in zip(*split_taxa):

        # ignorar niveles vacíos
        levels = [x for x in levels if x]

        if len(set(levels)) == 1:
            lca.append(levels[0])
        else:
            break

    if not lca:
        return "Unclassified"

    return ";".join(lca) + ";"
