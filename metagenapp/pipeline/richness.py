import pandas as pd


def calcular_riqueza(otu_path):
    df = pd.read_csv(otu_path, sep="\t")
    df = df.set_index(df.columns[0])

    riqueza = (df > 0).sum(axis=0)

    riqueza_df = (
        riqueza.rename_axis("Muestra")
        .reset_index(name="ASV_Richness")
    )

    return riqueza_df
