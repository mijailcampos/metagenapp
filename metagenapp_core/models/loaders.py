import pickle


def load_model(path):

    with open(path, "rb") as f:
        model = pickle.load(f)

    # Normalizar estructura del modelo
    if isinstance(model, tuple) and len(model) == 3:

        model_dict = {
            "kmer_counts": model[0],
            "taxonomy_labels": model[1],
            "kmer_size": model[2],
        }

        return model_dict

    return model
