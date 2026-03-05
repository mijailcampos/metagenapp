import importlib
import pkgutil
import metagenapp_core.models as models_pkg


def discover_models():

    models = {}

    for _, module_name, _ in pkgutil.iter_modules(models_pkg.__path__):

        if module_name.startswith("_"):
            continue

        module = importlib.import_module(
            f"metagenapp_core.models.{module_name}"
        )

        if hasattr(module, "MODEL_NAME"):
            models[module.MODEL_NAME] = module

    return models


MODELS = discover_models()
