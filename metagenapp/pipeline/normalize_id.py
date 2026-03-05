def normalize_id(seq_id):
    """
    Normaliza los identificadores de secuencia para una comparación consistente.
    Esta función intenta extraer la parte principal del ID, manejando los patrones
    más comunes que se han observado en los datos del usuario.

    Args:
        seq_id (str): El identificador de secuencia original.

    Returns:
        str: El identificador de secuencia normalizado.
    """
    # Patrón para IDs con sufijo de versión (ej. AF515816.1)
    if '.' in seq_id:
        return seq_id.split('.')[0]
    # Patrón para IDs con sufijo de conjunto de datos (ej. EU861894_S001148199)
    # Esto asume que la parte relevante es antes de '_S00'
    if '_S00' in seq_id:
        return seq_id.split('_S00')[0]
    # Si no coincide con los patrones específicos, devuelve el ID original.
    return seq_id
