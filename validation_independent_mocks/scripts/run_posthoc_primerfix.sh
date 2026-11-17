#!/bin/bash
# 12 nov 2026: ANÁLISIS POSTERIOR (no preregistrado; plan original = e353b80).
# Código = MetagenApp v1.1 (tag 0a4c842) + ÚNICO cambio: trim_primers.py sin anclar (BIOS-481), min_overlap = len - floor(0.1 len) (2a corrida; la 1a, estricta, en posthoc_primerfix_estricto).
# 1) regresión con los mocks del paper (Aachen, Barcelona); 2) conjuntos A/D/E de la validación.
# Parámetros idénticos a las corridas originales.
set -e
V=/data/validation_independent_mocks_20261111; M=/data/mock_community_20261001; O=$V/posthoc_primerfix
export PATH=/data/envs/metagenapp_env_v2/bin:$PATH
export PYTHONPATH=/data/projects/paper_metagenapp_benchmark/metagenapp_v11_primerfix
export METAGENAPP_REFS=/data/databases/metagenapp_refs
cd $O
python -c "import metagenapp; print(\"USANDO\", metagenapp.__file__)"
run(){ # nombre entrada
  echo "=== $1 $(date -u)"
  /usr/bin/time -v /data/envs/metagenapp_env_v2/bin/metagenapp -i $2 -o $O/$1 \
    --mode ref --marker 16S --classifier wang \
    --primer-f CCTACGGGNGGCWGCAG --primer-r GACTACHVGGGTATCTAATCC \
    --min-length 400 --max-length 440 --max-ambigs 0 --max-poly 8 -t 16 > $O/logs/$1.log 2>&1
  echo "=== $1 OK $(date -u)"
}
run regresion_aachen $M/input/aachen
run regresion_barcelona $M/input/barcelona
run validacion_ADE $V/input
echo "EXIT_OK_0 $(date -u)"
