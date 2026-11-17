#!/bin/bash
# 12 nov 2026: REEVALUACION POSTERIOR, NO CIEGA (Amendment 1, commit 14c6a1c; plan original e353b80).
# Codigo = MetagenApp v1.2.3 (tag, commit b97a9c6), copia verificada byte a byte (250 archivos).
# Valores predeterminados de v1.2.3 (chimera-mode per-sample). Solo conjuntos A/D/E. Mismos parametros de linea de comandos.
set -e
V=/data/validation_independent_mocks_20261111; O=$V/reeval_v123
export PATH=/data/envs/metagenapp_env_v2/bin:$PATH
export PYTHONPATH=/data/projects/paper_metagenapp_benchmark/metagenapp_v123_tag
export METAGENAPP_REFS=/data/databases/metagenapp_refs
cd $O
python -c "import metagenapp; print(\"USANDO\", metagenapp.__file__)"
echo "=== validacion_ADE $(date -u)"
/usr/bin/time -v /data/envs/metagenapp_env_v2/bin/metagenapp -i $V/input -o $O/validacion_ADE \
  --mode ref --marker 16S --classifier wang \
  --primer-f CCTACGGGNGGCWGCAG --primer-r GACTACHVGGGTATCTAATCC \
  --min-length 400 --max-length 440 --max-ambigs 0 --max-poly 8 -t 16 > $O/logs/validacion_ADE.log 2>&1
echo "=== validacion_ADE OK $(date -u)"
echo "EXIT_OK_0 $(date -u)"
