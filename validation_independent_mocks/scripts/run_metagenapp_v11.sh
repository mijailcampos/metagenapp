#!/bin/bash
# 11 nov 2026: validación independiente (plan congelado e353b80). MetagenApp v1.1 = tag 0a4c842 (verificado byte a byte).
set -e
B=/data/validation_independent_mocks_20261111
export PATH=/data/envs/metagenapp_env_v2/bin:$PATH
export PYTHONPATH=/data/projects/paper_metagenapp_benchmark/metagenapp_v11_tag
export METAGENAPP_REFS=/data/databases/metagenapp_refs
cd $B
python -c "import metagenapp; print(metagenapp.__file__)"
/usr/bin/time -v /data/envs/metagenapp_env_v2/bin/metagenapp \
  -i $B/input -o $B/outputs/metagenapp_v11 \
  --mode ref --marker 16S --classifier wang \
  --primer-f CCTACGGGNGGCWGCAG --primer-r GACTACHVGGGTATCTAATCC \
  --min-length 400 --max-length 440 \
  --max-ambigs 0 --max-poly 8 -t 16
echo "EXIT_OK $(date -u)"
