#!/bin/bash
# 11 nov 2026: validación independiente (plan e353b80). Batch = mothur_aachen_seed.batch con solo rutas cambiadas; entradas descomprimidas.
cd /data/validation_independent_mocks_20261111
/usr/bin/time -v /data/envs/mothur_env/bin/mothur /data/validation_independent_mocks_20261111/mothur_validation.batch
echo "EXIT_OK_$? $(date -u)"
