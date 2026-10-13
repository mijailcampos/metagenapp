#!/bin/bash
# Regenerates the statistics, tables and Figure of the MetagenApp v1.1 benchmark
# (36 libraries) from the pipeline outputs deposited on Zenodo
# (https://doi.org/10.5281/zenodo.23349702).
#
# Requirements: R >= 4.3 with vegan (tested: R 4.3.3, vegan 2.6-8);
#               Python 3 with pandas.
# Usage:        bash run_all.sh            (downloads the outputs if data/ is missing)
#               RSCRIPT=/path/to/Rscript PYTHON=/path/to/python bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
RSCRIPT=${RSCRIPT:-Rscript}
PYTHON=${PYTHON:-python3}
TARBALL=benchmark_outputs_v1.1.tar.gz
SHA256=c2f08c9772665f2e799d6624cf764375eacb120a6a461996323d0251fa7e03ee

if [ ! -d data/samples36 ]; then
  if [ ! -f "$TARBALL" ]; then
    curl -fL -o "$TARBALL" "https://zenodo.org/records/23349702/files/$TARBALL?download=1"
  fi
  echo "$SHA256  $TARBALL" | sha256sum -c -
  mkdir -p data
  tar xzf "$TARBALL" -C data
fi

mkdir -p logs tables figures rarefaction permanova_strata site_permanova
i=1
for s in rarefaction_analysis.R rarefied_diversity_all36.R alpha_diversity_stats_all36.R \
         beta_diversity_pcoa_permanova.R pcoa_variance.R taxonomic_composition.R \
         permanova_strata.R site_permanova.R; do
  log=logs/$(printf "%02d" $i)_${s%.R}.log
  echo "[$i/9] $s"
  "$RSCRIPT" scripts/$s > "$log" 2>&1
  i=$((i+1))
done
echo "[9/9] genus_resolution.py"
"$PYTHON" scripts/genus_resolution.py data/samples36/metagenapp tables/genus_resolution_table6.tsv \
  > logs/09_genus_resolution.log 2>&1
echo "Done. Tables in tables/, figure in figures/, logs in logs/."
