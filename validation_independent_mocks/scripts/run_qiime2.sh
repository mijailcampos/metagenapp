#!/bin/bash
# 11 nov 2026: validación independiente (plan e353b80). Copia de run_qiime2_aachen.sh; cambian OUT, manifiesto y trunc-len-f 250 -> 280 (valor del benchmark, como dice el plan).
set -e
export PATH=/data/envs/qiime2/bin:$PATH
OUT=/data/validation_independent_mocks_20261111/outputs/qiime2
mkdir -p $OUT
cd $OUT

echo "[$(date)] 1/5 import manifest"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /data/validation_independent_mocks_20261111/qiime2_manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path $OUT/demux.qza

echo "[$(date)] 2/5 cutadapt trim-paired (remove 341F/805R)"
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $OUT/demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-discard-untrimmed \
  --p-cores 16 \
  --o-trimmed-sequences $OUT/demux_trimmed.qza \
  --verbose

echo "[$(date)] 3/5 dada2 denoise-paired"
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $OUT/demux_trimmed.qza \
  --p-trunc-len-f 280 --p-trunc-len-r 220 \
  --p-n-threads 16 \
  --o-representative-sequences $OUT/rep-seqs.qza \
  --o-table $OUT/table.qza \
  --o-denoising-stats $OUT/stats.qza \
  --verbose

echo "[$(date)] 4/5 classify-sklearn"
qiime feature-classifier classify-sklearn \
  --i-classifier /data/results_v2/qiime2_classifier_138.2/silva-138.2-99-nb-classifier.qza \
  --i-reads $OUT/rep-seqs.qza \
  --p-n-jobs 16 \
  --o-classification $OUT/taxonomy.qza \
  --verbose

echo "[$(date)] 5/5 export"
qiime tools export --input-path $OUT/table.qza --output-path $OUT/table_export
qiime tools export --input-path $OUT/taxonomy.qza --output-path $OUT/taxonomy_export
qiime tools export --input-path $OUT/rep-seqs.qza --output-path $OUT/repseqs_export

echo "[$(date)] DONE"
