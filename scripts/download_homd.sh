#!/usr/bin/env bash
# ============================================================
# Descarga datos HOMD 16S rRNA para entrenamiento oral model
# Human Oral Microbiome Database — homd.org
# ============================================================
# USO:
#   bash scripts/download_homd.sh
#
# Descarga en: /data/databases/metagenapp_refs/16S/homd_raw/
# ============================================================

set -euo pipefail

OUTDIR="/data/databases/metagenapp_refs/16S/homd_raw"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

HOMD_VERSION="V16.03"
HOMD_BASE="https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/${HOMD_VERSION}"

FASTA_URL="$HOMD_BASE/HOMD_16S_rRNA_RefSeq_${HOMD_VERSION}.fasta"
TAX_URL="$HOMD_BASE/HOMD_16S_rRNA_RefSeq_${HOMD_VERSION}.mothur.taxonomy"

echo "=== Descargando HOMD ${HOMD_VERSION} ==="
echo ""
echo "FASTA:    $FASTA_URL"
echo "Taxonomy: $TAX_URL"
echo ""

wget -q --show-progress -c "$FASTA_URL" -O "$OUTDIR/homd_16S.fasta" &
wget -q --show-progress -c "$TAX_URL"   -O "$OUTDIR/homd_16S.taxonomy" &
wait

echo ""
echo "=== Verificación ==="
echo "FASTA:    $(grep -c '^>' homd_16S.fasta) secuencias"
echo "Taxonomy: $(wc -l < homd_16S.taxonomy) líneas"
echo ""
echo "Siguiente paso:"
echo "  python3 scripts/build_homd_trainset.py"
