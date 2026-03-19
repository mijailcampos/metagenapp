# =====================================================
# CONFIGURACIÓN GLOBAL DE RUTAS DE MetagenApp
# Arquitectura portable y desacoplada del working dir
# =====================================================

from pathlib import Path
import os

# =====================================================
# 1️⃣ DIRECTORIO DEL PAQUETE (siempre estable)
# =====================================================

PACKAGE_DIR = Path(__file__).resolve().parent

# =====================================================
# 2️⃣ DIRECTORIO DE DATOS DEL USUARIO (configurable)
# =====================================================

# Permite override por variable de entorno
DEFAULT_DATA_ROOT = PACKAGE_DIR / "user_data"
DATA_ROOT = Path(os.getenv("METAGENAPP_DATA", DEFAULT_DATA_ROOT))

INPUT_DIR = DATA_ROOT / "inputs"
OUTPUT_DIR = DATA_ROOT / "outputs"

INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# 3️⃣ BASES EXTERNAS CENTRALIZADAS (PRODUCCIÓN)
# =====================================================

REFERENCE_ROOT = Path("/data/databases")
METAGEN_REFS = REFERENCE_ROOT / "metagenapp_refs"

# =====================================================
# 16S
# =====================================================

REF_16S_ROOT = METAGEN_REFS / "16S"

NAIVE_MODEL_PATH = "/data/databases/metagenapp_refs/16S/naive_model_v4.pkl"

# Índice kraken-lite — SILVA 138.2 NR99 Bacteria+Archaea (451k seqs, 83k taxa)
KRAKEN_INDEX_PATH = REF_16S_ROOT / "kraken_index_silva.pkl"

# Índice anterior (V4 only, ~7500 taxa) — backup
# KRAKEN_INDEX_PATH = REF_16S_ROOT / "kraken_index_v4.pkl"

HIERARCHY_ROOT_16S = REF_16S_ROOT / "hierarchy"
PHYLUM_MODEL_PATH = HIERARCHY_ROOT_16S / "model_phylum_k7_opt.joblib"
GENUS_MODELS_DIR = HIERARCHY_ROOT_16S / "genus_opt"

# 🔬 Referencias 16S
SILVA_REFERENCE_RAW = REF_16S_ROOT / "silva_reference.fasta"
SILVA_REFERENCE_ALN = REF_16S_ROOT / "silva_reference_aligned.fasta"

# =====================================================
# 18S
# =====================================================

REF_18S_ROOT = METAGEN_REFS / "18S"

PR2_NAIVE_MODEL_PATH = REF_18S_ROOT / "naive_pr2_model.pkl"
PR2_REFERENCE = REF_18S_ROOT / "pr2_reference.fasta"

HIERARCHY_ROOT_18S = REF_18S_ROOT / "hierarchy"  # futuro

# =====================================================
# TRAINSET (si aún lo usas)
# =====================================================

TRAINSET_FASTA = REFERENCE_ROOT / "16S_refseq_rdp" / "trainset9_032012.pds.fasta"
TRAINSET_TAX = REFERENCE_ROOT / "16S_refseq_rdp" / "trainset9_032012.pds.tax"

SINTAX_TRAINSET = TRAINSET_FASTA
TAX_FILE = TRAINSET_TAX

# =====================================================
# 4️⃣ ARCHIVOS DE ENTRADA / SALIDA
# =====================================================

FILES_OUTPUT = INPUT_DIR / "input_samples.files"

CONTIGS_FASTA = OUTPUT_DIR / "assembled_contigs.fasta"
CONTIGS_COUNT_TABLE = OUTPUT_DIR / "assembled_contigs.count_table"

SCREENED_FASTA = OUTPUT_DIR / "screened_contigs.fasta"
SCREENED_COUNT_TABLE = OUTPUT_DIR / "screened_contigs.count_table"

UNIQUE_FASTA = OUTPUT_DIR / "unique_contigs.fasta"
UNIQUE_COUNT_TABLE = OUTPUT_DIR / "unique_contigs.count_table"

ALIGNED_VSEARCH_OUTPUT = OUTPUT_DIR / "aligned_vsearch.aln"
ALIGNED_TRIMMED_FASTA = OUTPUT_DIR / "aligned_trimmed.fasta"

CLUSTERED_FASTA_97 = OUTPUT_DIR / "clustered_97.fasta"
CLUSTERED_COUNT_TABLE_97 = OUTPUT_DIR / "clustered_97.count_table"
UC_CLUSTER_FILE = OUTPUT_DIR / "recortadas_97.uc"

CLASSIFICATION_TAX = OUTPUT_DIR / "classification.taxonomy"

FILTERED_MAFFT_UNIQUE_FASTA = OUTPUT_DIR / "filtered_mafft_unique.fasta"
FILTERED_MAFFT_UNIQUE_COUNT_TABLE = OUTPUT_DIR / "filtered_mafft_unique.count_table"

NON_CHIMERAS_FASTA = OUTPUT_DIR / "non_chimeras.fasta"
NON_CHIMERAS_COUNT_TABLE = OUTPUT_DIR / "non_chimeras.count_table"

FINAL_CLEAN_FASTA = OUTPUT_DIR / "final_clean.fasta"
FINAL_CLEAN_COUNT_TABLE = OUTPUT_DIR / "final_clean.count_table"
FINAL_CLEAN_TAXONOMY = OUTPUT_DIR / "final_clean.taxonomy"

OTU_TABLE_0_03 = OUTPUT_DIR / "otu_table_0_03.tsv"
FINAL_OTU_TAXONOMY = OUTPUT_DIR / "final.otu.taxonomy"

# =====================================================
# 5️⃣ PARÁMETROS GENERALES
# =====================================================

DEFAULT_VSEARCH_THREADS = 16
DEFAULT_CLUSTER_IDENTITY = 0.97
DEFAULT_ALIGN_IDENTITY = 0.97
