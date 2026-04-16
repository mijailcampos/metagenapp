# MetagenApp

Pipeline de **metabarcoding 16S/18S** desarrollado como alternativa a Mothur y QIIME2.
Clasifica taxonómicamente secuencias de ADN ribosomal desde muestras microbianas,
produciendo tablas de OTUs/ASVs con anotación taxonómica.

Desarrollado por Mijail Campos.

---

## Ventajas frente a Mothur / QIIME2

- Clasificador propio (`naive-v2`) basado en índice invertido de k-mers con confianza estilo Wang
- Extremadamente eficiente en recursos: **9 min y 1.67 GB RAM** vs 97 min / 49.9 GB de Mothur (mismo dataset)
- Soporte para múltiples modelos de referencia: general, oral, gut, skin, env
- CLI moderna con perfiles configurables

---

## Estructura del Proyecto

```
MetagenApp/
├── metagenapp/                    # Capa CLI y orquestación del pipeline
│   ├── cli/main.py                # Entry point: CLI con Typer
│   ├── pipeline/
│   │   ├── run_pipeline.py        # Orquestador: 26 pasos secuenciales
│   │   ├── clasificacion_tax.py   # Coordinación de clasificadores
│   │   ├── aggregate_cluster_counts.py  # Suma reads UC → centroides
│   │   └── [módulos por paso]
│   └── metagen_config.py          # Configuración global y rutas
│
├── metagenapp_core/               # Motores de clasificación
│   ├── models/
│   │   ├── naive_v2.py            # Clasificador naive-v2 (paralelismo fork)
│   │   ├── naive_v2_engine.py     # Motor: Wang-style bootstrap confidence
│   │   ├── kraken_lite.py         # Wrapper paralelo kraken-lite
│   │   └── kraken_lite_engine.py  # Motor: votación k-mer + LCA
│   └── utils/
│       ├── kmers.py               # Generación de k-mers
│       └── taxonomy.py            # Cálculo de LCA
│
└── scripts/                       # Entrenamiento, benchmarks, figuras
    ├── train_naive_v2.py
    ├── build_homd_trainset.py
    ├── Figure1_phylum_comparison_naivev2.R
    └── Figure2_genus_bacteroidetes.R
```

---

## Pipeline: 26 pasos

Toma **pares de FASTQ** (paired-end Illumina) y produce **tablas de abundancia con taxonomía**.

```
FASTQ R1 + R2
    ↓
[01] Ensamblado (VSEARCH mergepairs)
    ↓
[02] Filtrado (largo 250–600 bp, sin ambigüedades)
    ↓
[03] Unicización (deduplicación + conteo por muestra)
    ↓
[04] Alineamiento a referencia (VSEARCH, 70% identidad)
    ↓
[09] Recorte de la región de interés
    ↓
[10] Clustering 97% → centroides + archivo UC
    ↓
[11] Extracción de centroides (modo-dependiente)
    ↓
[12–18] Alineamiento MAFFT + filtrado + chimeras
         (solo modes student/premium/turbo)
    ↓
[20] Clasificación taxonómica  ← corazón del pipeline
    ↓
[21] Remoción de linajes (cloroplastos, mitocondrias)
    ↓
[22] Resumen por filo (summary_tax_phylum.tsv)
    ↓
[23–26] Tablas OTU/ASV finales
```

### Modos de ejecución

| Modo    | Descripción                          | Centroides máx |
|---------|--------------------------------------|----------------|
| student | Rápido, bajo consumo                 | 10,000         |
| premium | Balanceado                           | Todos          |
| turbo   | Máximo rendimiento                   | Todos          |
| ref     | Grado publicación (EDLib, sin MAFFT) | Todos          |
| qa      | Solo control de calidad              | —              |

---

## Clasificador naive-v2

### Algoritmo (Wang-style bootstrap)

1. Extracción de k-mers de la secuencia query
2. Filtrado: descartar k-mers con `psize ≥ 20` (demasiado genéricos)
3. Clasificación real: sumar pesos por taxón en cada nivel taxonómico
4. Bootstrap (100 iteraciones): remuestrear pool de k-mers con reemplazo → contar qué taxón gana
5. Asignación al nivel taxonómico más profundo con confianza ≥ 80%

**Parámetros clave:**
- `psize_max=20` — ignorar k-mers presentes en más de 20 taxa
- `confidence=0.80` — umbral de confianza para asignar
- `n_bootstrap=100` — iteraciones de remuestreo

### Modelos disponibles

| Modelo | Taxa | K-mers únicos | Ruta |
|--------|------|---------------|------|
| General v4 | 1,949 | 1,606,675 | `/data/databases/metagenapp_refs/16S/naive_model_v4.pkl` |
| Oral v1 (HOMD) | 802 | 343,773 | `/data/databases/metagenapp_refs/16S/models/naive_model_oral_v1.pkl` |

### Resolución a nivel de género

| Método | Resolución a género |
|--------|---------------------|
| MetagenApp LCA (anterior) | 42.4% |
| MetagenApp Wang bootstrap (actual) | **74.5%** |
| QIIME2 | ~98% |
| Mothur | ~100% |

La brecha se debe a cobertura del reference (1,949 taxa vs ~50,000 en SILVA).
Con modelo SILVA-scale se espera >95%.

---

## Clasificador Kraken-lite

Implementación propia del algoritmo de Kraken:

1. Genera k-mers (k=31, codificados como enteros 2-bit)
2. Busca en índice → votos por taxón
3. Filtra ruido (votos < 0.5% del total)
4. Voto mayoritario jerárquico: Order → Phylum
5. Si confianza ≥ 0.30 → asignar; si no, LCA de top 10 taxa

**Modelo:** SILVA 138.2 NR99, 451k secuencias, 83k taxa
**Ruta:** `/data/databases/metagenapp_refs/16S/kraken_index_silva.pkl`

---

## Uso (CLI)

```bash
metagenapp --input ./fastq --outdir ./results [opciones]

# Opciones principales:
--mode        student|premium|turbo|ref|qa
--classifier  flat|naive-v2|kraken-lite|pro-engine
--marker      16S|18S
--model-type  general|oral|gut|skin|env
--threads     N (default: 8)
--from-step   NOMBRE_PASO  # reanudar desde un paso específico

# Con perfil (auto-genera directorio de salida con timestamp):
metagenapp --profile ref --input ./fastq
```

### profiles.yaml

```yaml
profiles:
  ref:
    mode: ref
    threads: 16
    output_base: /data/results/runs/ref_runs
```

---

## Benchmark: MetagenApp vs Mothur vs QIIME2

Dataset: faringe (36 muestras, 1,519,852 reads entrada).

### Comparación de pipelines

| Parámetro | MetagenApp | Mothur | QIIME2 |
|-----------|------------|--------|--------|
| Algoritmo clustering | VSEARCH 97% | OptiClust 97% | DADA2 (ASV) |
| Reads retenidos | 1,147,591 (75.5%) | 1,145,172 (75.4%) | 191,568 (12.6%) |
| OTUs / ASVs | 15,080 | 10,473 | 2,098 |
| DB taxonómica | SILVA 138 (k-mer) | SILVA 138 (Wang) | SILVA 138 (sklearn NB) |

### Recursos computacionales

| Recurso | MetagenApp | Mothur | QIIME2 |
|---------|------------|--------|--------|
| Tiempo total | **9 min 03s** | 1h 37m 21s | 27 min 24s |
| RAM máxima | **1.67 GB** | 49.9 GB | 14.0 GB |
| CPUs | 16 | multicore | multicore |

### Distribución de filos (read-weighted)

| Filo | QIIME2 | Mothur | MetagenApp oral | MetagenApp general |
|------|--------|--------|-----------------|-------------------|
| Proteobacteria | 35.9% | 36.4% | 34.7% | 34.4% |
| Firmicutes | 38.6% | 37.3% | 33.5% | 37.2% |
| Bacteroidetes | 14.2% | 15.5% | 18.5% | 18.5% |
| Actinobacteria | 5.0% | 5.5% | 4.6% | 4.5% |
| Fusobacteria | 4.2% | 4.2% | 4.1% | 4.3% |

Los tres pipelines concuerdan a nivel de filo. No hay sesgo en MetagenApp.
