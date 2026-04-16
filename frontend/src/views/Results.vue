<script setup>
import { ref, onMounted } from 'vue'
import { useRouter } from 'vue-router'
import { api } from '../api/client.js'

const props  = defineProps({ jobId: String })
const router = useRouter()

const results = ref(null)
const job     = ref(null)
const error   = ref(null)

// ── EcoLab ──────────────────────────────────────────────────────────────────
const ANALYSES = [
  { key: 'alpha',       label: 'Diversidad Alfa',    icon: '📊', desc: 'Shannon · Simpson · Riqueza · Pielou' },
  { key: 'rarefaction', label: 'Curvas Rarefacción', icon: '📈', desc: 'Profundidad vs riqueza por muestra' },
  { key: 'nmds',        label: 'Beta NMDS',          icon: '🔵', desc: 'Bray-Curtis ordenación en 2D' },
  { key: 'pcoa',        label: 'PCoA',               icon: '🔷', desc: 'Coordenadas principales Bray-Curtis' },
  { key: 'permanova',   label: 'PERMANOVA',          icon: '📐', desc: 'Prueba estadística entre grupos (adonis2)', needsGroups: true },
  { key: 'phylum',      label: 'Abundancia Filo',    icon: '🌿', desc: 'Top 10 phyla barplot apilado' },
  { key: 'genus',       label: 'Abundancia Género',  icon: '🔬', desc: 'Top 15 géneros barplot apilado' },
  { key: 'rarefy',      label: 'Rarefy Tabla',       icon: '⚖',  desc: 'Normaliza al mínimo de reads' },
]

const ecoState = ref(
  Object.fromEntries(ANALYSES.map(a => [a.key, {
    loading: false,
    result:  null,
    error:   null,
    phase:   'idle',   // solo para needsGroups: 'idle' | 'setup' | 'done'
    samples: [],
    groups:  {},
  }]))
)

function uniqueGroups(key) {
  return [...new Set(Object.values(ecoState.value[key].groups).map(v => v.trim()).filter(Boolean))]
}

async function runAnalysis(key) {
  const s    = ecoState.value[key]
  const meta = ANALYSES.find(a => a.key === key)

  if (meta.needsGroups) {
    if (s.phase === 'idle') {
      // Paso 1: cargar nombres de muestras
      s.loading = true
      s.error   = null
      try {
        const r  = await api.getEcolabSamples(props.jobId)
        s.samples = r.samples
        s.groups  = Object.fromEntries(r.samples.map(n => [n, '']))
        s.phase   = 'setup'
      } catch (e) {
        s.error = e.message
      } finally {
        s.loading = false
      }
    } else if (s.phase === 'setup') {
      // Paso 2: validar y ejecutar
      const vals = s.samples.map(n => s.groups[n]?.trim() ?? '')
      if (vals.some(v => !v)) { s.error = 'Asigna un grupo a cada muestra'; return }
      if (new Set(vals).size < 2) { s.error = 'Se necesitan al menos 2 grupos distintos'; return }
      s.loading = true
      s.error   = null
      try {
        s.result = await api.runEcolab(props.jobId, key, { groups: vals })
        s.phase  = 'done'
      } catch (e) {
        s.error = e.message
      } finally {
        s.loading = false
      }
    }
    return
  }

  // Análisis normal (sin grupos)
  s.loading = true
  s.error   = null
  s.result  = null
  try {
    s.result = await api.runEcolab(props.jobId, key)
  } catch (e) {
    s.error = e.message
  } finally {
    s.loading = false
  }
}

// Metadatos de cada archivo para el usuario
const FILE_META = {
  'final_asv_table.tsv': {
    label: 'Tabla ASV',
    desc:  'Matriz abundancia ASV × muestras',
    tools: ['R', 'STAMP'],
    icon:  '📊',
  },
  'final_asv.taxonomy': {
    label: 'Taxonomía ASV',
    desc:  'Asignación taxonómica por ASV',
    tools: ['R'],
    icon:  '🏷',
  },
  'asv_resumen.tsv': {
    label: 'Resumen ASV',
    desc:  'Conteos por muestra desglosados',
    tools: ['R', 'STAMP'],
    icon:  '📋',
  },
  'otu_table_0_03.tsv': {
    label: 'Tabla OTU (97%)',
    desc:  'OTUs agrupados al 97% de identidad',
    tools: ['R', 'STAMP'],
    icon:  '🔬',
  },
  'summary_tax_phylum.tsv': {
    label: 'Resumen Phylum',
    desc:  'Abundancia relativa por phylum',
    tools: ['R', 'STAMP'],
    icon:  '🌿',
  },
  'final_clean.fasta': {
    label: 'Secuencias finales',
    desc:  'FASTA con los centroides clasificados',
    tools: ['QIIME2', 'BLAST'],
    icon:  '🧬',
  },
  'final_clean.taxonomy': {
    label: 'Taxonomía final',
    desc:  'Formato compatible con Mothur/QIIME2',
    tools: ['Mothur', 'QIIME2'],
    icon:  '📄',
  },
}

function humanSize(bytes) {
  if (!bytes) return '—'
  if (bytes < 1024)       return bytes + ' B'
  if (bytes < 1024 ** 2)  return (bytes / 1024).toFixed(1) + ' KB'
  return (bytes / 1024 ** 2).toFixed(1) + ' MB'
}

onMounted(async () => {
  try {
    const [r, j] = await Promise.all([
      api.getResults(props.jobId),
      api.getJob(props.jobId),
    ])
    results.value = r
    job.value     = j
  } catch (e) {
    error.value = e.message
  }
})
</script>

<template>
  <div>
    <div class="page-header">
      <h1 class="page-title">Resultados</h1>
      <span class="job-id">ID: {{ jobId }}</span>
    </div>

    <p v-if="error" class="error-msg">{{ error }}</p>

    <template v-if="results">
      <!-- Métricas rápidas -->
      <section v-if="job" class="metrics-row">
        <div class="metric">
          <span class="metric-val">{{ job.elapsed_s ? Math.round(job.elapsed_s / 60) + 'm' : '—' }}</span>
          <span class="metric-label">Tiempo total</span>
        </div>
        <div class="metric">
          <span class="metric-val">{{ job.config?.mode }}</span>
          <span class="metric-label">Modo</span>
        </div>
        <div class="metric">
          <span class="metric-val">{{ job.config?.classifier }}</span>
          <span class="metric-label">Clasificador</span>
        </div>
        <div class="metric">
          <span class="metric-val">{{ job.config?.marker }}</span>
          <span class="metric-label">Marcador</span>
        </div>
      </section>

      <!-- Archivos -->
      <section class="card">
        <h2 class="card-title">Archivos disponibles</h2>

        <p v-if="!Object.keys(results.files).length" class="empty">
          No hay archivos de resultado disponibles.
        </p>

        <div class="files-grid">
          <div
            v-for="(info, filename) in results.files"
            :key="filename"
            class="file-card"
          >
            <div class="file-icon">{{ FILE_META[filename]?.icon || '📁' }}</div>
            <div class="file-info">
              <p class="file-label">{{ FILE_META[filename]?.label || filename }}</p>
              <p class="file-desc">{{ FILE_META[filename]?.desc || '' }}</p>
              <p class="file-size">{{ humanSize(info.size_bytes) }}</p>
              <div class="file-tools">
                <span
                  v-for="tool in (FILE_META[filename]?.tools || [])"
                  :key="tool"
                  class="tool-badge"
                >{{ tool }}</span>
              </div>
            </div>
            <a
              class="btn-download"
              :href="api.downloadUrl(jobId, filename)"
              :download="filename"
            >
              ⬇ Descargar
            </a>
          </div>
        </div>
      </section>

      <!-- EcoLab -->
      <section class="card ecolab-card">
        <h2 class="card-title">Analizar con EcoLab</h2>
        <p class="ecolab-subtitle">Análisis ecológico directo sobre los resultados — powered by R/vegan</p>

        <div class="ecolab-grid">
          <div v-for="a in ANALYSES" :key="a.key" class="ecolab-item">
            <button
              class="btn-ecolab"
              :class="{ 'btn-ecolab--loading': ecoState[a.key].loading }"
              :disabled="ecoState[a.key].loading"
              @click="runAnalysis(a.key)"
            >
              <span class="eco-btn-icon">{{ a.icon }}</span>
              <span class="eco-btn-content">
                <span class="eco-btn-label">{{ a.label }}</span>
                <span class="eco-btn-desc">{{ a.desc }}</span>
              </span>
              <span v-if="ecoState[a.key].loading" class="eco-spinner">⏳</span>
              <span v-else-if="ecoState[a.key].result" class="eco-done">✓</span>
            </button>

            <!-- Error -->
            <p v-if="ecoState[a.key].error" class="eco-error">
              {{ ecoState[a.key].error }}
            </p>

            <!-- PERMANOVA: formulario de grupos -->
            <div v-if="a.needsGroups && ecoState[a.key].phase === 'setup'" class="perm-setup">
              <p class="perm-hint">
                Asigna un grupo a cada muestra. Las muestras del mismo grupo deben tener el mismo nombre.
              </p>
              <div class="perm-samples">
                <div v-for="sample in ecoState[a.key].samples" :key="sample" class="perm-row">
                  <span class="perm-sample">{{ sample }}</span>
                  <input
                    v-model="ecoState[a.key].groups[sample]"
                    class="perm-input"
                    placeholder="ej: Control"
                    :list="'perm-groups-' + a.key"
                  />
                </div>
              </div>
              <datalist :id="'perm-groups-' + a.key">
                <option v-for="g in uniqueGroups(a.key)" :key="g" :value="g" />
              </datalist>
              <button
                class="btn-perm-run"
                :disabled="ecoState[a.key].loading"
                @click="runAnalysis(a.key)"
              >
                {{ ecoState[a.key].loading ? 'Calculando…' : 'Ejecutar PERMANOVA' }}
              </button>
            </div>

            <!-- Resultado tabla simple -->
            <div v-if="ecoState[a.key].result?.type === 'table'" class="eco-table-wrap">
              <div class="eco-table-scroll">
                <table class="eco-table">
                  <thead><tr><th v-for="col in ecoState[a.key].result.columns" :key="col">{{ col }}</th></tr></thead>
                  <tbody>
                    <tr v-for="(row, i) in ecoState[a.key].result.rows" :key="i">
                      <td v-for="(val, j) in row" :key="j">{{ val }}</td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </div>

            <!-- Resultado PERMANOVA (dos tablas) -->
            <template v-if="ecoState[a.key].result?.type === 'permanova'">
              <div class="eco-table-wrap">
                <p class="eco-table-title">PERMANOVA (adonis2)</p>
                <div class="eco-table-scroll">
                  <table class="eco-table">
                    <thead><tr><th v-for="col in ecoState[a.key].result.permanova.columns" :key="col">{{ col }}</th></tr></thead>
                    <tbody>
                      <tr v-for="(row, i) in ecoState[a.key].result.permanova.rows" :key="i">
                        <td v-for="(val, j) in row" :key="j"
                          :class="{ 'pval-sig': j === 5 && parseFloat(val) < 0.05 }">{{ val }}</td>
                      </tr>
                    </tbody>
                  </table>
                </div>
              </div>
              <div class="eco-table-wrap" style="margin-top:0.6rem">
                <p class="eco-table-title">Homogeneidad de dispersión (betadisper)</p>
                <div class="eco-table-scroll">
                  <table class="eco-table">
                    <thead><tr><th v-for="col in ecoState[a.key].result.betadisper.columns" :key="col">{{ col }}</th></tr></thead>
                    <tbody>
                      <tr v-for="(row, i) in ecoState[a.key].result.betadisper.rows" :key="i">
                        <td v-for="(val, j) in row" :key="j"
                          :class="{ 'pval-sig': j === 5 && parseFloat(val) < 0.05 }">{{ val }}</td>
                      </tr>
                    </tbody>
                  </table>
                </div>
              </div>
            </template>

            <!-- Resultado imagen -->
            <div v-if="ecoState[a.key].result?.type === 'image'" class="eco-img-wrap">
              <img
                :src="api.ecolabImageUrl(jobId, ecoState[a.key].result.filename)"
                class="eco-img"
                alt="Gráfica EcoLab"
              />
            </div>
          </div>
        </div>
      </section>

      <button class="btn-secondary" @click="router.push('/')">
        ← Nueva corrida
      </button>
    </template>

    <div v-else-if="!error" class="loading">Cargando resultados…</div>
  </div>
</template>

<style scoped>
.page-header { display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1.5rem; }
.page-title  { font-size: 1.6rem; font-weight: 700; color: #f1f5f9; }
.job-id      { font-size: 0.8rem; color: #475569; font-family: monospace; }

/* Métricas */
.metrics-row {
  display: flex;
  gap: 1rem;
  margin-bottom: 1.5rem;
  flex-wrap: wrap;
}
.metric {
  flex: 1;
  min-width: 100px;
  background: #1a1d27;
  border: 1px solid #2d3148;
  border-radius: 10px;
  padding: 1rem;
  text-align: center;
}
.metric-val   { display: block; font-size: 1.4rem; font-weight: 700; color: #7c83fd; }
.metric-label { display: block; font-size: 0.75rem; color: #64748b; margin-top: 0.2rem; }

.card {
  background: #1a1d27;
  border: 1px solid #2d3148;
  border-radius: 12px;
  padding: 1.5rem;
  margin-bottom: 1.5rem;
}
.card-title {
  font-size: 0.8rem;
  font-weight: 600;
  color: #64748b;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  margin-bottom: 1.2rem;
}

/* Files */
.files-grid { display: flex; flex-direction: column; gap: 0.8rem; }

.file-card {
  display: flex;
  align-items: center;
  gap: 1rem;
  background: #12151f;
  border: 1px solid #1e2235;
  border-radius: 8px;
  padding: 1rem;
}
.file-icon  { font-size: 1.8rem; flex-shrink: 0; }
.file-info  { flex: 1; }
.file-label { font-size: 0.95rem; font-weight: 600; color: #e2e8f0; }
.file-desc  { font-size: 0.8rem; color: #64748b; margin: 0.15rem 0; }
.file-size  { font-size: 0.75rem; color: #475569; }

.file-tools { display: flex; gap: 0.3rem; margin-top: 0.4rem; flex-wrap: wrap; }
.tool-badge {
  font-size: 0.7rem;
  background: rgba(124,131,253,0.15);
  color: #7c83fd;
  padding: 0.1rem 0.4rem;
  border-radius: 4px;
  font-weight: 600;
}

.btn-download {
  flex-shrink: 0;
  background: #7c83fd;
  color: #fff;
  text-decoration: none;
  padding: 0.5rem 1rem;
  border-radius: 7px;
  font-size: 0.85rem;
  font-weight: 600;
  transition: background 0.15s;
  white-space: nowrap;
}
.btn-download:hover { background: #6366f1; }

.btn-secondary {
  background: #1e2235;
  border: 1px solid #2d3148;
  color: #94a3b8;
  border-radius: 6px;
  padding: 0.5rem 1rem;
  font-size: 0.85rem;
  cursor: pointer;
}
.btn-secondary:hover { border-color: #7c83fd; color: #7c83fd; }

.error-msg { color: #f87171; font-size: 0.85rem; }
.empty     { color: #475569; font-style: italic; }
.loading   { color: #475569; padding: 2rem; text-align: center; }

/* ── EcoLab ────────────────────────────────────────────────────────────────── */
.ecolab-card     { border-color: #1e3a3a; }
.ecolab-subtitle { font-size: 0.8rem; color: #4b7a7a; margin: -0.6rem 0 1.2rem; }

.ecolab-grid {
  display: flex;
  flex-direction: column;
  gap: 1.2rem;
}

.ecolab-item { display: flex; flex-direction: column; gap: 0.6rem; }

.btn-ecolab {
  display: flex;
  align-items: center;
  gap: 0.8rem;
  background: #0f1f1f;
  border: 1px solid #1e3a3a;
  border-radius: 8px;
  padding: 0.75rem 1rem;
  cursor: pointer;
  transition: border-color 0.15s, background 0.15s;
  text-align: left;
  width: 100%;
  color: inherit;
}
.btn-ecolab:hover:not(:disabled) { border-color: #00bfbf; background: #0a1a1a; }
.btn-ecolab--loading { opacity: 0.7; cursor: wait; }
.btn-ecolab:disabled { cursor: not-allowed; }

.eco-btn-icon    { font-size: 1.4rem; flex-shrink: 0; }
.eco-btn-content { flex: 1; display: flex; flex-direction: column; gap: 0.1rem; }
.eco-btn-label   { font-size: 0.9rem; font-weight: 600; color: #e2e8f0; }
.eco-btn-desc    { font-size: 0.75rem; color: #4b7a7a; }
.eco-spinner     { font-size: 1rem; animation: spin 1.2s linear infinite; }
.eco-done        { font-size: 0.85rem; color: #22c55e; font-weight: 700; }

@keyframes spin { to { transform: rotate(360deg); } }

.eco-error {
  font-size: 0.8rem;
  color: #f87171;
  background: rgba(248,113,113,0.08);
  border: 1px solid rgba(248,113,113,0.2);
  border-radius: 6px;
  padding: 0.5rem 0.8rem;
}

/* Tabla */
.eco-table-wrap  { background: #0a0d14; border: 1px solid #1e2235; border-radius: 8px; overflow: hidden; }
.eco-table-scroll { overflow-x: auto; }

.eco-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.8rem;
  font-family: monospace;
}
.eco-table th {
  background: #0f1520;
  color: #00bfbf;
  padding: 0.5rem 0.8rem;
  text-align: left;
  font-weight: 600;
  border-bottom: 1px solid #1e2235;
  white-space: nowrap;
}
.eco-table td {
  padding: 0.4rem 0.8rem;
  color: #94a3b8;
  border-bottom: 1px solid #0f1520;
  white-space: nowrap;
}
.eco-table tr:last-child td { border-bottom: none; }
.eco-table tr:hover td { background: #0f1520; color: #cbd5e1; }

/* Imagen */
.eco-img-wrap { border: 1px solid #1e2235; border-radius: 8px; overflow: hidden; }
.eco-img      { width: 100%; display: block; }

/* Tabla title */
.eco-table-title {
  font-size: 0.75rem;
  font-weight: 600;
  color: #00bfbf;
  padding: 0.5rem 0.8rem 0;
  text-transform: uppercase;
  letter-spacing: 0.04em;
}
.pval-sig { color: #22c55e; font-weight: 700; }

/* PERMANOVA setup */
.perm-setup {
  background: #0a0f1a;
  border: 1px solid #1e3a3a;
  border-radius: 8px;
  padding: 1rem;
  display: flex;
  flex-direction: column;
  gap: 0.8rem;
}
.perm-hint { font-size: 0.78rem; color: #4b7a7a; }

.perm-samples {
  display: flex;
  flex-direction: column;
  gap: 0.4rem;
  max-height: 260px;
  overflow-y: auto;
}
.perm-row {
  display: flex;
  align-items: center;
  gap: 0.8rem;
}
.perm-sample {
  font-size: 0.8rem;
  font-family: monospace;
  color: #94a3b8;
  flex: 1;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
.perm-input {
  width: 130px;
  flex-shrink: 0;
  background: #0f1520;
  border: 1px solid #1e3a3a;
  border-radius: 5px;
  color: #e2e8f0;
  padding: 0.3rem 0.6rem;
  font-size: 0.8rem;
}
.perm-input:focus { outline: none; border-color: #00bfbf; }

.btn-perm-run {
  align-self: flex-start;
  background: #004d4d;
  border: 1px solid #00bfbf;
  color: #00ffff;
  border-radius: 6px;
  padding: 0.45rem 1.1rem;
  font-size: 0.85rem;
  font-weight: 600;
  cursor: pointer;
  transition: background 0.15s;
}
.btn-perm-run:hover:not(:disabled) { background: #006666; }
.btn-perm-run:disabled { opacity: 0.5; cursor: wait; }
</style>
