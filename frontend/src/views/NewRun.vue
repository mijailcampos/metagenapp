<script setup>
import { ref, computed } from 'vue'
import { useRouter } from 'vue-router'
import { api } from '../api/client.js'

const router = useRouter()

// ── Upload ────────────────────────────────────────────────────────────────
const dragOver   = ref(false)
const files      = ref([])
const uploading  = ref(false)
const uploadDone = ref(false)
const uploadId   = ref(null)
const uploadErr  = ref(null)

function addFiles(newFiles) {
  const valid = Array.from(newFiles).filter(f => f.name.match(/\.(fastq|fq)(\.gz)?$/i))
  const existing = new Set(files.value.map(f => f.name))
  const toAdd = valid.filter(f => !existing.has(f.name))
  files.value = [...files.value, ...toAdd]
}

function onDrop(e) {
  dragOver.value = false
  addFiles(e.dataTransfer.files)
}

function onFileInput(e) {
  addFiles(e.target.files)
  e.target.value = ''   // permite volver a seleccionar los mismos archivos
}

function removeFile(i) {
  files.value = files.value.filter((_, idx) => idx !== i)
  if (files.value.length === 0) { uploadDone.value = false; uploadId.value = null }
}

function humanSize(bytes) {
  if (bytes < 1024)      return bytes + ' B'
  if (bytes < 1024 ** 2) return (bytes / 1024).toFixed(1) + ' KB'
  return (bytes / 1024 ** 2).toFixed(1) + ' MB'
}

async function doUpload() {
  uploadErr.value = null
  uploading.value = true
  try {
    const res = await api.uploadFastq(files.value)
    uploadId.value   = res.upload_id
    uploadDone.value = true
  } catch (e) {
    uploadErr.value = e.message
  } finally {
    uploading.value = false
  }
}

function resetUpload() {
  files.value      = []
  uploadDone.value = false
  uploadId.value   = null
  uploadErr.value  = null
  qaStatus.value   = 'idle'
  qaLogs.value     = []
  qaJobId.value    = null
}

// ── QA ────────────────────────────────────────────────────────────────────
const qaStatus = ref('idle')   // idle | running | done | failed
const qaLogs   = ref([])
const qaJobId  = ref(null)
const qaResult = ref(null)
let   qaTimer  = null

async function runQA() {
  qaStatus.value = 'running'
  qaLogs.value   = []
  qaResult.value = null
  try {
    const res = await api.createJob({
      upload_id:  uploadId.value,
      mode:       'qa',
      classifier: 'flat',
      marker:     params.value.marker,
      threads:    params.value.threads,
      min_length: params.value.min_length,
      max_length: params.value.max_length,
      max_ambigs: params.value.max_ambigs,
      max_poly:   params.value.max_poly,
      extract_centroids: 'full',
    })
    qaJobId.value = res.job_id
    qaTimer = setInterval(pollQA, 2000)
  } catch (e) {
    qaStatus.value = 'failed'
    qaLogs.value   = ['ERROR: ' + e.message]
  }
}

async function pollQA() {
  try {
    const [job, logRes] = await Promise.all([
      api.getJob(qaJobId.value),
      api.getLogs(qaJobId.value, 80),
    ])
    qaLogs.value = logRes.logs || []
    if (job.status === 'done' || job.status === 'failed') {
      clearInterval(qaTimer)
      qaStatus.value = job.status
      if (job.status === 'failed') {
        qaLogs.value.push('ERROR: ' + job.error)
      } else {
        // Obtener resultado estructurado
        try {
          qaResult.value = await api.getQa(qaJobId.value)
          // Pre-rellenar params con valores sugeridos del QA
          if (qaResult.value?.stats?.median_bp) {
            const med = qaResult.value.stats.median_bp
            params.value.min_length = Math.max(100, Math.round(med * 0.85))
            params.value.max_length = Math.min(1000, Math.round(med * 1.15))
          }
        } catch (_) {}
      }
    }
  } catch (e) {
    clearInterval(qaTimer)
    qaStatus.value = 'failed'
  }
}

// ── Parámetros ────────────────────────────────────────────────────────────
const params = ref({
  mode:              'student',
  classifier:        'flat',
  marker:            '16S',
  model_type:        'general',
  threads:           8,
  min_length:        250,
  max_length:        600,
  max_ambigs:        0,
  max_poly:          8,
  extract_centroids: 'full',
})

// ── Submit ────────────────────────────────────────────────────────────────
const submitting = ref(false)
const submitErr  = ref(null)
const canSubmit  = computed(() => uploadDone.value && !submitting.value)

async function launch() {
  submitErr.value = null
  submitting.value = true
  try {
    const res = await api.createJob({ upload_id: uploadId.value, ...params.value })
    router.push(`/monitor/${res.job_id}`)
  } catch (e) {
    submitErr.value = e.message
  } finally {
    submitting.value = false
  }
}
</script>

<template>
  <div>
    <h1 class="page-title">Nueva corrida</h1>
    <p class="page-sub">Sube tus archivos FASTQ y configura el pipeline</p>

    <!-- ── 1. Upload ──────────────────────────────────────────────────── -->
    <section class="card">
      <h2 class="card-title">1 — Archivos FASTQ</h2>

      <div
        class="dropzone"
        :class="{ active: dragOver, done: uploadDone }"
        @dragover.prevent="dragOver = true"
        @dragleave="dragOver = false"
        @drop.prevent="onDrop"
      >
        <template v-if="!uploadDone">
          <span class="drop-icon">📂</span>
          <p>Arrastra tus archivos FASTQ aquí</p>
          <p class="drop-hint">(.fastq, .fq, .fastq.gz) — puedes soltar varias veces para agregar más</p>
          <label class="btn-secondary">
            + Agregar archivos
            <input
              type="file"
              multiple
              accept=".fastq,.fq,.gz"
              hidden
              @change="onFileInput"
            />
          </label>
        </template>
        <template v-else>
          <span class="drop-icon">✅</span>
          <p>{{ files.length }} archivo(s) subidos al servidor</p>
          <p class="drop-hint">upload_id: <code>{{ uploadId }}</code></p>
          <button class="btn-secondary" style="margin-top:0.5rem" @click="resetUpload">
            ✕ Cambiar archivos
          </button>
        </template>
      </div>

      <!-- Lista de archivos pendientes -->
      <template v-if="files.length && !uploadDone">
        <p class="files-count">{{ files.length }} archivo(s) seleccionado(s)</p>
        <ul class="file-list">
          <li v-for="(f, i) in files" :key="i" class="file-item">
            <span class="file-name">{{ f.name }}</span>
            <span class="file-size">{{ humanSize(f.size) }}</span>
            <button class="btn-remove" @click="removeFile(i)">✕</button>
          </li>
        </ul>
        <p v-if="uploadErr" class="error-msg">{{ uploadErr }}</p>
        <button class="btn-primary" :disabled="uploading" @click="doUpload">
          {{ uploading ? 'Subiendo…' : `Subir ${files.length} archivo(s)` }}
        </button>
      </template>
    </section>

    <!-- ── 2. QA ───────────────────────────────────────────────────────── -->
    <section class="card" v-if="uploadDone">
      <h2 class="card-title">2 — Análisis de calidad (QA) <span class="optional">recomendado</span></h2>
      <p class="card-desc">
        Ensambla y filtra las lecturas con los parámetros actuales para ver cuántas pasan.
        Úsalo para ajustar <strong>longitud mínima / máxima</strong> antes de correr el pipeline completo.
      </p>

      <button
        class="btn-qa"
        :disabled="qaStatus === 'running'"
        @click="runQA"
      >
        {{ qaStatus === 'running' ? '⏳ Analizando…' : '🔍 Correr QA' }}
      </button>

      <!-- Resultados QA -->
      <div v-if="qaStatus !== 'idle'" class="qa-results">
        <div class="qa-status-row">
          <span class="status-badge" :class="qaStatus">
            {{ qaStatus === 'running' ? 'corriendo' : qaStatus }}
          </span>
          <span v-if="qaJobId" class="qa-job-id">job: {{ qaJobId }}</span>
        </div>

        <!-- Resultado estructurado -->
        <div v-if="qaResult" class="qa-panel">

          <!-- ── Sección 1: Ensamblado por muestra ── -->
          <div class="qa-section" v-if="qaResult.assembly?.samples?.length">
            <p class="qa-section-title">Ensamblado de contigs</p>
            <table class="qa-table">
              <thead>
                <tr>
                  <th>Muestra</th>
                  <th class="num">Pares</th>
                  <th class="num">Mergeados</th>
                  <th class="num">%</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="s in qaResult.assembly.samples" :key="s.sample">
                  <td class="mono">{{ s.sample }}</td>
                  <td class="num">{{ s.pairs.toLocaleString() }}</td>
                  <td class="num">{{ s.merged.toLocaleString() }}</td>
                  <td class="num" :class="s.pct >= 70 ? 'ok' : s.pct >= 50 ? 'warn' : 'bad'">{{ s.pct }}%</td>
                </tr>
              </tbody>
              <tfoot>
                <tr class="total-row">
                  <td>Total</td>
                  <td class="num">{{ qaResult.assembly.total_pairs.toLocaleString() }}</td>
                  <td class="num">{{ qaResult.assembly.total_merged.toLocaleString() }}</td>
                  <td class="num" :class="qaResult.assembly.overall_pct >= 70 ? 'ok' : 'warn'">{{ qaResult.assembly.overall_pct }}%</td>
                </tr>
              </tfoot>
            </table>
            <p class="qa-section-note">
              {{ qaResult.stats.unique_contigs?.toLocaleString() }} contigs únicos &nbsp;·&nbsp;
              {{ qaResult.retained?.toLocaleString() }} pasan el filtro
            </p>
          </div>

          <!-- ── Sección 2: Stats de contigs ── -->
          <div class="qa-section" v-if="qaResult.stats_table && Object.keys(qaResult.stats_table).length">
            <p class="qa-section-title">Estadísticas de contigs (tras filtrado)</p>
            <table class="qa-table">
              <thead>
                <tr>
                  <th>Métrica</th>
                  <th class="num">Start</th>
                  <th class="num">End</th>
                  <th class="num nbases">NBases</th>
                  <th class="num">Ambigs</th>
                  <th class="num">Polymer</th>
                  <th class="num">NumSeqs</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="(vals, stat) in qaResult.stats_table" :key="stat">
                  <td class="metric-label">{{ stat }}</td>
                  <td class="num">{{ vals.Start }}</td>
                  <td class="num">{{ vals.End }}</td>
                  <td class="num nbases">{{ vals.NBases }}</td>
                  <td class="num" :class="Number(vals.Ambigs) > 0 ? 'warn' : 'dim'">{{ vals.Ambigs }}</td>
                  <td class="num" :class="Number(vals.Polymer) >= 8 ? 'bad' : 'dim'">{{ vals.Polymer }}</td>
                  <td class="num">{{ vals.NumSeqs }}</td>
                </tr>
              </tbody>
            </table>
            <p class="qa-section-note dim">NBases = longitud de contig &nbsp;·&nbsp; Ambigs = bases ambiguas &nbsp;·&nbsp; Polymer ≥ 8 indica homopolímero largo</p>
          </div>

          <!-- ── Sección 3: Marcadores inferidos ── -->
          <div class="qa-section">
            <p class="qa-section-title">
              Marcador inferido &nbsp;
              <span class="median-badge">mediana {{ qaResult.stats.median_bp }} bp</span>
            </p>
            <table class="qa-table">
              <thead>
                <tr>
                  <th>Región</th>
                  <th>Primers</th>
                  <th class="num">Rango típico</th>
                  <th></th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="m in qaResult.markers" :key="m.region" :class="{ 'match-row': m.match }">
                  <td class="metric-label">{{ m.region }}</td>
                  <td class="mono dim">{{ m.primers }}</td>
                  <td class="num">{{ m.min_bp }}–{{ m.max_bp }} bp</td>
                  <td class="match-tag">{{ m.match ? '← esta muestra' : '' }}</td>
                </tr>
              </tbody>
            </table>
            <p v-if="qaResult.markers.filter(m=>m.match).length > 1" class="qa-warn">
              ⚠ Múltiples regiones coinciden — confirma con tus primers para desambiguar.
            </p>
            <p v-else-if="!qaResult.markers.some(m=>m.match)" class="qa-warn bad">
              ✗ Longitud mediana fuera de rango para cualquier marcador 16S conocido.
            </p>
          </div>

          <!-- Parámetros de filtrado → van al pipeline completo -->
          <div class="qa-cutoffs">
            <p class="qa-cutoffs-title">Confirmar parámetros de filtrado</p>
            <p class="qa-cutoffs-hint">
              Valores precargados desde el QA (±15% de la mediana). Ajústalos si es necesario —
              se usarán directamente en el pipeline completo.
            </p>
            <div class="cutoff-row">
              <div class="param-inline">
                <label>Marcador</label>
                <select v-model="params.marker">
                  <option value="16S">16S</option>
                  <option value="18S">18S</option>
                </select>
              </div>
              <div class="param-inline">
                <label>Long. mínima (bp)</label>
                <input type="number" v-model.number="params.min_length" min="100" max="1000" />
              </div>
              <div class="param-inline">
                <label>Long. máxima (bp)</label>
                <input type="number" v-model.number="params.max_length" min="100" max="1000" />
              </div>
              <div class="param-inline">
                <label>Max ambiguas</label>
                <input type="number" v-model.number="params.max_ambigs" min="0" max="10" />
              </div>
              <div class="param-inline">
                <label>Max homopolímero</label>
                <input type="number" v-model.number="params.max_poly" min="4" max="20" />
              </div>
            </div>
            <p class="qa-cutoffs-arrow">↓ Continúa en los parámetros del pipeline</p>
          </div>
        </div>

        <!-- Log durante la ejecución -->
        <div class="log-box" v-if="qaStatus === 'running' || !qaResult">
          <p
            v-for="(line, i) in qaLogs"
            :key="i"
            class="log-line"
            :class="{
              'log-ok':  line.startsWith('✅'),
              'log-err': line.startsWith('❌') || line.startsWith('ERROR'),
              'log-run': line.startsWith('▶'),
            }"
          >{{ line }}</p>
          <p v-if="!qaLogs.length" class="log-empty">Esperando logs…</p>
        </div>
      </div>
    </section>

    <!-- ── 3. Parámetros completos ────────────────────────────────────── -->
    <section class="card" v-if="qaResult">
      <h2 class="card-title">{{ uploadDone ? '3' : '2' }} — Parámetros del pipeline</h2>
      <div class="params-grid">

        <div class="param">
          <label>Modo</label>
          <select v-model="params.mode">
            <option value="student">student — rápido (10k centroids)</option>
            <option value="premium">premium — completo</option>
            <option value="turbo">turbo — máximo CPU</option>
            <option value="ref">ref — calidad publicación</option>
          </select>
        </div>

        <div class="param">
          <label>Clasificador</label>
          <select v-model="params.classifier">
            <option value="flat">flat — Naive Bayes básico</option>
            <option value="naive-v2">naive-v2 — Wang bootstrap ⭐</option>
            <option value="kraken-lite">kraken-lite — k-mer + LCA</option>
            <option value="pro-engine">pro-engine — híbrido alta res.</option>
          </select>
        </div>

        <div class="param">
          <label>Marcador</label>
          <select v-model="params.marker">
            <option value="16S">16S rRNA</option>
            <option value="18S">18S rRNA</option>
          </select>
        </div>

        <div class="param">
          <label>Tipo de modelo</label>
          <select v-model="params.model_type">
            <option value="general">general</option>
            <option value="oral">oral</option>
            <option value="gut">gut</option>
            <option value="skin">skin</option>
            <option value="env">env</option>
          </select>
        </div>

        <div class="param">
          <label>Threads</label>
          <input type="number" v-model.number="params.threads" min="1" max="64" />
        </div>

        <div class="param">
          <label>Long. mínima (bp)</label>
          <input type="number" v-model.number="params.min_length" min="100" max="1000" />
        </div>

        <div class="param">
          <label>Long. máxima (bp)</label>
          <input type="number" v-model.number="params.max_length" min="100" max="1000" />
        </div>

        <div class="param">
          <label>Max ambiguas</label>
          <input type="number" v-model.number="params.max_ambigs" min="0" max="10" />
        </div>

        <div class="param">
          <label>Max homopolímero</label>
          <input type="number" v-model.number="params.max_poly" min="4" max="20" />
        </div>

        <div class="param">
          <label>Centroids</label>
          <select v-model="params.extract_centroids">
            <option value="test">test (1k)</option>
            <option value="student">student (10k)</option>
            <option value="full">full (todos)</option>
          </select>
        </div>

      </div>
    </section>

    <!-- ── Lanzar ─────────────────────────────────────────────────────── -->
    <template v-if="qaResult">
      <p v-if="submitErr" class="error-msg">{{ submitErr }}</p>
      <button class="btn-primary btn-large" :disabled="!canSubmit" @click="launch">
        {{ submitting ? 'Lanzando…' : 'Lanzar pipeline completo' }}
      </button>
    </template>
  </div>
</template>

<style scoped>
.page-title { font-size: 1.6rem; font-weight: 700; color: #f1f5f9; }
.page-sub   { color: #64748b; margin-top: 0.3rem; margin-bottom: 1.5rem; }

.card {
  background: #1a1d27;
  border: 1px solid #2d3148;
  border-radius: 12px;
  padding: 1.5rem;
  margin-bottom: 1.5rem;
}

.card-title {
  font-size: 1rem;
  font-weight: 600;
  color: #94a3b8;
  margin-bottom: 0.5rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  display: flex;
  align-items: center;
  gap: 0.6rem;
}

.optional {
  font-size: 0.7rem;
  background: rgba(124,131,253,0.15);
  color: #7c83fd;
  padding: 0.15rem 0.5rem;
  border-radius: 4px;
  text-transform: none;
  letter-spacing: 0;
}

.card-desc { font-size: 0.85rem; color: #64748b; margin-bottom: 1rem; }

/* Dropzone */
.dropzone {
  border: 2px dashed #2d3148;
  border-radius: 10px;
  padding: 2rem;
  text-align: center;
  transition: all 0.2s;
  color: #64748b;
}
.dropzone.active { border-color: #7c83fd; background: rgba(124,131,253,0.05); }
.dropzone.done   { border-color: #22c55e; background: rgba(34,197,94,0.05); color: #4ade80; }
.drop-icon  { font-size: 2.5rem; display: block; margin-bottom: 0.5rem; }
.drop-hint  { font-size: 0.8rem; color: #475569; margin: 0.3rem 0 0.8rem; }

.files-count { font-size: 0.82rem; color: #64748b; margin: 0.8rem 0 0.4rem; }

/* File list */
.file-list  { list-style: none; margin-bottom: 0.8rem; display: flex; flex-direction: column; gap: 0.3rem; max-height: 220px; overflow-y: auto; }
.file-item  { display: flex; align-items: center; gap: 0.5rem; background: #12151f; padding: 0.4rem 0.8rem; border-radius: 6px; }
.file-name  { flex: 1; font-size: 0.82rem; color: #cbd5e1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-family: monospace; }
.file-size  { font-size: 0.72rem; color: #475569; white-space: nowrap; }
.btn-remove { background: none; border: none; color: #ef4444; cursor: pointer; font-size: 0.8rem; padding: 0 0.2rem; }

/* QA */
.qa-params-row {
  display: flex;
  flex-wrap: wrap;
  gap: 0.8rem;
  margin-bottom: 1rem;
}
.param-inline { display: flex; flex-direction: column; gap: 0.25rem; }
.param-inline label { font-size: 0.75rem; color: #64748b; }
.param-inline select,
.param-inline input {
  background: #12151f;
  border: 1px solid #2d3148;
  border-radius: 6px;
  color: #e2e8f0;
  padding: 0.35rem 0.5rem;
  font-size: 0.82rem;
  width: 130px;
}
.param-inline select:focus,
.param-inline input:focus { outline: none; border-color: #7c83fd; }

.btn-qa {
  background: #1e2235;
  border: 1px solid #475569;
  color: #94a3b8;
  border-radius: 8px;
  padding: 0.5rem 1.2rem;
  font-size: 0.9rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.15s;
}
.btn-qa:hover:not(:disabled) { border-color: #7c83fd; color: #7c83fd; }
.btn-qa:disabled { opacity: 0.5; cursor: not-allowed; }

.qa-results { margin-top: 1rem; }
.qa-status-row { display: flex; align-items: center; gap: 0.6rem; margin-bottom: 0.5rem; }
.qa-job-id { font-size: 0.75rem; color: #475569; font-family: monospace; }
.qa-done-hint { font-size: 0.82rem; color: #4ade80; margin-top: 0.6rem; }

.status-badge {
  font-size: 0.7rem; font-weight: 700; text-transform: uppercase;
  padding: 0.2rem 0.5rem; border-radius: 999px;
}
.status-badge.running { background: #1c3a1c; color: #4ade80; }
.status-badge.done    { background: #14532d; color: #86efac; }
.status-badge.failed  { background: #450a0a; color: #fca5a5; }

/* Log */
.log-box {
  background: #0f1117;
  border: 1px solid #1e2235;
  border-radius: 8px;
  padding: 0.8rem 1rem;
  max-height: 200px;
  overflow-y: auto;
  font-family: 'JetBrains Mono', 'Fira Code', monospace;
  font-size: 0.76rem;
}
.log-line  { color: #64748b; line-height: 1.7; }
.log-ok    { color: #4ade80; font-weight: 600; }
.log-err   { color: #f87171; }
.log-run   { color: #7c83fd; }
.log-empty { color: #334155; font-style: italic; }

/* Params grid */
.params-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 1rem;
}
.param label { display: block; font-size: 0.78rem; color: #64748b; margin-bottom: 0.3rem; }
.param select,
.param input {
  width: 100%;
  background: #12151f;
  border: 1px solid #2d3148;
  border-radius: 6px;
  color: #e2e8f0;
  padding: 0.45rem 0.6rem;
  font-size: 0.85rem;
}
.param select:focus,
.param input:focus { outline: none; border-color: #7c83fd; }

/* Buttons */
.btn-primary {
  background: #7c83fd;
  color: #fff;
  border: none;
  border-radius: 8px;
  padding: 0.55rem 1.2rem;
  font-size: 0.9rem;
  font-weight: 600;
  cursor: pointer;
  margin-top: 1rem;
  transition: background 0.15s;
}
.btn-primary:hover:not(:disabled) { background: #6366f1; }
.btn-primary:disabled { opacity: 0.4; cursor: not-allowed; }
.btn-primary.btn-large { padding: 0.8rem 2rem; font-size: 1rem; }

.btn-secondary {
  display: inline-block;
  background: #1e2235;
  border: 1px solid #2d3148;
  color: #94a3b8;
  border-radius: 6px;
  padding: 0.4rem 0.9rem;
  font-size: 0.85rem;
  cursor: pointer;
  margin-top: 0.5rem;
}
.btn-secondary:hover { border-color: #7c83fd; color: #7c83fd; }

.error-msg { color: #f87171; font-size: 0.85rem; margin-top: 0.5rem; }
.hint      { color: #475569; font-size: 0.8rem; margin-top: 0.5rem; }

/* QA panel */
.qa-panel { margin-top: 1rem; display: flex; flex-direction: column; gap: 1rem; }

/* QA Sections */
.qa-section {
  background: #12151f;
  border: 1px solid #1e2235;
  border-radius: 8px;
  padding: 1rem;
}
.qa-section-title {
  font-size: 0.78rem;
  font-weight: 700;
  color: #64748b;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  margin-bottom: 0.7rem;
  display: flex;
  align-items: center;
  gap: 0.5rem;
}
.qa-section-note { font-size: 0.75rem; color: #475569; margin-top: 0.5rem; }
.median-badge {
  font-size: 0.72rem;
  background: rgba(124,131,253,0.15);
  color: #7c83fd;
  padding: 0.1rem 0.5rem;
  border-radius: 4px;
  font-weight: 600;
  text-transform: none;
  letter-spacing: 0;
}
.qa-warn { font-size: 0.8rem; color: #facc15; margin-top: 0.5rem; }
.qa-warn.bad { color: #f87171; }

/* QA Table */
.qa-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.8rem;
}
.qa-table th {
  color: #64748b;
  font-weight: 600;
  text-align: left;
  padding: 0.25rem 0.5rem;
  border-bottom: 1px solid #1e2235;
  white-space: nowrap;
}
.qa-table td {
  padding: 0.2rem 0.5rem;
  color: #94a3b8;
  border-bottom: 1px solid #0f1117;
}
.qa-table tfoot td { border-top: 1px solid #2d3148; border-bottom: none; }
.qa-table .num { text-align: right; font-family: monospace; }
.qa-table .mono { font-family: monospace; font-size: 0.75rem; }
.qa-table .dim { color: #334155; }
.qa-table .ok   { color: #4ade80; }
.qa-table .warn { color: #facc15; }
.qa-table .bad  { color: #f87171; }
.qa-table .nbases { color: #4ade80; font-weight: 600; }
.qa-table .metric-label { font-weight: 600; color: #cbd5e1; }
.qa-table .total-row td { font-weight: 700; color: #e2e8f0; }
.qa-table .match-row td { color: #4ade80; background: rgba(34,197,94,0.05); }
.qa-table .match-tag { color: #4ade80; font-size: 0.75rem; font-weight: 600; }

.qa-cutoffs {
  background: rgba(124,131,253,0.06);
  border: 1px solid rgba(124,131,253,0.2);
  border-radius: 8px;
  padding: 1rem;
}
.qa-cutoffs-title  { font-size: 0.85rem; font-weight: 600; color: #a5b4fc; margin-bottom: 0.3rem; }
.qa-cutoffs-hint   { font-size: 0.78rem; color: #64748b; margin-bottom: 0.8rem; }
.cutoff-row        { display: flex; gap: 1rem; align-items: flex-end; flex-wrap: wrap; }
.qa-cutoffs-arrow  { font-size: 0.8rem; color: #7c83fd; margin-top: 0.8rem; font-weight: 600; }
</style>
