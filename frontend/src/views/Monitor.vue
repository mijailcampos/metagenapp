<script setup>
import { ref, onMounted, onUnmounted, computed } from 'vue'
import { useRouter } from 'vue-router'
import { api } from '../api/client.js'

const props  = defineProps({ jobId: String })
const router = useRouter()

const job    = ref(null)
const logs   = ref([])
const error  = ref(null)
let   timer  = null

const STEPS = [
  '00_Generate_files','01_Assemble_contigs','02_Filter_contigs','03_Unique_contigs',
  '04_Alignment_vsearch','09_Trim_vsearch','10_Clustering_97','11_Extract_Centroids',
  '12_MAFFT_Alignment','14_MAFFT_Filter','15_Filter_Columns','16_Unique_MAFFT',
  '17_Precluster','18_ChimeraDetection','20_TaxonomicClassification','21_RemoveLineages',
  '22_TaxonomicSummary','23_Generate_OTUs','24_Final_ASV_Table',
  '25_Assign_Taxonomy_to_ASVs','26_ASV_Summary_Table',
]

const stepsDone = computed(() => new Set(job.value?.progress?.done
  ? STEPS.slice(0, job.value.progress.done)
  : []))

function stepState(s) {
  if (!job.value) return 'pending'
  if (job.value.current_step === s) return 'running'
  if (stepsDone.value.has(s))       return 'done'
  return 'pending'
}

function elapsed(seconds) {
  if (!seconds) return '—'
  const m = Math.floor(seconds / 60)
  const s = Math.floor(seconds % 60)
  return m ? `${m}m ${s}s` : `${s}s`
}

async function poll() {
  try {
    const [j, l] = await Promise.all([
      api.getJob(props.jobId),
      api.getLogs(props.jobId, 40),
    ])
    job.value  = j
    logs.value = l.logs || []

    if (j.status === 'done') {
      clearInterval(timer)
      setTimeout(() => router.push(`/results/${props.jobId}`), 1500)
    } else if (j.status === 'failed') {
      clearInterval(timer)
    }
  } catch (e) {
    error.value = e.message
    clearInterval(timer)
  }
}

onMounted(() => { poll(); timer = setInterval(poll, 2000) })
onUnmounted(() => clearInterval(timer))
</script>

<template>
  <div>
    <div class="page-header">
      <h1 class="page-title">Monitor de corrida</h1>
      <span class="job-id">ID: {{ jobId }}</span>
    </div>

    <p v-if="error" class="error-msg">{{ error }}</p>

    <div v-if="job" class="monitor-layout">

      <!-- ── Barra progreso ──────────────────────────────────────────── -->
      <section class="card">
        <div class="progress-header">
          <span class="status-badge" :class="job.status">{{ job.status }}</span>
          <span class="elapsed">⏱ {{ elapsed(job.elapsed_s) }}</span>
        </div>

        <div class="progress-bar-wrap">
          <div class="progress-bar-fill" :style="{ width: job.progress.percent + '%' }"></div>
        </div>

        <p class="progress-label">
          {{ job.progress.done }} / {{ job.progress.total }} steps
          ({{ job.progress.percent }}%)
        </p>

        <p v-if="job.current_step" class="current-step">
          ▶ {{ job.current_step.replace(/_/g, ' ') }}
        </p>
        <p v-if="job.status === 'done'" class="done-msg">
          ✅ Pipeline completado — redirigiendo a resultados…
        </p>
        <p v-if="job.status === 'failed'" class="error-msg">
          ❌ Error: {{ job.error }}
        </p>
      </section>

      <!-- ── Steps ──────────────────────────────────────────────────── -->
      <section class="card">
        <h2 class="card-title">Steps</h2>
        <ul class="steps-list">
          <li
            v-for="s in STEPS"
            :key="s"
            class="step-item"
            :class="stepState(s)"
          >
            <span class="step-icon">
              {{ stepState(s) === 'done' ? '✅' : stepState(s) === 'running' ? '⚙️' : '○' }}
            </span>
            <span class="step-name">{{ s.replace(/_/g, ' ') }}</span>
          </li>
        </ul>
      </section>

      <!-- ── Log ────────────────────────────────────────────────────── -->
      <section class="card">
        <h2 class="card-title">Log en vivo</h2>
        <div class="log-box">
          <p v-for="(line, i) in logs" :key="i" class="log-line" :class="{
            'log-ok':  line.startsWith('✅'),
            'log-err': line.startsWith('❌') || line.startsWith('ERROR'),
            'log-run': line.startsWith('▶'),
          }">{{ line }}</p>
          <p v-if="!logs.length" class="log-empty">Esperando logs…</p>
        </div>
      </section>

    </div>

    <div v-else-if="!error" class="loading">Cargando…</div>
  </div>
</template>

<style scoped>
.page-header { display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1.5rem; }
.page-title  { font-size: 1.6rem; font-weight: 700; color: #f1f5f9; }
.job-id      { font-size: 0.8rem; color: #475569; font-family: monospace; }

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
  margin-bottom: 1rem;
}

/* Status badge */
.status-badge {
  font-size: 0.75rem;
  font-weight: 700;
  text-transform: uppercase;
  padding: 0.25rem 0.6rem;
  border-radius: 999px;
  letter-spacing: 0.06em;
}
.status-badge.queued  { background: #1e3a5f; color: #60a5fa; }
.status-badge.running { background: #1c3a1c; color: #4ade80; }
.status-badge.done    { background: #14532d; color: #86efac; }
.status-badge.failed  { background: #450a0a; color: #fca5a5; }

.progress-header { display: flex; align-items: center; justify-content: space-between; margin-bottom: 0.8rem; }
.elapsed { font-size: 0.85rem; color: #64748b; }

.progress-bar-wrap {
  height: 8px;
  background: #12151f;
  border-radius: 999px;
  overflow: hidden;
}
.progress-bar-fill {
  height: 100%;
  background: linear-gradient(90deg, #7c83fd, #a78bfa);
  border-radius: 999px;
  transition: width 0.4s ease;
}

.progress-label { font-size: 0.8rem; color: #64748b; margin-top: 0.5rem; }
.current-step   { font-size: 0.9rem; color: #7c83fd; margin-top: 0.5rem; font-weight: 500; }
.done-msg       { color: #4ade80; margin-top: 0.5rem; font-size: 0.9rem; }

/* Steps */
.steps-list { list-style: none; display: flex; flex-direction: column; gap: 0.3rem; }
.step-item  { display: flex; align-items: center; gap: 0.5rem; padding: 0.3rem 0; font-size: 0.85rem; color: #475569; }
.step-item.done    { color: #4ade80; }
.step-item.running { color: #7c83fd; font-weight: 600; }
.step-icon  { width: 1.4rem; text-align: center; flex-shrink: 0; }

/* Log */
.log-box {
  background: #0f1117;
  border: 1px solid #1e2235;
  border-radius: 8px;
  padding: 1rem;
  max-height: 260px;
  overflow-y: auto;
  font-family: 'JetBrains Mono', 'Fira Code', monospace;
  font-size: 0.78rem;
}
.log-line   { color: #64748b; line-height: 1.6; }
.log-ok     { color: #4ade80; }
.log-err    { color: #f87171; }
.log-run    { color: #7c83fd; }
.log-empty  { color: #334155; font-style: italic; }

.error-msg { color: #f87171; font-size: 0.85rem; margin-top: 0.5rem; }
.loading   { color: #475569; padding: 2rem; text-align: center; }
</style>
