<script setup>
import { ref, onMounted } from 'vue'
import { useRouter } from 'vue-router'
import { api } from '../api/client.js'

const router = useRouter()
const jobs   = ref([])
const error  = ref(null)

function elapsed(s) {
  if (!s) return '—'
  const m = Math.floor(s / 60)
  const sec = Math.floor(s % 60)
  return m ? `${m}m ${sec}s` : `${sec}s`
}

function date(ts) {
  if (!ts) return '—'
  return new Date(ts * 1000).toLocaleString('es-ES', {
    month: 'short', day: 'numeric',
    hour: '2-digit', minute: '2-digit',
  })
}

async function load() {
  try {
    jobs.value = await api.listJobs()
  } catch (e) {
    error.value = e.message
  }
}

async function cancelJob(id) {
  try {
    await api.cancelJob(id)
    await load()
  } catch (e) {
    alert(e.message)
  }
}

onMounted(load)
</script>

<template>
  <div>
    <div class="page-header">
      <h1 class="page-title">Historial de corridas</h1>
      <button class="btn-refresh" @click="load">↻ Actualizar</button>
    </div>

    <p v-if="error" class="error-msg">{{ error }}</p>

    <p v-if="!jobs.length && !error" class="empty">
      Aún no hay corridas. <router-link to="/">Lanza la primera</router-link>.
    </p>

    <div class="jobs-list">
      <div v-for="job in jobs" :key="job.job_id" class="job-card">

        <div class="job-header">
          <div class="job-id-wrap">
            <span class="job-id">{{ job.job_id }}</span>
            <span class="status-badge" :class="job.status">{{ job.status }}</span>
          </div>
          <span class="job-date">{{ date(job.created_at) }}</span>
        </div>

        <!-- Barra de progreso -->
        <div class="progress-bar-wrap">
          <div
            class="progress-bar-fill"
            :class="job.status"
            :style="{ width: job.progress.percent + '%' }"
          ></div>
        </div>

        <div class="job-meta">
          <span>{{ job.progress.done }}/{{ job.progress.total }} steps</span>
          <span>⏱ {{ elapsed(job.elapsed_s) }}</span>
          <span class="pill">{{ job.config?.mode }}</span>
          <span class="pill">{{ job.config?.classifier }}</span>
          <span class="pill">{{ job.config?.marker }}</span>
        </div>

        <p v-if="job.status === 'failed'" class="job-error">
          ❌ {{ job.error }}
        </p>

        <div class="job-actions">
          <button
            v-if="job.status === 'running' || job.status === 'queued'"
            class="btn-action"
            @click="router.push(`/monitor/${job.job_id}`)"
          >Ver monitor</button>

          <button
            v-if="job.status === 'done'"
            class="btn-action primary"
            @click="router.push(`/results/${job.job_id}`)"
          >Ver resultados</button>

          <button
            v-if="job.status === 'queued'"
            class="btn-action danger"
            @click="cancelJob(job.job_id)"
          >Cancelar</button>
        </div>

      </div>
    </div>
  </div>
</template>

<style scoped>
.page-header { display: flex; align-items: center; justify-content: space-between; margin-bottom: 1.5rem; }
.page-title  { font-size: 1.6rem; font-weight: 700; color: #f1f5f9; }

.btn-refresh {
  background: #1e2235;
  border: 1px solid #2d3148;
  color: #94a3b8;
  border-radius: 6px;
  padding: 0.4rem 0.8rem;
  font-size: 0.85rem;
  cursor: pointer;
}
.btn-refresh:hover { border-color: #7c83fd; color: #7c83fd; }

.jobs-list { display: flex; flex-direction: column; gap: 1rem; }

.job-card {
  background: #1a1d27;
  border: 1px solid #2d3148;
  border-radius: 12px;
  padding: 1.2rem 1.5rem;
}

.job-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-bottom: 0.8rem;
}
.job-id-wrap { display: flex; align-items: center; gap: 0.6rem; }
.job-id   { font-family: monospace; font-size: 0.9rem; color: #7c83fd; font-weight: 700; }
.job-date { font-size: 0.78rem; color: #475569; }

.status-badge {
  font-size: 0.7rem;
  font-weight: 700;
  text-transform: uppercase;
  padding: 0.2rem 0.5rem;
  border-radius: 999px;
}
.status-badge.queued    { background: #1e3a5f; color: #60a5fa; }
.status-badge.running   { background: #1c3a1c; color: #4ade80; }
.status-badge.done      { background: #14532d; color: #86efac; }
.status-badge.failed    { background: #450a0a; color: #fca5a5; }
.status-badge.cancelled { background: #2d1f00; color: #fbbf24; }

.progress-bar-wrap {
  height: 5px;
  background: #12151f;
  border-radius: 999px;
  overflow: hidden;
  margin-bottom: 0.6rem;
}
.progress-bar-fill {
  height: 100%;
  border-radius: 999px;
  transition: width 0.4s ease;
  background: linear-gradient(90deg, #7c83fd, #a78bfa);
}
.progress-bar-fill.done   { background: #22c55e; }
.progress-bar-fill.failed { background: #ef4444; }

.job-meta {
  display: flex;
  align-items: center;
  gap: 0.6rem;
  font-size: 0.78rem;
  color: #64748b;
  flex-wrap: wrap;
}
.pill {
  background: #12151f;
  border: 1px solid #1e2235;
  border-radius: 4px;
  padding: 0.1rem 0.4rem;
  color: #94a3b8;
}

.job-error { font-size: 0.8rem; color: #f87171; margin-top: 0.5rem; }

.job-actions { display: flex; gap: 0.5rem; margin-top: 0.8rem; flex-wrap: wrap; }

.btn-action {
  background: #1e2235;
  border: 1px solid #2d3148;
  color: #94a3b8;
  border-radius: 6px;
  padding: 0.35rem 0.8rem;
  font-size: 0.82rem;
  cursor: pointer;
  transition: all 0.15s;
}
.btn-action:hover         { border-color: #7c83fd; color: #7c83fd; }
.btn-action.primary       { background: #7c83fd; border-color: #7c83fd; color: #fff; }
.btn-action.primary:hover { background: #6366f1; }
.btn-action.danger        { border-color: #ef4444; color: #f87171; }
.btn-action.danger:hover  { background: rgba(239,68,68,0.1); }

.error-msg { color: #f87171; font-size: 0.85rem; }
.empty     { color: #475569; }
.empty a   { color: #7c83fd; }
</style>
