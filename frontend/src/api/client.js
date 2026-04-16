const BASE = import.meta.env.VITE_API_URL || 'http://localhost:8000'

async function request(method, path, body = null) {
  const opts = { method, headers: {} }
  if (body && !(body instanceof FormData)) {
    opts.headers['Content-Type'] = 'application/json'
    opts.body = JSON.stringify(body)
  } else if (body instanceof FormData) {
    opts.body = body
  }
  const res = await fetch(`${BASE}${path}`, opts)
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: res.statusText }))
    throw new Error(err.detail || 'Error desconocido')
  }
  if (res.status === 204) return null
  return res.json()
}

export const api = {
  health: ()                        => request('GET',  '/health'),

  // Upload
  uploadFastq: (files)              => {
    const fd = new FormData()
    files.forEach(f => fd.append('files', f))
    return request('POST', '/upload', fd)
  },
  getUpload: (uploadId)             => request('GET',  `/upload/${uploadId}`),

  // Jobs
  createJob: (payload)              => request('POST', '/jobs', payload),
  listJobs: ()                      => request('GET',  '/jobs'),
  getJob: (jobId)                   => request('GET',  `/jobs/${jobId}`),
  getLogs: (jobId, last = 60)       => request('GET',  `/jobs/${jobId}/logs?last=${last}`),
  getQa: (jobId)                    => request('GET',  `/jobs/${jobId}/qa`),
  getResults: (jobId)               => request('GET',  `/jobs/${jobId}/results`),
  cancelJob: (jobId)                => request('DELETE', `/jobs/${jobId}`),

  downloadUrl: (jobId, filename)    => `${BASE}/jobs/${jobId}/download/${filename}`,

  // EcoLab
  runEcolab: (jobId, analysis, body = null) => request('POST', `/jobs/${jobId}/ecolab/${analysis}`, body),
  getEcolabSamples: (jobId)                 => request('GET',  `/jobs/${jobId}/ecolab/samples`),
  ecolabImageUrl: (jobId, filename)         => `${BASE}/jobs/${jobId}/ecolab/image/${filename}`,
}
