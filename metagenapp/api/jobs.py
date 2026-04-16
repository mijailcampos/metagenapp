"""
Job store y ejecución del pipeline en background.
Un solo pipeline corre a la vez (pipeline_lock).
"""
import threading
import uuid
import time
from dataclasses import dataclass, field
from typing import Optional, List

# ─── Lista canónica de steps (misma que run_pipeline.py) ───────────────────
PIPELINE_STEPS = [
    "00_Generate_files",
    "01_Assemble_contigs",
    "02_Filter_contigs",
    "03_Unique_contigs",
    "04_Alignment_vsearch",
    "09_Trim_vsearch",
    "10_Clustering_97",
    "11_Extract_Centroids",
    "12_MAFFT_Alignment",
    "14_MAFFT_Filter",
    "15_Filter_Columns",
    "16_Unique_MAFFT",
    "17_Precluster",
    "18_ChimeraDetection",
    "20_TaxonomicClassification",
    "21_RemoveLineages",
    "22_TaxonomicSummary",
    "23_Generate_OTUs",
    "24_Final_ASV_Table",
    "25_Assign_Taxonomy_to_ASVs",
    "26_ASV_Summary_Table",
]

RESULT_FILES = [
    "final_asv_table.tsv",
    "final_asv.taxonomy",
    "asv_resumen.tsv",
    "otu_table_0_03.tsv",
    "summary_tax_phylum.tsv",
    "final_clean.fasta",
    "final_clean.taxonomy",
]


# ─── Modelo de Job ──────────────────────────────────────────────────────────
@dataclass
class Job:
    job_id: str
    status: str          # queued | running | done | failed
    config: dict
    current_step: Optional[str] = None
    steps_done: List[str] = field(default_factory=list)
    logs: List[str] = field(default_factory=list)
    error: Optional[str] = None
    created_at: float = field(default_factory=time.time)
    finished_at: Optional[float] = None
    output_dir: Optional[str] = None
    qa_result: Optional[dict] = None   # resultado estructurado del QA


# ─── Store en memoria ───────────────────────────────────────────────────────
_jobs: dict = {}
_lock = threading.Lock()
_pipeline_lock = threading.Lock()   # solo un pipeline a la vez


# ─── API pública ────────────────────────────────────────────────────────────
def create_job(config: dict, output_dir: str) -> str:
    job_id = str(uuid.uuid4())[:8]
    job = Job(job_id=job_id, status="queued", config=config, output_dir=output_dir)
    with _lock:
        _jobs[job_id] = job
    return job_id


def get_job(job_id: str) -> Optional[Job]:
    return _jobs.get(job_id)


def list_jobs() -> List[Job]:
    with _lock:
        return list(_jobs.values())


def start_job(job_id: str) -> None:
    thread = threading.Thread(
        target=_run_pipeline_thread,
        args=(job_id,),
        daemon=True,
        name=f"pipeline-{job_id}",
    )
    thread.start()


def cancel_job(job_id: str) -> bool:
    """Marca el job como cancelado si aún está en cola."""
    job = _jobs.get(job_id)
    if job and job.status == "queued":
        with _lock:
            job.status = "cancelled"
            job.finished_at = time.time()
        return True
    return False


# ─── Worker ─────────────────────────────────────────────────────────────────
def _run_pipeline_thread(job_id: str) -> None:
    from metagenapp.pipeline.run_pipeline import run_pipeline
    from metagenapp.pipeline import step_tracker

    job = _jobs.get(job_id)
    if job is None:
        return

    def on_step_event(event: str, step_name: str, success: bool = True, data: dict = None) -> None:
        with _lock:
            if event == "start":
                job.current_step = step_name
                job.logs.append(f"▶  {step_name}")
            elif event == "end":
                if success:
                    job.steps_done.append(step_name)
                    job.logs.append(f"✅ {step_name}")
                else:
                    job.logs.append(f"❌ {step_name}")
                job.current_step = None
            elif event == "qa_result" and data:
                job.qa_result = data
            # Mantener los últimos 200 logs
            job.logs = job.logs[-200:]

    # Esperar turno (un pipeline a la vez)
    with _pipeline_lock:
        if job.status == "cancelled":
            return

        step_tracker.set_step_callback(on_step_event)
        job.status = "running"

        try:
            run_pipeline(**job.config)
            job.status = "done"
        except Exception as exc:
            job.status = "failed"
            job.error = str(exc)
            with _lock:
                job.logs.append(f"ERROR: {exc}")
        finally:
            job.finished_at = time.time()
            step_tracker.set_step_callback(None)
