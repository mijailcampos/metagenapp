"""
MetagenApp — API REST (FastAPI)
Lanzar con:  uvicorn metagenapp.api.main:app --reload --port 8000
"""
import csv
import shutil
import subprocess
import tempfile
import time
import uuid
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, field_validator, model_validator

from metagenapp.metagen_config import OUTPUT_DIR, INPUT_DIR
from metagenapp.api import jobs as job_store

# ─── App ────────────────────────────────────────────────────────────────────
app = FastAPI(
    title="MetagenApp API",
    version="0.1",
    description="API REST para el pipeline de metabarcoding MetagenApp.",
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://localhost:3000",
        "http://192.168.100.200:5173",
        "https://bioagens.com",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Router con prefijo /api para todos los endpoints
api = APIRouter(prefix="/api")


# ─── Esquemas ────────────────────────────────────────────────────────────────
class PipelineRequest(BaseModel):
    # Fuente de datos: uno de los dos es obligatorio
    input_dir: Optional[str] = None   # ruta en el servidor
    upload_id: Optional[str] = None   # ID de un upload previo desde el navegador

    outdir: Optional[str] = None
    mode: str = "student"
    classifier: str = "flat"
    marker: str = "16S"
    model_type: str = "general"
    threads: int = 8
    min_length: int = 250
    max_length: int = 600
    max_ambigs: int = 0
    max_poly: int = 8
    extract_centroids: str = "full"
    from_step: Optional[str] = None

    @model_validator(mode="after")
    def check_input_source(self):
        if not self.input_dir and not self.upload_id:
            raise ValueError("Debes proporcionar 'input_dir' o 'upload_id'")
        return self

    @field_validator("mode")
    @classmethod
    def validate_mode(cls, v):
        valid = {"student", "premium", "turbo", "ref", "qa"}
        if v not in valid:
            raise ValueError(f"mode debe ser uno de: {valid}")
        return v

    @field_validator("classifier")
    @classmethod
    def validate_classifier(cls, v):
        valid = {"flat", "naive-v2", "kraken-lite", "pro-engine"}
        if v not in valid:
            raise ValueError(f"classifier debe ser uno de: {valid}")
        return v

    @field_validator("marker")
    @classmethod
    def validate_marker(cls, v):
        if v not in {"16S", "18S"}:
            raise ValueError("marker debe ser '16S' o '18S'")
        return v


# ─── Helpers ─────────────────────────────────────────────────────────────────
def _job_to_dict(job: job_store.Job) -> dict:
    total = len(job_store.PIPELINE_STEPS)
    done = len(job.steps_done)
    elapsed = None
    if job.created_at and job.finished_at:
        elapsed = round(job.finished_at - job.created_at, 1)
    elif job.status == "running":
        elapsed = round(time.time() - job.created_at, 1)

    return {
        "job_id": job.job_id,
        "status": job.status,
        "current_step": job.current_step,
        "progress": {
            "done": done,
            "total": total,
            "percent": round(done / total * 100) if total else 0,
        },
        "config": job.config,
        "output_dir": job.output_dir,
        "error": job.error,
        "elapsed_s": elapsed,
        "created_at": job.created_at,
        "finished_at": job.finished_at,
    }


# ─── Rutas ───────────────────────────────────────────────────────────────────
@api.get("/health")
def health():
    """Comprueba que la API está activa."""
    return {"status": "ok", "version": "0.1"}


# ── Upload de archivos FASTQ ──────────────────────────────────────────────────
@api.post("/upload", status_code=201)
async def upload_fastq(files: List[UploadFile] = File(...)):
    """
    Recibe archivos FASTQ (R1 + R2 por muestra) y los guarda en disco.
    Devuelve un upload_id para usar en POST /jobs.
    """
    if not files:
        raise HTTPException(status_code=400, detail="No se recibieron archivos")

    # Validar extensiones
    allowed = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
    for f in files:
        name = f.filename or ""
        if not any(name.endswith(ext) for ext in allowed):
            raise HTTPException(
                status_code=400,
                detail=f"Archivo no permitido: {name}. Solo FASTQ (.fastq, .fq, .fastq.gz, .fq.gz)",
            )

    upload_id = str(uuid.uuid4())[:8]
    upload_dir = INPUT_DIR / f"upload_{upload_id}"
    upload_dir.mkdir(parents=True, exist_ok=True)

    saved = []
    for f in files:
        dest = upload_dir / (f.filename or f"file_{len(saved)}.fastq")
        with dest.open("wb") as out:
            shutil.copyfileobj(f.file, out)
        saved.append({"name": f.filename, "size_bytes": dest.stat().st_size})

    return {
        "upload_id": upload_id,
        "directory": str(upload_dir),
        "files": saved,
        "count": len(saved),
    }


@api.get("/upload/{upload_id}")
def get_upload(upload_id: str):
    """Lista los archivos de un upload."""
    upload_dir = INPUT_DIR / f"upload_{upload_id}"
    if not upload_dir.exists():
        raise HTTPException(status_code=404, detail="Upload no encontrado")

    files = [
        {"name": f.name, "size_bytes": f.stat().st_size}
        for f in sorted(upload_dir.iterdir())
        if f.is_file()
    ]
    return {"upload_id": upload_id, "directory": str(upload_dir), "files": files}


# ── Jobs ──────────────────────────────────────────────────────────────────────
@api.post("/jobs", status_code=201)
def create_job(req: PipelineRequest):
    """Inicia un nuevo run del pipeline."""
    # Resolver directorio de entrada
    if req.upload_id:
        upload_dir = INPUT_DIR / f"upload_{req.upload_id}"
        if not upload_dir.exists():
            raise HTTPException(status_code=404, detail=f"upload_id '{req.upload_id}' no encontrado")
        resolved_input = str(upload_dir)
    else:
        p = Path(req.input_dir)
        if not p.exists() or not p.is_dir():
            raise HTTPException(status_code=400, detail=f"input_dir no existe: {req.input_dir}")
        resolved_input = req.input_dir

    outdir = req.outdir or str(OUTPUT_DIR / f"run_{int(time.time())}")

    config = req.model_dump(exclude={"outdir", "upload_id", "input_dir"})
    config["input_dir"] = resolved_input
    config["outdir"] = outdir

    job_id = job_store.create_job(config=config, output_dir=outdir)
    job_store.start_job(job_id)

    return {"job_id": job_id, "output_dir": outdir, "status": "queued"}


@api.get("/jobs")
def list_jobs():
    """Lista todos los jobs (más recientes primero)."""
    jobs = sorted(job_store.list_jobs(), key=lambda j: j.created_at, reverse=True)
    return [_job_to_dict(j) for j in jobs]


@api.get("/jobs/{job_id}")
def get_job(job_id: str):
    """Estado detallado de un job."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    return _job_to_dict(job)


@api.get("/jobs/{job_id}/qa")
def get_qa_result(job_id: str):
    """Resultado estructurado del análisis QA."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    if job.qa_result is None:
        raise HTTPException(status_code=404, detail="Este job no tiene resultado QA aún")
    return job.qa_result


@api.get("/jobs/{job_id}/logs")
def get_logs(job_id: str, last: int = 50):
    """Últimas N líneas del log del pipeline."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    return {"job_id": job_id, "logs": job.logs[-last:]}


@api.get("/jobs/{job_id}/results")
def get_results(job_id: str):
    """Archivos de resultado disponibles para un job completado."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    if job.status not in {"done", "failed"}:
        raise HTTPException(status_code=400, detail=f"Job en estado: {job.status}")

    outdir = Path(job.output_dir)
    files = {}
    for fname in job_store.RESULT_FILES:
        fpath = outdir / fname
        if fpath.exists():
            size = fpath.stat().st_size
            files[fname] = {
                "size_bytes": size,
                "size_human": _human_size(size),
            }
    return {"job_id": job_id, "output_dir": str(outdir), "files": files}


@api.get("/jobs/{job_id}/download/{filename}")
def download_file(job_id: str, filename: str):
    """Descarga un archivo de resultado."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")

    # Solo permitir archivos conocidos (seguridad: evitar path traversal)
    if filename not in job_store.RESULT_FILES:
        raise HTTPException(status_code=400, detail="Archivo no permitido")

    fpath = Path(job.output_dir) / filename
    if not fpath.exists():
        raise HTTPException(status_code=404, detail="Archivo no encontrado")

    return FileResponse(path=str(fpath), filename=filename)


@api.delete("/jobs/{job_id}", status_code=204)
def cancel_job(job_id: str):
    """Cancela un job en cola (no puede detener uno que ya corre)."""
    if not job_store.cancel_job(job_id):
        raise HTTPException(
            status_code=400,
            detail="Solo se pueden cancelar jobs en estado 'queued'",
        )


# ── EcoLab ────────────────────────────────────────────────────────────────────

_ECOLAB_ANALYSES = {"alpha", "rarefaction", "nmds", "pcoa", "permanova", "phylum", "genus", "rarefy"}

_R_ALPHA = r"""
suppressPackageStartupMessages(library(vegan))
args    <- commandArgs(trailingOnly=TRUE)
mat     <- read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE)
mat     <- t(mat)
shannon <- round(diversity(mat, "shannon"), 4)
simpson <- round(diversity(mat, "simpson"), 4)
invsim  <- round(diversity(mat, "invsimpson"), 4)
rich    <- specnumber(mat)
total   <- rowSums(mat)
pielou  <- round(shannon / log(rich), 4)
df      <- data.frame(
    Muestra    = rownames(mat),
    Reads      = total,
    Riqueza    = rich,
    Shannon    = shannon,
    Simpson    = simpson,
    InvSimpson = invsim,
    Pielou_J   = pielou,
    stringsAsFactors = FALSE
)
write.table(df, args[2], sep="\t", row.names=FALSE, quote=FALSE)
"""

_R_RAREFACTION = r"""
suppressPackageStartupMessages(library(vegan))
args  <- commandArgs(trailingOnly=TRUE)
mat   <- t(read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE))
n     <- nrow(mat)
cols  <- rainbow(n, alpha=0.85)
min_r <- min(rowSums(mat))
png(args[2], width=1400, height=900, res=130)
par(bg="#0d0d0d", col.axis="#00ffff", col.lab="#ff00ff",
    col.main="#00ffff", fg="#444466", mar=c(5,5,5,2))
rarecurve(mat, step=300, col=cols, lwd=2,
          main="Curvas de Rarefaccion",
          xlab="Reads por muestra", ylab="Riqueza ASV",
          label=FALSE, cex.lab=1.1, cex.main=1.3)
abline(v=min_r, col="#ff00ff", lty=2, lwd=2)
text(min_r, par("usr")[4]*0.93,
     labels=paste0("min = ", min_r),
     col="#ff00ff", cex=0.85, pos=4)
legend("bottomright", legend=rownames(mat), col=cols, lty=1, lwd=2,
       bg="#1a1a2e", text.col="#aaccff", bty="o", cex=0.6, ncol=3)
dev.off()
"""

_R_NMDS = r"""
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly=TRUE)
mat  <- t(read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE))
set.seed(42)
nmds   <- metaMDS(mat, distance="bray", k=2, trymax=100, trace=FALSE)
stress <- round(nmds$stress, 4)
pts    <- as.data.frame(nmds$points)
n      <- nrow(pts)
cols   <- rainbow(n, s=0.9, v=0.95, alpha=0.9)
qual   <- ifelse(stress < 0.1, "excelente", ifelse(stress < 0.2, "aceptable", "poco fiable"))
main_txt <- paste0("NMDS  Bray-Curtis\nStress = ", stress, "  (", qual, ")")
png(args[2], width=1300, height=1000, res=130)
par(bg="#0d0d0d", col.axis="#00ffff", col.lab="#ff00ff",
    col.main="#00ffff", fg="#444466", mar=c(5,5,5,2))
plot(pts[,1], pts[,2], pch=21, bg=cols, col="#ffffff", cex=2.5,
     main=main_txt, xlab="MDS1", ylab="MDS2", cex.main=1.1, cex.lab=1.1)
text(pts[,1], pts[,2], labels=rownames(pts), cex=0.6, col="#00ffff", pos=3, offset=0.4)
abline(h=0, v=0, col="#222244", lty=3)
dev.off()
"""

_R_PHYLUM = r"""
args  <- commandArgs(trailingOnly=TRUE)
mat   <- read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rel   <- sweep(mat, 2, colSums(mat), "/") * 100
means <- sort(rowMeans(rel), decreasing=TRUE)
top_n    <- 10
top_nm   <- names(means)[1:min(top_n, length(means))]
rest_idx <- !(rownames(rel) %in% top_nm)
if(any(rest_idx)) {
    rel_plot <- rbind(rel[top_nm,,drop=FALSE], Otros=colSums(rel[rest_idx,,drop=FALSE]))
} else {
    rel_plot <- rel[top_nm,,drop=FALSE]
}
cols <- c("#00ffff","#ff00ff","#ffff00","#00ff88","#ff6622",
          "#8844ff","#ff0066","#44ccff","#aaff88","#ff9900","#888888")
cols <- cols[seq_len(nrow(rel_plot))]
png(args[2], width=1600, height=860, res=120)
layout(matrix(c(1,2), 1, 2), widths=c(3.2, 1))
# Panel 1: barplot sin leyenda
par(bg="#0d0d0d", col.axis="#00ffff", col.lab="#ff00ff",
    col.main="#00ffff", fg="#444466", mar=c(12,5,4,1))
barplot(as.matrix(rel_plot), col=cols, las=2, border=NA,
        main="Abundancia Relativa por Filo  (Top 10)",
        ylab="Abundancia relativa (%)", cex.names=0.65, cex.main=1.2, axes=TRUE)
# Panel 2: solo la leyenda
par(bg="#0d0d0d", mar=c(1,0,1,1))
plot.new()
legend("topleft", legend=rownames(rel_plot), fill=cols,
       text.col="#aaccff", bg="#1a1a2e", bty="o",
       cex=0.82, title="Filo", title.col="#00ffff", xpd=TRUE)
dev.off()
"""

_R_RAREFY = r"""
suppressPackageStartupMessages(library(vegan))
args  <- commandArgs(trailingOnly=TRUE)
mat   <- read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE)
mat_t <- t(mat)
min_r <- min(rowSums(mat_t))
set.seed(42)
rare  <- t(rrarefy(mat_t, sample=min_r))
rare  <- rare[rowSums(rare) > 0, ]
write.table(data.frame(Sequence_ID=rownames(rare), rare, check.names=FALSE),
            args[2], sep="\t", row.names=FALSE, quote=FALSE)
"""

_R_GENUS = r"""
args <- commandArgs(trailingOnly=TRUE)
df   <- read.table(args[1], header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
meta <- c("Representative_Sequence","ASV","FullTaxonomy","Kingdom","Phylum","Class","Order","Family","Genus","Species")
sample_cols <- setdiff(colnames(df), meta)
# Quitar no clasificados a nivel género
df <- df[!is.na(df$Genus) & trimws(df$Genus) != "" & df$Genus != "Unclassified", ]
if(nrow(df) == 0) stop("No hay ASVs clasificados a nivel genero")
# Agregar por género
agg <- aggregate(df[, sample_cols, drop=FALSE],
                 by=list(Genus=df$Genus), FUN=sum)
rownames(agg) <- agg$Genus
agg$Genus <- NULL
mat <- as.matrix(agg)
# Abundancia relativa respecto al total del dataset
total_reads <- colSums(mat)
rel  <- sweep(mat, 2, total_reads, "/") * 100
# Top 15 por media
means   <- sort(rowMeans(rel), decreasing=TRUE)
top_n   <- 15
top_nm  <- names(means)[1:min(top_n, length(means))]
rest_i  <- !(rownames(rel) %in% top_nm)
if(any(rest_i)) {
    rel_plot <- rbind(rel[top_nm,,drop=FALSE], Otros=colSums(rel[rest_i,,drop=FALSE]))
} else {
    rel_plot <- rel[top_nm,,drop=FALSE]
}
n_col <- nrow(rel_plot)
pal <- colorRampPalette(c("#00ffff","#ff00ff","#ffff00","#00ff88",
                           "#ff6622","#8844ff","#ff0066","#44ccff",
                           "#aaff88","#ff9900","#22dd88","#dd2266",
                           "#aaaaff","#ffaa44","#888888"))(n_col)
png(args[2], width=1600, height=900, res=120)
layout(matrix(c(1,2), 1, 2), widths=c(3.2, 1))
par(bg="#0d0d0d", col.axis="#00ffff", col.lab="#ff00ff",
    col.main="#00ffff", fg="#444466", mar=c(12,5,4,1))
barplot(as.matrix(rel_plot), col=pal, las=2, border=NA,
        main="Abundancia Relativa por Genero  (Top 15)",
        ylab="Abundancia relativa (%)", cex.names=0.65, cex.main=1.2, axes=TRUE)
par(bg="#0d0d0d", mar=c(1,0,1,1))
plot.new()
legend("topleft", legend=rownames(rel_plot), fill=pal,
       text.col="#aaccff", bg="#1a1a2e", bty="o",
       cex=0.75, title="Genero", title.col="#00ffff", xpd=TRUE)
dev.off()
"""

_R_PCOA = r"""
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly=TRUE)
mat  <- t(read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE))
dist_mat <- vegdist(mat, method="bray")
pc  <- cmdscale(dist_mat, k=2, eig=TRUE)
eig <- pc$eig
pct <- round(eig / sum(eig[eig > 0]) * 100, 1)
pts <- as.data.frame(pc$points)
colnames(pts) <- c("PC1","PC2")
n    <- nrow(pts)
cols <- rainbow(n, s=0.9, v=0.95, alpha=0.9)
png(args[2], width=1300, height=1000, res=130)
par(bg="#0d0d0d", col.axis="#00ffff", col.lab="#ff00ff",
    col.main="#00ffff", fg="#444466", mar=c(5,5,5,2))
plot(pts$PC1, pts$PC2, pch=21, bg=cols, col="#ffffff", cex=2.5,
     main="PCoA  Bray-Curtis",
     xlab=paste0("PC1  (", pct[1], "%)"),
     ylab=paste0("PC2  (", pct[2], "%)"),
     cex.main=1.2, cex.lab=1.1)
text(pts$PC1, pts$PC2, labels=rownames(pts), cex=0.6, col="#00ffff", pos=3, offset=0.4)
abline(h=0, v=0, col="#222244", lty=3)
legend("bottomright", legend=rownames(pts), pch=21, pt.bg=cols,
       col="#ffffff", text.col="#aaccff", bg="#1a1a2e", bty="o", cex=0.6, ncol=2)
dev.off()
"""

_R_PERMANOVA = r"""
suppressPackageStartupMessages(library(vegan))
args   <- commandArgs(trailingOnly=TRUE)
mat    <- t(read.table(args[1], header=TRUE, sep="\t", row.names=1, check.names=FALSE))
groups <- factor(strsplit(args[3], ",")[[1]])
if(length(groups) != nrow(mat))
    stop(paste("Grupos:", length(groups), "pero muestras:", nrow(mat)))
set.seed(42)
res <- adonis2(mat ~ groups, method="bray", permutations=999)
# Homogeneidad de dispersión (betadisper)
bd  <- betadisper(vegdist(mat, method="bray"), groups)
bd_test <- permutest(bd, permutations=999)
# Tabla PERMANOVA
perm_df <- data.frame(
    Termino  = rownames(res),
    Df       = res$Df,
    SumOfSqs = round(res$SumOfSqs, 5),
    R2       = round(res$R2, 5),
    F        = round(res$F, 4),
    p_valor  = res$`Pr(>F)`,
    stringsAsFactors = FALSE
)
write.table(perm_df, args[2], sep="\t", row.names=FALSE, quote=FALSE, na="")
# Tabla betadisper aparte (append)
bd_df <- data.frame(
    Termino  = rownames(bd_test$tab),
    Df       = bd_test$tab$Df,
    SumOfSqs = round(bd_test$tab$`Sum Sq`, 5),
    R2       = NA,
    F        = round(bd_test$tab$F, 4),
    p_valor  = bd_test$tab$`Pr(>F)`,
    stringsAsFactors = FALSE
)
cat("\n##BETADISPER\n", file=args[2], append=TRUE)
write.table(bd_df, args[2], sep="\t", row.names=FALSE, quote=FALSE, na="", append=TRUE)
"""

_ECOLAB_SCRIPTS = {
    "alpha":       (_R_ALPHA,       "alpha_diversity.tsv",    "table"),
    "rarefaction": (_R_RAREFACTION, "rarefaction_curves.png", "image"),
    "nmds":        (_R_NMDS,        "beta_nmds.png",          "image"),
    "pcoa":        (_R_PCOA,        "pcoa.png",               "image"),
    "permanova":   (_R_PERMANOVA,   "permanova.tsv",          "permanova"),
    "phylum":      (_R_PHYLUM,      "phylum_abundance.png",   "image"),
    "genus":       (_R_GENUS,       "genus_abundance.png",    "image"),
    "rarefy":      (_R_RAREFY,      "asv_table_rarefied.tsv", "table"),
}


class EcolabRequest(BaseModel):
    groups: Optional[List[str]] = None


def _tsv_to_json(path: Path) -> dict:
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)
    if not rows:
        return {"columns": [], "rows": []}
    return {"columns": rows[0], "rows": rows[1:]}


def _permanova_tsv_to_json(path: Path) -> dict:
    """Parsea el TSV de PERMANOVA que tiene dos secciones separadas por ##BETADISPER."""
    text = path.read_text()
    parts = text.split("##BETADISPER")
    perm = _tsv_to_json_from_text(parts[0].strip())
    beta = _tsv_to_json_from_text(parts[1].strip()) if len(parts) > 1 else {"columns": [], "rows": []}
    return {"permanova": perm, "betadisper": beta}


def _tsv_to_json_from_text(text: str) -> dict:
    reader = csv.reader(text.splitlines(), delimiter="\t")
    rows = list(reader)
    if not rows:
        return {"columns": [], "rows": []}
    return {"columns": rows[0], "rows": rows[1:]}


@api.get("/jobs/{job_id}/ecolab/samples")
def get_ecolab_samples(job_id: str):
    """Devuelve los nombres de las muestras del job (para asignar grupos en PERMANOVA)."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    if job.status != "done":
        raise HTTPException(status_code=400, detail=f"El job no ha terminado (estado: {job.status})")

    asv_table = Path(job.output_dir) / "final_asv_table.tsv"
    if not asv_table.exists():
        raise HTTPException(status_code=404, detail="final_asv_table.tsv no encontrado")

    with open(asv_table) as f:
        header = f.readline().strip().split("\t")
    # Primera columna es Sequence_ID, las demás son muestras
    samples = header[1:] if header else []
    return {"samples": samples}


@api.post("/jobs/{job_id}/ecolab/{analysis}")
def run_ecolab(job_id: str, analysis: str, body: Optional[EcolabRequest] = None):
    """Ejecuta un análisis EcoLab (R/vegan) sobre los resultados del job."""
    if analysis not in _ECOLAB_ANALYSES:
        raise HTTPException(status_code=400, detail=f"Análisis no válido. Opciones: {sorted(_ECOLAB_ANALYSES)}")

    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")
    if job.status != "done":
        raise HTTPException(status_code=400, detail=f"El job no ha terminado (estado: {job.status})")

    outdir = Path(job.output_dir)
    ecolab_dir = outdir / "ecolab_output"
    ecolab_dir.mkdir(exist_ok=True)

    script, out_filename, result_type = _ECOLAB_SCRIPTS[analysis]

    # Seleccionar archivo de entrada según el análisis
    if analysis == "phylum":
        input_file = outdir / "summary_tax_phylum.tsv"
        if not input_file.exists():
            raise HTTPException(status_code=404, detail="summary_tax_phylum.tsv no encontrado en este run")
    elif analysis == "genus":
        input_file = outdir / "asv_resumen.tsv"
        if not input_file.exists():
            raise HTTPException(status_code=404, detail="asv_resumen.tsv no encontrado en este run")
    else:
        input_file = outdir / "final_asv_table.tsv"
        if not input_file.exists():
            raise HTTPException(status_code=404, detail="final_asv_table.tsv no encontrado en este run")

    # Args extra para PERMANOVA
    extra_args: list[str] = []
    if analysis == "permanova":
        if not body or not body.groups:
            raise HTTPException(status_code=400, detail="PERMANOVA requiere 'groups' en el cuerpo de la petición")
        extra_args = [",".join(body.groups)]

    out_path = ecolab_dir / out_filename

    try:
        result = subprocess.run(
            ["Rscript", "--vanilla", "-", str(input_file), str(out_path)] + extra_args,
            input=script,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except FileNotFoundError:
        raise HTTPException(status_code=503, detail="Rscript no está instalado en el servidor")
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="El análisis tardó demasiado (timeout 5 min)")

    if result.returncode != 0 or not out_path.exists():
        detail = result.stderr.strip()[-400:] if result.stderr else "Error desconocido en R"
        raise HTTPException(status_code=500, detail=f"Error en R: {detail}")

    if result_type == "table":
        return {"type": "table", **_tsv_to_json(out_path)}
    elif result_type == "permanova":
        return {"type": "permanova", **_permanova_tsv_to_json(out_path)}
    else:
        return {"type": "image", "filename": out_filename}


@api.get("/jobs/{job_id}/ecolab/image/{filename}")
def get_ecolab_image(job_id: str, filename: str):
    """Sirve una imagen PNG generada por EcoLab."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job no encontrado")

    allowed_images = {
        "rarefaction_curves.png", "beta_nmds.png",
        "phylum_abundance.png", "genus_abundance.png", "pcoa.png",
    }
    if filename not in allowed_images:
        raise HTTPException(status_code=400, detail="Imagen no permitida")

    img_path = Path(job.output_dir) / "ecolab_output" / filename
    if not img_path.exists():
        raise HTTPException(status_code=404, detail="Imagen no encontrada — ejecuta el análisis primero")

    return FileResponse(path=str(img_path), media_type="image/png")


# ─── Helpers ─────────────────────────────────────────────────────────────────
# ─── Registrar router y frontend estático ────────────────────────────────────
app.include_router(api)

_DIST = Path(__file__).parent.parent.parent / "frontend" / "dist"
if _DIST.exists():
    # Archivos estáticos (js, css, imágenes, etc.)
    app.mount("/metagenapp/assets", StaticFiles(directory=str(_DIST / "assets")), name="assets")

    # Fallback SPA: cualquier ruta desconocida devuelve index.html
    @app.get("/metagenapp/{full_path:path}")
    def spa_fallback(full_path: str):
        return FileResponse(str(_DIST / "index.html"))

    @app.get("/metagenapp")
    def spa_root():
        return FileResponse(str(_DIST / "index.html"))


# ─── Helpers ─────────────────────────────────────────────────────────────────
def _human_size(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def serve(host: str = "0.0.0.0", port: int = 8000, reload: bool = False):
    """Punto de entrada para lanzar el servidor."""
    import uvicorn
    uvicorn.run("metagenapp.api.main:app", host=host, port=port, reload=reload)


def main():
    """Entry point para el script metagenapp-api."""
    import argparse
    parser = argparse.ArgumentParser(description="Lanza el servidor API de MetagenApp")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()
    serve(host=args.host, port=args.port, reload=args.reload)
