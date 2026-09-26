import typer
import yaml
from pathlib import Path
from datetime import datetime

# ============================================================
# Create Typer app
# ============================================================
app = typer.Typer(
    help=(
        "MetagenApp — Fast 16S/18S metabarcoding pipeline.\n\n"
        "Classifies microbial communities from paired-end Illumina reads,\n"
        "producing OTU tables with taxonomic annotation.\n\n"
        "[dim]Tip: run [bold]metagenapp-launch[/bold] for the interactive menu.[/dim]"
    ),
    invoke_without_command=True,
    rich_markup_mode="rich",
)

# ============================================================
# Profile Loader
# ============================================================
def load_profile(profile_name: str):
    config_path = Path.home() / "profiles.yaml"

    if not config_path.exists():
        typer.echo("❌ profiles.yaml not found in home directory.")
        raise typer.Exit(code=1)

    with open(config_path) as f:
        data = yaml.safe_load(f)

    profiles = data.get("profiles", {})

    if profile_name not in profiles:
        typer.echo(f"❌ Profile '{profile_name}' not found.")
        raise typer.Exit(code=1)

    return profiles[profile_name]


# ============================================================
# Main CLI (default command)
# ============================================================
@app.callback()
def main(
    profile: str = typer.Option(
        None,
        "--profile",
        help="Load predefined profile configuration"
    ),

    input: Path = typer.Option(
        ..., "--input", "-i",
        exists=True,
        file_okay=False,
        readable=True,
        help="Input directory containing FASTQ files"
    ),

    outdir: Path = typer.Option(
        None, "--outdir", "-o",
        help="Output directory where results will be stored"
    ),

    threads: int = typer.Option(
        8, "--threads", "-t",
        help="Number of CPU threads to use"
    ),

    mode: str = typer.Option(
        "student",
        "--mode",
        help=(
            "[bold]student[/bold]  quick, low compute (≤10k centroids)\n"
            "[bold]premium[/bold]  balanced, all centroids\n"
            "[bold]turbo[/bold]    maximum speed\n"
            "[bold]ref[/bold]      publication-grade (EDLib, no MAFFT)\n"
            "[bold]qa[/bold]       quality control only"
        ),
        case_sensitive=False,
    ),

    classifier: str = typer.Option(
        "naive-v2",
        "--classifier",
        help=(
            "[bold]naive-v2[/bold]    Wang bootstrap k-mer [dim](recommended)[/dim]\n"
            "[bold]kraken-lite[/bold] k-mer + LCA (SILVA 138)\n"
            "[bold]flat[/bold]        simple assigner\n"
            "[bold]pro-engine[/bold]  high-resolution engine\n"
            "[bold]metaspecies[/bold] SSI+ANI species-level"
        ),
    ),

    marker: str = typer.Option(
        "16S",
        "--marker",
        help="[bold]16S[/bold] bacteria/archaea  [bold]18S[/bold] eukaryotes",
        case_sensitive=False,
    ),

    min_length: int = typer.Option(
        None,
        "--min-length",
        help="Minimum contig length after merging [dim](default: 250)[/dim]",
    ),

    max_length: int = typer.Option(
        None,
        "--max-length",
        help="Maximum contig length after merging [dim](default: 600)[/dim]",
    ),

    max_ambigs: int = typer.Option(
        None,
        "--max-ambigs",
        help="Maximum ambiguous bases allowed [dim](default: 0)[/dim]",
    ),

    max_poly: int = typer.Option(
        None,
        "--max-poly",
        help="Maximum homopolymer run length [dim](default: 8)[/dim]",
    ),

    primer_f: str = typer.Option(
        None,
        "--primer-f",
        help="Forward primer to trim from the 5' end (e.g. 341F) [dim](requires --primer-r)[/dim]",
    ),

    primer_r: str = typer.Option(
        None,
        "--primer-r",
        help="Reverse primer; its reverse complement is trimmed from the 3' end [dim](requires --primer-f)[/dim]",
    ),

    extract_centroids: str = typer.Option(
        "full",
        "--extract-centroids",
        help="Centroid extraction mode: [bold]test[/bold] | [bold]student[/bold] | [bold]full[/bold]",
    ),

    model_type: str = typer.Option(
        "general",
        "--model-type",
        help=(
            "[bold]silva[/bold]    SILVA 138 NR99 — 83K taxa [dim](publication)[/dim]\n"
            "[bold]general[/bold]  1,949 taxa [dim](fast)[/dim]\n"
            "[bold]oral[/bold]     HOMD — oral/pharyngeal\n"
            "[bold]gut[/bold]      gut microbiome\n"
            "[bold]skin[/bold]     skin microbiome\n"
            "[bold]env[/bold]      environmental samples"
        ),
    ),

    from_step: str = typer.Option(
        None,
        "--from-step",
        help="Resume from a specific pipeline step",
    ),
):
    """Run the MetagenApp metabarcoding pipeline."""

    # ========================================================
    # Load profile if provided
    # ========================================================
    if profile:

        profile_data = load_profile(profile)

        mode = profile_data.get("mode", mode)
        threads = profile_data.get("threads", threads)

        if min_length is None:
            min_length = profile_data.get("min_length", 250)

        if max_length is None:
            max_length = profile_data.get("max_length", 600)

        if max_ambigs is None:
            max_ambigs = profile_data.get("max_ambigs", 0)

        if max_poly is None:
            max_poly = profile_data.get("max_poly", 8)

        # auto output folder
        if outdir is None:

            input_resolved = input.resolve()
            sample_name = input_resolved.name
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")

            base_output = Path(profile_data["output_base"])
            outdir = base_output / f"{sample_name}_{timestamp}"

            outdir.mkdir(parents=True, exist_ok=False)

            typer.echo(f"📁 Auto-generated output: {outdir}")

    # ========================================================
    # Global defaults
    # ========================================================
    if min_length is None:
        min_length = 250

    if max_length is None:
        max_length = 600

    if max_ambigs is None:
        max_ambigs = 0

    if max_poly is None:
        max_poly = 8

    # ========================================================
    # Auto output dir (no profile, no --outdir given)
    # ========================================================
    if outdir is None:

        input_resolved = input.resolve()
        sample_name = input_resolved.name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        outdir = Path.cwd() / f"{sample_name}_{timestamp}"
        outdir.mkdir(parents=True, exist_ok=False)
        typer.echo(f"📁 Output: {outdir}")

    # ========================================================
    # Normalize case
    # ========================================================
    mode       = mode.lower()
    classifier = classifier.lower()
    marker     = marker.upper()
    model_type = model_type.lower()

    # ========================================================
    # Validate mode
    # ========================================================
    valid_modes = {"student", "premium", "turbo", "ref", "qa"}

    if mode not in valid_modes:

        typer.echo(f"❌ Invalid mode: {mode}")
        typer.echo(f"Valid modes are: {', '.join(sorted(valid_modes))}")

        raise typer.Exit(code=1)

    # ========================================================
    # Validate classifier
    # ========================================================
    valid_classifiers = {
        "flat",
        "naive-v2",
        "kraken-lite",
        "pro-engine",
        "metaspecies",
    }

    if classifier not in valid_classifiers:

        typer.echo(f"❌ Invalid classifier: {classifier}")
        typer.echo("Valid options: flat | naive-v2 | kraken-lite | pro-engine | metaspecies")

        raise typer.Exit(code=1)

    # ========================================================
    # Validate model_type (only relevant for naive-v2)
    # ========================================================
    from metagenapp.metagen_config import VALID_MODEL_TYPES

    if model_type not in VALID_MODEL_TYPES:
        typer.echo(f"❌ Invalid --model-type: {model_type}")
        typer.echo(f"Valid options: {' | '.join(sorted(VALID_MODEL_TYPES))}")
        raise typer.Exit(code=1)

    if classifier == "naive-v2" and model_type != "general":
        from metagenapp.metagen_config import NAIVE_MODELS
        marker_models = NAIVE_MODELS.get(marker, {})
        model_path = marker_models.get(model_type)
        if model_path and not model_path.exists():
            typer.echo(f"⚠️  Modelo '{model_type}' para {marker} aún no entrenado.")
            typer.echo(f"   Esperado en: {model_path}")
            typer.echo(f"   Usa: metagenapp-train --marker {marker} --model-type {model_type}")
            raise typer.Exit(code=1)

    # ========================================================
    # Validate marker
    # ========================================================
    valid_markers = {"16S", "18S"}

    if marker not in valid_markers:

        typer.echo(f"❌ Invalid marker: {marker}")
        typer.echo("Valid options: 16S | 18S")

        raise typer.Exit(code=1)

    # ========================================================
    # Validate primers (both or none)
    # ========================================================
    if bool(primer_f) != bool(primer_r):
        typer.echo("❌ --primer-f and --primer-r must be given together")
        raise typer.Exit(code=1)

    # ========================================================
    # QA MODE
    # ========================================================
    if mode == "qa":

        from metagenapp.pipeline.qa import run_qa_contigs

        run_qa_contigs(
            input_dir=input,
            outdir=outdir,
            threads=threads
        )

        typer.echo("🛑 QA mode finished. Pipeline stopped before filtering.")

        raise typer.Exit(code=0)

    # ========================================================
    # NORMAL PIPELINE
    # ========================================================
    from metagenapp.pipeline.run_pipeline import run_pipeline

    run_pipeline(
        input_dir=input,
        outdir=outdir,
        threads=threads,
        mode=mode,
        classifier=classifier,
        marker=marker,
        model_type=model_type,
        min_length=min_length,
        max_length=max_length,
        max_ambigs=max_ambigs,
        max_poly=max_poly,
        extract_centroids=extract_centroids,
        from_step=from_step,
        primer_f=primer_f,
        primer_r=primer_r,
    )


if __name__ == "__main__":
    app()