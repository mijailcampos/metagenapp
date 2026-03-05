import typer
import yaml
from pathlib import Path
from datetime import datetime

# ============================================================
# Create Typer app
# ============================================================
app = typer.Typer(
    help="MetagenApp CLI — reproducible metabarcoding pipeline",
    invoke_without_command=True
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
        help="Input directory with FASTQ files"
    ),
    outdir: Path = typer.Option(
        None, "--outdir", "-o",
        help="Output directory for results"
    ),
    threads: int = typer.Option(
        8, "--threads", "-t",
        help="Number of CPU threads to use"
    ),
    mode: str = typer.Option(
        "student",
        "--mode",
        help="Pipeline mode: student | premium | turbo | ref | qa",
        case_sensitive=False
    ),
    classifier: str = typer.Option(
        "flat",
        "--classifier",
        help="Taxonomic classifier engine: flat | pro-engine",
        case_sensitive=False
    ),
    marker: str = typer.Option(
        "16S",
        "--marker",
        help="Marker type: 16S | 18S",
        case_sensitive=False
    ),
    min_length: int = typer.Option(None, "--min-length"),
    max_length: int = typer.Option(None, "--max-length"),
    max_ambigs: int = typer.Option(None, "--max-ambigs"),
    max_poly: int = typer.Option(None, "--max-poly"),
    extract_centroids: str = typer.Option(
        "full",
        "--extract-centroids",
        help="Centroid extraction mode: test | student | full"
    ),
    from_step: str = typer.Option(
        None,
        "--from-step",
        help="Resume pipeline from a specific step"
    ),
):
    """
    Run MetagenApp metabarcoding pipeline
    """

    # ========================================================
    # Load profile if provided
    # ========================================================
    if profile:
        profile_data = load_profile(profile)

        # Mode & threads from profile
        mode = profile_data.get("mode", mode)
        threads = profile_data.get("threads", threads)

        # Filters (CLI overrides profile)
        if min_length is None:
            min_length = profile_data.get("min_length", 250)

        if max_length is None:
            max_length = profile_data.get("max_length", 600)

        if max_ambigs is None:
            max_ambigs = profile_data.get("max_ambigs", 0)

        if max_poly is None:
            max_poly = profile_data.get("max_poly", 8)

        # Auto-generate outdir if not provided
        if outdir is None:
            input_resolved = input.resolve()
            sample_name = input_resolved.name
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")

            base_output = Path(profile_data["output_base"])
            outdir = base_output / f"{sample_name}_{timestamp}"
            outdir.mkdir(parents=True, exist_ok=False)

            typer.echo(f"📁 Auto-generated output: {outdir}")

    # ========================================================
    # Global defaults if still None
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
    # Require outdir if no profile
    # ========================================================
    if outdir is None:
        typer.echo("❌ --outdir is required if no profile is used.")
        raise typer.Exit(code=1)

    # ========================================================
    # Normalize case
    # ========================================================
    mode = mode.lower()
    classifier = classifier.lower()
    marker = marker.upper()

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
    valid_classifiers = {"flat", "pro-engine"}
    if classifier not in valid_classifiers:
        typer.echo(f"❌ Invalid classifier: {classifier}")
        typer.echo("Valid options: flat | pro-engine")
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
        min_length=min_length,
        max_length=max_length,
        max_ambigs=max_ambigs,
        max_poly=max_poly,
        extract_centroids=extract_centroids,
        from_step=from_step
    )


if __name__ == "__main__":
    app()