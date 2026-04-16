import typer
import yaml
from pathlib import Path
from datetime import datetime

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.rule import Rule
from rich import box

console = Console()

CRACK_BANNER = """\
  \u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2557  \u2588\u2588\u2557
  \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255d\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255d\u2588\u2588\u2551 \u2588\u2588\u2554\u255d
  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2551\u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2554\u255d
  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2551\u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2588\u2588\u2557
  \u255a\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2551  \u2588\u2588\u2557\u2588\u2588\u2551  \u2588\u2588\u2551\u255a\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2551  \u2588\u2588\u2557
   \u255a\u2550\u2550\u2550\u2550\u2550\u255d\u255a\u2550\u255d  \u255a\u2550\u255d\u255a\u2550\u255d  \u255a\u2550\u255d \u255a\u2550\u2550\u2550\u2550\u2550\u255d\u255a\u2550\u255d  \u255a\u2550\u255d"""

# ============================================================
# App
# ============================================================
app = typer.Typer(
    help="CRACK — MetagenApp TUI runner",
    invoke_without_command=True
)

# ============================================================
# Profile Loader
# ============================================================
def load_profile(profile_name: str):
    config_path = Path.home() / "profiles.yaml"

    if not config_path.exists():
        console.print("[bold red]✗[/] profiles.yaml not found in home directory.")
        raise typer.Exit(code=1)

    with open(config_path) as f:
        data = yaml.safe_load(f)

    profiles = data.get("profiles", {})

    if profile_name not in profiles:
        console.print(f"[bold red]✗[/] Profile '[bold]{profile_name}[/]' not found.")
        raise typer.Exit(code=1)

    return profiles[profile_name]


# ============================================================
# Banner + header
# ============================================================
def print_banner():
    console.print()
    console.print(Text(CRACK_BANNER, style="bold cyan"))
    console.print(
        Text("  Metabarcoding pipeline · MetagenApp engine", style="dim cyan"),
        justify="left"
    )
    console.print()


def print_config_panel(cfg: dict):
    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 2))
    table.add_column("Key",   style="bold cyan",  no_wrap=True)
    table.add_column("Value", style="white")

    order = [
        ("Input",       cfg["input"]),
        ("Output",      cfg["outdir"]),
        ("Mode",        cfg["mode"]),
        ("Classifier",  cfg["classifier"]),
        ("Marker",      cfg["marker"]),
        ("Model",       cfg["model_type"]),
        ("Threads",     str(cfg["threads"])),
        ("Length",      f"{cfg['min_length']}–{cfg['max_length']} bp"),
        ("Max ambigs",  str(cfg["max_ambigs"])),
        ("Max poly",    str(cfg["max_poly"])),
        ("Centroids",   cfg["extract_centroids"]),
    ]
    if cfg.get("from_step"):
        order.append(("Resume from", cfg["from_step"]))

    for key, val in order:
        table.add_row(key, str(val))

    console.print(Panel(table, title="[bold cyan]Run configuration[/]", border_style="cyan"))


# ============================================================
# Main command
# ============================================================
@app.callback()
def main(
    profile: str = typer.Option(
        None, "--profile",
        help="Load predefined profile configuration"
    ),

    input: Path = typer.Option(
        None, "--input", "-i",
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
        "student", "--mode",
        help="student | premium | turbo | ref | qa",
        case_sensitive=False
    ),

    classifier: str = typer.Option(
        "flat", "--classifier",
        help="flat | naive-v2 | kraken-lite | pro-engine"
    ),

    marker: str = typer.Option(
        "16S", "--marker",
        help="16S | 18S",
        case_sensitive=False
    ),

    min_length: int = typer.Option(None, "--min-length"),
    max_length:  int = typer.Option(None, "--max-length"),
    max_ambigs:  int = typer.Option(None, "--max-ambigs"),
    max_poly:    int = typer.Option(None, "--max-poly"),

    extract_centroids: str = typer.Option(
        "full", "--extract-centroids",
        help="test | student | full"
    ),

    model_type: str = typer.Option(
        "general", "--model-type",
        help="general | oral | gut | skin | env"
    ),

    from_step: str = typer.Option(None, "--from-step"),

    yes: bool = typer.Option(
        False, "--yes", "-y",
        help="Skip confirmation prompt and run immediately"
    ),
):
    """
    Run MetagenApp metabarcoding pipeline with TUI display
    """

    print_banner()

    # --------------------------------------------------------
    # Require --input
    # --------------------------------------------------------
    if input is None:
        console.print("[bold red]✗[/]  Missing option [bold cyan]--input[/] / [bold cyan]-i[/]")
        console.print("   Usage: [bold]crack -i <FASTQ_DIR> [OPTIONS][/]")
        console.print()
        raise typer.Exit(code=1)

    if not input.exists() or not input.is_dir():
        console.print(f"[bold red]✗[/]  Directory not found: [bold]{input}[/]")
        raise typer.Exit(code=1)

    # --------------------------------------------------------
    # Profile
    # --------------------------------------------------------
    if profile:
        profile_data = load_profile(profile)

        mode    = profile_data.get("mode",    mode)
        threads = profile_data.get("threads", threads)

        if min_length is None:
            min_length = profile_data.get("min_length", 250)
        if max_length is None:
            max_length = profile_data.get("max_length", 600)
        if max_ambigs is None:
            max_ambigs = profile_data.get("max_ambigs", 0)
        if max_poly is None:
            max_poly = profile_data.get("max_poly", 8)

        if outdir is None:
            input_resolved = input.resolve()
            sample_name    = input_resolved.name
            timestamp      = datetime.now().strftime("%Y%m%d_%H%M")
            base_output    = Path(profile_data["output_base"])
            outdir         = base_output / f"{sample_name}_{timestamp}"
            outdir.mkdir(parents=True, exist_ok=False)
            console.print(f"[dim]Auto-generated output:[/] [bold]{outdir}[/]")

    # --------------------------------------------------------
    # Defaults
    # --------------------------------------------------------
    if min_length is None: min_length = 250
    if max_length  is None: max_length  = 600
    if max_ambigs  is None: max_ambigs  = 0
    if max_poly    is None: max_poly    = 8

    if outdir is None:
        input_resolved = input.resolve()
        sample_name    = input_resolved.name
        timestamp      = datetime.now().strftime("%Y%m%d_%H%M")
        outdir         = Path.cwd() / f"{sample_name}_{timestamp}"
        outdir.mkdir(parents=True, exist_ok=False)

    # --------------------------------------------------------
    # Normalize
    # --------------------------------------------------------
    mode       = mode.lower()
    classifier = classifier.lower()
    marker     = marker.upper()
    model_type = model_type.lower()

    # --------------------------------------------------------
    # Validate
    # --------------------------------------------------------
    valid_modes = {"student", "premium", "turbo", "ref", "qa"}
    if mode not in valid_modes:
        console.print(f"[bold red]✗[/] Invalid mode: [bold]{mode}[/]")
        console.print(f"   Valid: {', '.join(sorted(valid_modes))}")
        raise typer.Exit(code=1)

    valid_classifiers = {"flat", "naive-v2", "kraken-lite", "pro-engine"}
    if classifier not in valid_classifiers:
        console.print(f"[bold red]✗[/] Invalid classifier: [bold]{classifier}[/]")
        console.print(f"   Valid: flat | naive-v2 | kraken-lite | pro-engine")
        raise typer.Exit(code=1)

    from metagenapp.metagen_config import VALID_MODEL_TYPES
    if model_type not in VALID_MODEL_TYPES:
        console.print(f"[bold red]✗[/] Invalid --model-type: [bold]{model_type}[/]")
        console.print(f"   Valid: {' | '.join(sorted(VALID_MODEL_TYPES))}")
        raise typer.Exit(code=1)

    if classifier == "naive-v2" and model_type != "general":
        from metagenapp.metagen_config import NAIVE_MODELS
        marker_models = NAIVE_MODELS.get(marker, {})
        model_path    = marker_models.get(model_type)
        if model_path and not model_path.exists():
            console.print(f"[bold yellow]⚠[/]  Model '[bold]{model_type}[/]' for {marker} not trained yet.")
            console.print(f"   Expected: {model_path}")
            console.print(f"   Run: metagenapp-train --marker {marker} --model-type {model_type}")
            raise typer.Exit(code=1)

    valid_markers = {"16S", "18S"}
    if marker not in valid_markers:
        console.print(f"[bold red]✗[/] Invalid marker: [bold]{marker}[/]")
        console.print("   Valid: 16S | 18S")
        raise typer.Exit(code=1)

    # --------------------------------------------------------
    # Config panel
    # --------------------------------------------------------
    cfg = dict(
        input=input, outdir=outdir, mode=mode, classifier=classifier,
        marker=marker, model_type=model_type, threads=threads,
        min_length=min_length, max_length=max_length,
        max_ambigs=max_ambigs, max_poly=max_poly,
        extract_centroids=extract_centroids, from_step=from_step,
    )
    print_config_panel(cfg)

    # --------------------------------------------------------
    # Confirm (skip with -y)
    # --------------------------------------------------------
    if not yes:
        console.print(Rule(style="dim cyan"))
        launch = typer.confirm("  Launch pipeline?", default=True)
        if not launch:
            console.print("[dim]Aborted.[/]")
            raise typer.Exit(code=0)
        console.print()

    # --------------------------------------------------------
    # QA MODE
    # --------------------------------------------------------
    if mode == "qa":
        from metagenapp.pipeline.qa import run_qa_contigs
        run_qa_contigs(input_dir=input, outdir=outdir, threads=threads)
        console.print(Rule("[bold cyan]QA mode complete[/]", style="cyan"))
        raise typer.Exit(code=0)

    # --------------------------------------------------------
    # FULL PIPELINE
    # --------------------------------------------------------
    console.print(Rule("[bold cyan]Pipeline running[/]", style="cyan"))

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
    )

    console.print()
    console.print(Rule("[bold cyan]Done[/]", style="cyan"))
    console.print(f"  Results → [bold]{outdir}[/]")
    console.print()


if __name__ == "__main__":
    app()
