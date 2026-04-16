from pathlib import Path
import typer
from rich.console import Console
from rich.table import Table
from rich import box
from rich.text import Text

from metagenapp.pipeline.assemble_contigs import assemble_contigs
from metagenapp.pipeline.summarize_contigs import summarize_contigs

console = Console()

# (region, primers, min_bp, max_bp)
_MARKERS = [
    ("V1-V2",  "27F / 338R",   250, 340),
    ("V3-V4",  "341F / 806R",  420, 500),
    ("V4",     "515F / 806R",  230, 270),
    ("V4-V5",  "515F / 926R",  350, 430),
    ("V5-V6",  "799F / 1193R", 380, 470),
]

# Margen de tolerancia (bp) para advertencia cuando la longitud queda cerca del rango
_TOLERANCE_BP = 20


def _detect_markers(median_bp: float) -> tuple[list[str], list[str]]:
    """
    Devuelve (matches, near_matches).
    matches      — regiones cuyo rango incluye median_bp exactamente.
    near_matches — regiones dentro del margen de tolerancia (pero fuera del rango).
    """
    matches, near = [], []
    for name, _, lo, hi in _MARKERS:
        if lo <= median_bp <= hi:
            matches.append(name)
        elif (lo - _TOLERANCE_BP) <= median_bp <= (hi + _TOLERANCE_BP):
            near.append(name)
    return matches, near


def _print_marker_section(median_bp: float) -> bool:
    """
    Imprime la tabla de marcadores.
    Retorna True si hay al menos una coincidencia exacta o cercana (tolerancia).
    """
    matches, near_matches = _detect_markers(median_bp)
    all_highlighted = set(matches) | set(near_matches)

    console.print()
    console.rule("[bold cyan]MARCADOR INFERIDO[/bold cyan]")
    console.print()
    console.print(
        f"  Longitud mediana : [bold white]{median_bp:.0f} bp[/bold white]"
    )
    console.print()

    table = Table(
        box=box.SIMPLE_HEAD,
        show_header=True,
        header_style="bold cyan",
        min_width=52,
    )
    table.add_column("Región",  style="bold", no_wrap=True)
    table.add_column("Primers", no_wrap=True)
    table.add_column("Rango típico", justify="right")
    table.add_column("",        no_wrap=True)

    for name, primers, lo, hi in _MARKERS:
        if name in matches:
            region_text  = Text(name,             style="bold green")
            primers_text = Text(primers,           style="bold green")
            range_text   = Text(f"{lo}–{hi} bp",  style="bold green")
            mark_text    = Text("← esta muestra", style="bold green")
        elif name in near_matches:
            region_text  = Text(name,                      style="bold yellow")
            primers_text = Text(primers,                   style="bold yellow")
            range_text   = Text(f"{lo}–{hi} bp",           style="bold yellow")
            mark_text    = Text("← cerca del rango",       style="bold yellow")
        else:
            region_text  = Text(name,             style="dim")
            primers_text = Text(primers,           style="dim")
            range_text   = Text(f"{lo}–{hi} bp",  style="dim")
            mark_text    = Text("")
        table.add_row(region_text, primers_text, range_text, mark_text)

    console.print(table)

    if matches:
        if len(matches) > 1:
            console.print(
                "  [yellow]⚠[/yellow]  Múltiples regiones coinciden — "
                "confirma con tus primers para desambiguar."
            )
    elif near_matches:
        console.print(
            "  [bold yellow]⚠  Longitud mediana ({:.0f} bp) ligeramente fuera del rango "
            "típico de {} — posiblemente correcto dependiendo del organismo y "
            "primers exactos.[/bold yellow]".format(median_bp, ", ".join(near_matches))
        )
        console.print(
            "  [yellow]Continuando pipeline. Verifica tus primers si los resultados "
            "parecen inesperados.[/yellow]"
        )
    else:
        console.print(
            "  [bold red]✗  Longitud mediana ({:.0f} bp) fuera de rango para "
            "cualquier marcador 16S conocido.[/bold red]".format(median_bp)
        )
        console.print(
            "  [red]Verifica que tus datos sean amplicones 16S antes de continuar.[/red]"
        )

    console.print()
    return bool(matches or near_matches)


def _fmt(value, col):
    """Colorea valores según la columna y su importancia."""
    text = Text(str(value))
    if col == "NBases":
        text.stylize("bold green")
    elif col == "Ambigs":
        if float(value) > 0:
            text.stylize("bold yellow")
        else:
            text.stylize("dim")
    elif col == "Polymer":
        if float(value) >= 8:
            text.stylize("bold red")
        else:
            text.stylize("dim")
    else:
        text.stylize("white")
    return text


def _print_summary_table(summary):
    console.print()
    console.rule("[bold cyan]RESULTADO DEL ENSAMBLADO[/bold cyan]")
    console.print()

    # Totales
    console.print(
        f"  Contigs únicos : [bold white]{summary['# of contigs']:,}[/bold white]   "
        f"Secuencias totales : [bold white]{int(summary['total sequences']):,}[/bold white]"
    )
    console.print()

    columns = ["Start", "End", "NBases", "Ambigs", "Polymer", "NumSeqs"]

    table = Table(
        box=box.ROUNDED,
        show_header=True,
        header_style="bold cyan",
        title="[bold]Estadísticas de contigs[/bold]",
        title_style="bold white",
        min_width=60,
    )

    table.add_column("Métrica", style="bold", no_wrap=True)
    for col in columns:
        justify = "right"
        if col == "NBases":
            table.add_column(f"[bold green]{col}[/bold green]", justify=justify)
        else:
            table.add_column(col, justify=justify)

    for stat in ("Minimum", "Median", "Mean"):
        if stat not in summary:
            continue
        vals = summary[stat]
        row = [stat]
        for col in columns:
            row.append(_fmt(vals[col], col))
        table.add_row(*row)

    console.print(table)
    console.print()
    console.print(
        "  [dim]Verde = longitud de contig (NBases) · "
        "Amarillo = bases ambiguas presentes · "
        "Rojo = homopolímero largo[/dim]"
    )
    console.print()


def run_qa_contigs(
    input_dir: Path,
    outdir: Path,
    threads: int = 1,
):
    """
    QA mode: assemble contigs only and stop.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files_path = outdir / "input_samples.files"
    output_fasta = outdir / "assembled_contigs.fasta"
    output_count = outdir / "assembled_contigs.count_table"

    # --------------------------------------------------
    # Generate input_samples.files
    # --------------------------------------------------
    r1_files = sorted([*input_dir.glob("*_R1_*.fastq"), *input_dir.glob("*_R1_*.fastq.gz")])
    r2_files = sorted([*input_dir.glob("*_R2_*.fastq"), *input_dir.glob("*_R2_*.fastq.gz")])

    if not r1_files or not r2_files:
        console.print("  [red]✗[/red] No se encontraron archivos FASTQ R1/R2")
        raise typer.Exit(code=1)

    samples = {}
    for r1 in r1_files:
        sample = r1.name.split("_R1_")[0]
        samples[sample] = {"R1": r1.name}

    for r2 in r2_files:
        sample = r2.name.split("_R2_")[0]
        if sample in samples:
            samples[sample]["R2"] = r2.name

    with open(files_path, "w") as f:
        for sample, files in samples.items():
            if "R1" in files and "R2" in files:
                f.write(f"{sample}\t{files['R1']}\t{files['R2']}\n")
            else:
                console.print(f"  [yellow]⚠[/yellow] Muestra '{sample}' sin par R2 — omitida")

    # --------------------------------------------------
    # Assemble contigs
    # --------------------------------------------------
    output, error, _ = assemble_contigs(
        files_path=files_path,
        input_dir=input_dir,
        output_fasta=output_fasta,
        output_count=output_count,
        threads=threads,
    )

    if error:
        console.print(f"  [red]✗[/red] Error en ensamblado: {error}")
        raise typer.Exit(code=1)

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------
    summary = summarize_contigs(
        fasta_path=output_fasta,
        count_table_path=output_count
    )

    _print_summary_table(summary)

    median_bp = float(summary["Median"]["NBases"])
    marker_ok = _print_marker_section(median_bp)

    if not marker_ok:
        console.print(
            "  [bold red]🛑 Pipeline detenido:[/bold red] los datos no corresponden "
            "a ningún marcador 16S soportado por MetagenApp."
        )
        raise typer.Exit(code=1)

    return output_fasta
