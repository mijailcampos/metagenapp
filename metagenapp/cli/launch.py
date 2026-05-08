"""
metagenapp-launch — Interactive TUI launcher for MetagenApp.

Usage:
    metagenapp-launch
"""

import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml
from rich.console import Console
from rich.prompt import Prompt

console = Console()

BANNER_ART = """\
███╗   ███╗███████╗████████╗ █████╗  ██████╗ ███████╗███╗   ██╗
████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║
██╔████╔██║█████╗     ██║   ███████║██║  ███╗█████╗  ██╔██╗ ██║
██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██╔══╝  ██║╚██╗██║
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝███████╗██║ ╚████║
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝"""

BANNER_SUBTITLE = "  Pipeline Launcher — Metabarcoding 16S/18S · Mijail Campos"


def banner():
    console.print()
    for line in BANNER_ART.splitlines():
        console.print(f"[bold cyan]{line}[/bold cyan]")
    console.print(f"\n[dim]{BANNER_SUBTITLE}[/dim]\n")


def _sep():
    console.print("  [cyan]──────────────────────────────────────────────[/cyan]")


def menu_select(title: str, options: list[tuple[str, str, str]]) -> str:
    """Show a numbered menu and return the selected value."""
    console.print(f"\n  [bold cyan]── {title} ──────────────────────────────────────────[/bold cyan]\n")
    for i, (_, label, desc) in enumerate(options, 1):
        console.print(f"  [magenta][{i}][/magenta] {label}  [dim]{desc}[/dim]")
    console.print()
    while True:
        choice = Prompt.ask("  Elige opción", default="1")
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(options):
                return options[idx][0]
        except ValueError:
            pass
        console.print("  [red]Opción inválida.[/red]")


def _pedir(prompt: str, default: str = "") -> str:
    return Prompt.ask(f"  [cyan]{prompt}[/cyan]", default=default)


# ---------------------------------------------------------------------------
# Parameter selectors
# ---------------------------------------------------------------------------

def elegir_modo() -> str:
    return menu_select("MODO DE EJECUCIÓN", [
        ("student", "student", "Rápido, bajo consumo  (≤10k centroides)"),
        ("premium", "premium", "Balanceado, todos los centroides"),
        ("turbo",   "turbo",   "Máximo rendimiento"),
        ("ref",     "ref",     "Grado publicación  (EDLib, sin MAFFT)"),
        ("qa",      "qa",      "Solo control de calidad"),
    ])


def elegir_clasificador() -> str:
    return menu_select("CLASIFICADOR TAXONÓMICO", [
        ("naive-v2",    "naive-v2",    "Wang bootstrap k-mer  (recomendado)"),
        ("kraken-lite", "kraken-lite", "k-mer + LCA jerárquico  (SILVA 138)"),
        ("flat",        "flat",        "Asignador simple"),
        ("pro-engine",  "pro-engine",  "Motor alta resolución"),
    ])


def elegir_marcador() -> str:
    return menu_select("MARCADOR RIBOSOMAL", [
        ("16S", "16S", "Bacterias/Arqueas  (default)"),
        ("18S", "18S", "Eucariotas"),
    ])


def elegir_modelo(clasificador: str, marcador: str) -> str:
    if clasificador != "naive-v2":
        return "general"
    if marcador == "18S":
        return menu_select("MODELO DE REFERENCIA  (naive-v2)", [
            ("general", "general", "PR2 database — eucariotas  (default)"),
        ])
    return menu_select("MODELO DE REFERENCIA  (naive-v2)", [
        ("silva",   "silva",   "SILVA 138 NR99 — 83K taxa  (publicación)"),
        ("general", "general", "1,949 taxa SILVA  (rápido)"),
        ("oral",    "oral",    "HOMD, 802 taxa  (microbioma oral/faringe)"),
        ("gut",     "gut",     "Subconjunto intestinal"),
        ("skin",    "skin",    "Microbioma de piel"),
        ("env",     "env",     "Muestras ambientales"),
    ])


def elegir_threads() -> int:
    val = _pedir("Número de CPUs a usar", default="8")
    try:
        t = int(val)
        return t if t > 0 else 8
    except ValueError:
        return 8


def elegir_longitudes(min_obs: int | None, med_obs: int | None) -> tuple[int, int]:
    default_min = min_obs if min_obs else 250
    default_max = (med_obs + 3) if med_obs else 600

    if min_obs and med_obs:
        console.print("\n  [green]Sugerencia basada en el QA:[/green]")
        console.print(f"  [cyan]  mínima :[/cyan] [bold]{default_min} bp[/bold]  [dim](mínimo observado)[/dim]")
        console.print(f"  [cyan]  máxima :[/cyan] [bold]{default_max} bp[/bold]  [dim](mediana + 3 bp)[/dim]\n")

    val = _pedir("Longitud mínima", default=str(default_min))
    min_l = int(val) if val.isdigit() and int(val) > 0 else default_min
    val = _pedir("Longitud máxima", default=str(default_max))
    max_l = int(val) if val.isdigit() and int(val) > 0 else default_max
    return min_l, max_l


# ---------------------------------------------------------------------------
# QA step
# ---------------------------------------------------------------------------

def fn_paso_qa(input_dir: str, threads: int) -> tuple[int | None, int | None]:
    qa_dir = Path(tempfile.mkdtemp(prefix="metagen_qa_"))
    console.print(f"\n  [bold cyan]── CONTROL DE CALIDAD (QA) ─────────────────────────[/bold cyan]\n")
    result = subprocess.run([
        "metagenapp", "-i", input_dir, "-o", str(qa_dir),
        "--mode", "qa", "--threads", str(threads),
    ])
    if result.returncode != 0:
        console.print("\n  [red]✗ Error en el QA. Continuando sin sugerencias.[/red]")
        return None, None

    contigs = qa_dir / "assembled_contigs.fasta"
    if not contigs.exists():
        return None, None

    try:
        from Bio import SeqIO
        lens = [len(r.seq) for r in SeqIO.parse(str(contigs), "fasta")]
        if lens:
            return min(lens), int(statistics.median(lens))
    except Exception:
        pass
    return None, None


# ---------------------------------------------------------------------------
# Flows
# ---------------------------------------------------------------------------

def fn_nuevo_analisis():
    console.print("\n  [bold cyan]── NUEVO ANÁLISIS ──────────────────────────────────[/bold cyan]\n")

    input_dir = _pedir("Ruta de la carpeta con los FASTQs")
    if not input_dir:
        return

    output_dir = _pedir("Ruta de la carpeta de salida")
    if not output_dir:
        return

    threads = elegir_threads()

    run_qa = Prompt.ask(
        "\n  [yellow]¿Correr QA primero para sugerencia de longitudes?[/yellow]",
        choices=["s", "n"], default="s",
    )

    min_obs = med_obs = None
    if run_qa == "s":
        min_obs, med_obs = fn_paso_qa(input_dir, threads)

    min_l, max_l = elegir_longitudes(min_obs, med_obs)
    modo = elegir_modo()
    clasificador = elegir_clasificador()
    marcador = elegir_marcador()
    modelo = elegir_modelo(clasificador, marcador)

    console.print()
    _sep()
    console.print("  [bold]Resumen del análisis:[/bold]\n")
    console.print(f"  [cyan]Entrada      :[/cyan] {input_dir}")
    console.print(f"  [cyan]Salida       :[/cyan] {output_dir}")
    console.print(f"  [cyan]Modo         :[/cyan] {modo}")
    console.print(f"  [cyan]Clasificador :[/cyan] {clasificador}")
    console.print(f"  [cyan]Marcador     :[/cyan] {marcador}")
    console.print(f"  [cyan]Modelo       :[/cyan] {modelo}")
    console.print(f"  [cyan]Threads      :[/cyan] {threads}")
    console.print(f"  [cyan]Longitudes   :[/cyan] {min_l}–{max_l} bp")
    _sep()

    confirm = Prompt.ask("\n  [yellow]¿Ejecutar?[/yellow]", choices=["s", "n"], default="n")
    if confirm != "s":
        console.print("  [yellow]Cancelado.[/yellow]\n")
        return

    subprocess.run([
        "metagenapp",
        "-i", input_dir, "-o", output_dir,
        "--mode", modo,
        "--classifier", clasificador,
        "--marker", marcador,
        "--model-type", modelo,
        "--threads", str(threads),
        "--min-length", str(min_l),
        "--max-length", str(max_l),
    ])


def fn_reanudar():
    console.print("\n  [bold cyan]── REANUDAR PIPELINE ───────────────────────────────[/bold cyan]\n")

    input_dir = _pedir("Ruta de la carpeta de entrada (misma del análisis original)")
    if not input_dir:
        return
    output_dir = _pedir("Ruta de la carpeta del análisis original (salida)")
    if not output_dir:
        return
    from_step = _pedir("Nombre del paso desde donde reanudar")

    modo = elegir_modo()
    clasificador = elegir_clasificador()
    marcador = elegir_marcador()
    modelo = elegir_modelo(clasificador, marcador)
    threads = elegir_threads()

    console.print()
    _sep()
    console.print(f"  [bold]Reanudando desde:[/bold] [magenta]{from_step}[/magenta]")
    console.print(f"  [cyan]Entrada      :[/cyan] {input_dir}")
    console.print(f"  [cyan]Salida       :[/cyan] {output_dir}")
    console.print(f"  [cyan]Modo         :[/cyan] {modo}")
    console.print(f"  [cyan]Clasificador :[/cyan] {clasificador}")
    _sep()

    confirm = Prompt.ask("\n  [yellow]¿Ejecutar?[/yellow]", choices=["s", "n"], default="n")
    if confirm != "s":
        console.print("  [yellow]Cancelado.[/yellow]\n")
        return

    subprocess.run([
        "metagenapp",
        "-i", input_dir, "-o", output_dir,
        "--mode", modo,
        "--classifier", clasificador,
        "--marker", marcador,
        "--model-type", modelo,
        "--threads", str(threads),
        "--from-step", from_step,
    ])


def fn_con_perfil():
    console.print("\n  [bold cyan]── EJECUTAR CON PERFIL ─────────────────────────────[/bold cyan]\n")

    profiles_path = Path.home() / "profiles.yaml"
    if profiles_path.exists():
        try:
            with open(profiles_path) as f:
                data = yaml.safe_load(f)
            perfiles = list(data.get("profiles", {}).keys())
            if perfiles:
                console.print("  [cyan]Perfiles disponibles:[/cyan]\n")
                for p in perfiles:
                    console.print(f"  [magenta]  {p}[/magenta]")
                console.print()
        except Exception:
            console.print("  [yellow]No se pudo leer profiles.yaml[/yellow]\n")
    else:
        console.print("  [yellow]profiles.yaml no encontrado en el home.[/yellow]\n")

    input_dir = _pedir("Ruta de entrada (directorio con FASTQs)")
    perfil = _pedir("Nombre del perfil")

    console.print()
    _sep()
    console.print(f"  [cyan]Perfil  :[/cyan] {perfil}")
    console.print(f"  [cyan]Entrada :[/cyan] {input_dir}")
    _sep()

    confirm = Prompt.ask("\n  [yellow]¿Ejecutar?[/yellow]", choices=["s", "n"], default="n")
    if confirm != "s":
        console.print("  [yellow]Cancelado.[/yellow]\n")
        return

    subprocess.run(["metagenapp", "-i", input_dir, "--profile", perfil])


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main():
    try:
        while True:
            banner()
            console.print("  [bold cyan]── OPCIONES ─────────────────────────────────────────[/bold cyan]\n")
            console.print("  [magenta][1][/magenta] Nuevo análisis")
            console.print("  [magenta][2][/magenta] Reanudar desde un paso específico")
            console.print("  [magenta][3][/magenta] Ejecutar con perfil [cyan](profiles.yaml)[/cyan]")
            console.print("  [magenta][0][/magenta] Salir\n")

            op = Prompt.ask("  Elige opción")

            if op == "1":
                fn_nuevo_analisis()
            elif op == "2":
                fn_reanudar()
            elif op == "3":
                fn_con_perfil()
            elif op == "0":
                console.print()
                break
            else:
                console.print("  [red]Opción inválida.[/red]\n")
    except KeyboardInterrupt:
        console.print("\n\n  [dim]Saliendo...[/dim]\n")


if __name__ == "__main__":
    main()
