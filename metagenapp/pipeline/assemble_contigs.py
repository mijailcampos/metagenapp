import subprocess
import re
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn
from rich.table import Table
from rich import box

from metagenapp.pipeline.step_tracker import start_step, end_step

try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False

console = Console()


def _parse_vsearch_stats(stderr_text):
    """Extrae pares y reads mergeados del stderr de VSEARCH."""
    pairs = merged = 0
    pct = 0.0
    m = re.search(r'(\d+)\s+Pairs', stderr_text)
    if m:
        pairs = int(m.group(1))
    m = re.search(r'(\d+)\s+Merged\s+\(([0-9.]+)%\)', stderr_text)
    if m:
        merged = int(m.group(1))
        pct = float(m.group(2))
    return pairs, merged, pct


def assemble_contigs(
    files_path,
    input_dir,
    output_fasta,
    output_count,
    threads: int = 1,
):
    """
    Assemble contigs from paired-end FASTQ files using VSEARCH mergepairs.
    Produces merged FASTQ, converts to FASTA, and generates count_table.
    """

    files_path = Path(files_path)
    input_dir = Path(input_dir)
    output_fasta = Path(output_fasta)
    output_count = Path(output_count)

    if not files_path.exists():
        return None, f"File not found: {files_path}", []

    # Leer lista de muestras
    samples_list = []
    with open(files_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3:
                samples_list.append(tuple(parts))

    all_samples_ordered = [s[0] for s in samples_list]
    count_table_data = defaultdict(lambda: defaultdict(int))
    merge_results = []  # (sample, pairs, merged, pct, fastq_path)

    console.print()
    console.rule("[bold cyan]ENSAMBLADO DE CONTIGS[/bold cyan]")
    console.print()

    # ── Fase 1: VSEARCH mergepairs ────────────────────────────────
    with Progress(
        SpinnerColumn(),
        TextColumn("[cyan]{task.description}"),
        BarColumn(bar_width=30),
        MofNCompleteColumn(),
        TextColumn("[dim]{task.fields[info]}"),
        console=console,
        transient=False,
    ) as progress:
        task = progress.add_task(
            "Mergeando pares", total=len(samples_list), info=""
        )

        for sample, r1_file, r2_file in samples_list:
            r1_path = input_dir / r1_file
            r2_path = input_dir / r2_file

            progress.update(task, description=f"[cyan]{sample}", info="")

            if not r1_path.exists() or not r2_path.exists():
                console.print(f"  [yellow]⚠ FASTQ faltante para {sample}[/yellow]")
                progress.advance(task)
                continue

            merged_fastq = output_fasta.parent / f"{sample}_merged.fastq"

            cmd = [
                "vsearch",
                "--fastq_mergepairs", str(r1_path),
                "--reverse", str(r2_path),
                "--fastqout", str(merged_fastq),
                "--fastq_minovlen", "20",
                "--fastq_maxdiffs", "5",
                "--threads", str(threads),
            ]

            try:
                result = subprocess.run(
                    cmd, check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                pairs, merged, pct = _parse_vsearch_stats(result.stderr)
            except subprocess.CalledProcessError as e:
                return None, f"Error VSEARCH en {sample}: {e}", []

            merge_results.append((sample, pairs, merged, pct, merged_fastq))
            color = "green" if pct >= 70 else "yellow" if pct >= 50 else "red"
            progress.update(task, info=f"[{color}]{merged:,} / {pairs:,} ({pct:.1f}%)[/{color}]")
            progress.advance(task)

    # ── Tabla resumen de merge ────────────────────────────────────
    merge_table = Table(box=box.SIMPLE, show_header=True, header_style="bold cyan", pad_edge=False)
    merge_table.add_column("Muestra", style="white")
    merge_table.add_column("Pares", justify="right")
    merge_table.add_column("Mergeados", justify="right")
    merge_table.add_column("%", justify="right")

    total_pairs = total_merged = 0
    for sample, pairs, merged, pct, _ in merge_results:
        color = "green" if pct >= 70 else "yellow" if pct >= 50 else "red"
        merge_table.add_row(
            sample,
            f"{pairs:,}",
            f"{merged:,}",
            f"[{color}]{pct:.1f}%[/{color}]",
        )
        total_pairs += pairs
        total_merged += merged

    overall_pct = (total_merged / total_pairs * 100) if total_pairs else 0
    merge_table.add_section()
    merge_table.add_row(
        "[bold]Total[/bold]",
        f"[bold]{total_pairs:,}[/bold]",
        f"[bold]{total_merged:,}[/bold]",
        f"[bold]{overall_pct:.1f}%[/bold]",
    )
    console.print(merge_table)

    # ── Fase 2: Parsear FASTQs y construir count table ────────────
    merged_fastqs = [(s, fq) for s, _, _, _, fq in merge_results]

    with Progress(
        SpinnerColumn(),
        TextColumn("[cyan]{task.description}"),
        BarColumn(bar_width=30),
        MofNCompleteColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Procesando reads", total=len(merged_fastqs))

        try:
            for sample, merged_fastq in merged_fastqs:
                progress.update(task, description=f"[cyan]Procesando {sample}")
                for record in SeqIO.parse(merged_fastq, "fastq"):
                    count_table_data[str(record.seq)][sample] += 1
                progress.advance(task)
        except Exception as e:
            return None, f"Error parsing merged FASTQ: {e}", []

    # ── Fase 3: Escribir FASTA y count table ─────────────────────
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    with Progress(
        SpinnerColumn(),
        TextColumn("[cyan]Escribiendo contigs..."),
        BarColumn(bar_width=30),
        MofNCompleteColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Escribiendo", total=len(count_table_data))

        with open(output_fasta, "w") as f_out, open(output_count, "w") as f_count:
            f_count.write("ContigID\t" + "\t".join(all_samples_ordered) + "\n")

            for i, (seq, sample_counts) in enumerate(count_table_data.items(), start=1):
                contig_id = f"contig_{i}"
                f_out.write(f">{contig_id}\n{seq}\n")
                row = [str(int(sample_counts.get(s, 0))) for s in all_samples_ordered]
                f_count.write(f"{contig_id}\t" + "\t".join(row) + "\n")
                progress.advance(task)

    console.print(
        f"  [green]✔[/green] [bold]{len(count_table_data):,}[/bold] contigs únicos  "
        f"· [bold]{len(all_samples_ordered)}[/bold] muestras"
    )
    console.print()

    return output_fasta, None, merge_results
