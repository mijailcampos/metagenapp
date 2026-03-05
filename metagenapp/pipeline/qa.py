from pathlib import Path
import typer


def run_qa_contigs(
    input_dir: Path,
    outdir: Path,
    threads: int = 1,
):
    """
    QA mode: assemble contigs only and stop.
    """

    typer.echo("🔍 QA mode: running contig assembly only")

    from metagenapp.pipeline.assemble_contigs import assemble_contigs
    from metagenapp.pipeline.summarize_contigs import summarize_contigs

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files_path = outdir / "input_samples.files"
    output_fasta = outdir / "assembled_contigs.fasta"
    output_count = outdir / "assembled_contigs.count_table"

    # --------------------------------------------------
    # Generate input_samples.files
    # --------------------------------------------------
    typer.echo("📝 Generating input_samples.files")

    r1_files = sorted(input_dir.glob("*_R1_*.fastq"))
    r2_files = sorted(input_dir.glob("*_R2_*.fastq"))

    if not r1_files or not r2_files:
        typer.echo("❌ No FASTQ R1/R2 files found")
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

    # --------------------------------------------------
    # Assemble contigs
    # --------------------------------------------------
    typer.echo("🧬 Assembling contigs")

    output, error = assemble_contigs(
        files_path=files_path,
        input_dir=input_dir,
        output_fasta=output_fasta,
        output_count=output_count,
    )

    if error:
        typer.echo(f"❌ Error during contig assembly: {error}")
        raise typer.Exit(code=1)

    typer.echo(f"✅ Contigs assembled: {output_fasta}")
    typer.echo(f"📊 Count table: {output_count}")

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------
    typer.echo("\n📊 CONTIG SUMMARY (QA MODE)")
    typer.echo("-" * 40)

    summary = summarize_contigs(
        fasta_path=output_fasta,
        count_table_path=output_count
    )

    typer.echo(f"Number of contigs: {summary['# of contigs']}")
    typer.echo(f"Total sequences: {summary['total sequences']}")

    typer.echo("\nMinimum values:")
    for k, v in summary["Minimum"].items():
        typer.echo(f"  {k}: {v}")

    typer.echo("\nMedian values:")
    for k, v in summary["Median"].items():
        typer.echo(f"  {k}: {v}")

    typer.echo("\nMean values:")
    for k, v in summary["Mean"].items():
        typer.echo(f"  {k}: {v}")

    return output_fasta
