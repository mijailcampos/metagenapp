from pathlib import Path
import re


def generate_input_files(input_dir: Path, output_files: Path):
    """
    Generate a mothur-like .files file from paired-end FASTQ files.
    """

    # Accept .fastq/.fq and their .gz variants (vsearch reads them directly);
    # previously only ".fastq" was detected and compressed input yielded 0 samples
    r1_pattern = re.compile(r"(.*)_R1_.*\.(fastq|fq)(\.gz)?$")
    fastq_exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    samples = {}

    for fastq in input_dir.iterdir():
        if not fastq.name.endswith(fastq_exts):
            continue

        m = r1_pattern.match(fastq.name)
        if m:
            sample = m.group(1)
            samples.setdefault(sample, {})["R1"] = fastq.name
        elif "_R2_" in fastq.name:
            sample = fastq.name.split("_R2_")[0]
            samples.setdefault(sample, {})["R2"] = fastq.name

    with open(output_files, "w") as f:
        for sample, reads in sorted(samples.items()):
            if "R1" in reads and "R2" in reads:
                f.write(f"{sample}\t{reads['R1']}\t{reads['R2']}\n")

    return output_files
