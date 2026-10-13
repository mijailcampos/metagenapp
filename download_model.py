"""
Download MetagenApp reference files from Zenodo.

Usage:
    python download_model.py                  # v1.1 files to ~/.metagenapp/refs/
    python download_model.py --set v1.0       # legacy v1.0 naive-v2 model
    python download_model.py --dest /my/path  # custom destination

After downloading, set METAGENAPP_REFS if you used a custom destination:
    export METAGENAPP_REFS=/my/path
"""

import argparse
import gzip
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Model registry
# One entry per Zenodo record. "record" is the Zenodo record id (None until
# the record is published); "dest" is relative to METAGENAPP_REFS and matches
# the paths in metagenapp/metagen_config.py. Files with "gunzip_sha256" are
# downloaded gzipped and decompressed to "dest".
# ---------------------------------------------------------------------------

SETS = {
    # v1.1: --mode ref --marker 16S --classifier wang (version used in the manuscript)
    "v1.1": {
        "record": 23326733,
        "doi": "10.5281/zenodo.23326733",
        "files": {
            "wang_silva1382_full_species.npz": {
                "dest": "16S_silva1382_full/wang_silva1382_full_species.npz",
                "size_gb": 1.22,
                "sha256": "9b25c6708b8f4cc2c2e2a979ab6ce7b98c023f5acb27dbab500b3e2669950beb",
            },
            "silva1382_full_derep.fasta.gz": {
                "dest": "16S_silva1382_full/silva1382_full_derep.fasta",
                "size_gb": 0.10,
                "sha256": "e7a9871f5c09778ea82892405cf16acf75194bb97d6d4cf8cb7e5d384a36501d",
                # Stored gzipped on Zenodo; "dest" and this checksum are the decompressed file.
                "gunzip_sha256": "b153052f725e691b94c719ac0ee14b9301127c3b5a23c8e1dc8c8dc6d29f411c",
            },
            "silva1382_full_derep.tax": {
                "dest": "16S_silva1382_full/silva1382_full_derep.tax",
                "size_gb": 0.04,
                "sha256": "549f49e6d12aa16281b262ea2ddf3456fb0931adffd80c4046b587cbec8f401c",
            },
        },
    },
    # v1.0: legacy naive-v2 model trained on SILVA 138
    "v1.0": {
        "record": 20076968,
        "doi": "10.5281/zenodo.20076968",
        "files": {
            "naive_model_silva_v1.pkl": {
                "dest": "16S/models/naive_model_silva_v1.pkl",
                "size_gb": 2.31,
                "sha256": None,
            },
            "silva_reference.fasta": {
                "dest": "16S/silva_reference.fasta",
                "size_gb": 0.022,
                "sha256": None,
            },
        },
    },
}

DEFAULT_DEST = Path.home() / ".metagenapp" / "refs"


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while data := f.read(chunk):
            h.update(data)
    return h.hexdigest()


def _progress_hook(filename: str, size_gb: float):
    downloaded = [0]

    def hook(count, block_size, total_size):
        downloaded[0] += block_size
        if total_size > 0:
            pct = min(100, downloaded[0] * 100 / total_size)
            mb = downloaded[0] / 1e6
            bar = "#" * int(pct / 2)
            print(f"\r  [{bar:<50}] {pct:5.1f}%  {mb:,.0f} MB", end="", flush=True)
        else:
            mb = downloaded[0] / 1e6
            print(f"\r  {mb:,.0f} MB downloaded...", end="", flush=True)

    return hook


def download_file(url: str, dest: Path, filename: str, size_gb: float, sha256: str | None):
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists():
        print(f"  Already exists: {dest}")
        if sha256:
            print("  Verifying checksum...", end=" ")
            if _sha256(dest) == sha256:
                print("OK")
                return
            else:
                print("MISMATCH — re-downloading")
        else:
            return

    print(f"  Downloading {filename} ({size_gb:.1f} GB)...")
    try:
        urllib.request.urlretrieve(url, dest, reporthook=_progress_hook(filename, size_gb))
        print()  # newline after progress bar
    except Exception as e:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed: {e}") from e

    if sha256:
        print("  Verifying checksum...", end=" ")
        actual = _sha256(dest)
        if actual != sha256:
            dest.unlink()
            raise RuntimeError(f"Checksum mismatch for {filename}.\n  expected: {sha256}\n  got:      {actual}")
        print("OK")


def download_gzipped(url: str, dest: Path, filename: str, size_gb: float,
                     sha256_gz: str, sha256: str):
    if dest.exists():
        print(f"  Already exists: {dest}")
        print("  Verifying checksum...", end=" ")
        if _sha256(dest) == sha256:
            print("OK")
            return
        print("MISMATCH — re-downloading")

    gz_path = dest.with_name(dest.name + ".gz")
    download_file(url, gz_path, filename, size_gb, sha256_gz)

    print(f"  Decompressing to {dest.name}...", end=" ")
    tmp = dest.with_name(dest.name + ".part")
    with gzip.open(gz_path, "rb") as src, open(tmp, "wb") as out:
        shutil.copyfileobj(src, out, 1 << 20)
    actual = _sha256(tmp)
    if actual != sha256:
        tmp.unlink()
        raise RuntimeError(f"Checksum mismatch for {dest.name}.\n  expected: {sha256}\n  got:      {actual}")
    tmp.replace(dest)
    gz_path.unlink()
    print("OK")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", dest="file_set", choices=list(SETS), default="v1.1",
                        help="Which reference set to download (default: v1.1)")
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST,
                        help=f"Directory to download models into (default: {DEFAULT_DEST})")
    parser.add_argument("--only", metavar="FILE",
                        help="Download only this file of the selected set")
    args = parser.parse_args()

    file_set = SETS[args.file_set]
    if file_set["record"] is None:
        print(f"ERROR: the {args.file_set} reference files are not on Zenodo yet.")
        sys.exit(1)
    if args.only and args.only not in file_set["files"]:
        parser.error(f"--only must be one of: {', '.join(file_set['files'])}")

    dest_root = args.dest.resolve()
    print(f"MetagenApp model downloader (set {args.file_set})")
    print(f"Zenodo: https://doi.org/{file_set['doi']}")
    print(f"Destination: {dest_root}")
    print()

    targets = {args.only: file_set["files"][args.only]} if args.only else file_set["files"]

    for filename, info in targets.items():
        url = f"https://zenodo.org/records/{file_set['record']}/files/{filename}?download=1"
        dest_path = dest_root / info["dest"]
        print(f"[{filename}]")
        if "gunzip_sha256" in info:
            download_gzipped(url, dest_path, filename, info["size_gb"], info["sha256"], info["gunzip_sha256"])
        else:
            download_file(url, dest_path, filename, info["size_gb"], info["sha256"])
        print()

    print("Done.")
    if str(dest_root) != str(DEFAULT_DEST):
        print(f"\nAdd to your shell config:")
        print(f"  export METAGENAPP_REFS={dest_root}")


if __name__ == "__main__":
    main()
