"""
Download MetagenApp reference models from Zenodo.

Usage:
    python download_model.py                  # downloads to ~/.metagenapp/refs/
    python download_model.py --dest /my/path  # custom destination

After downloading, set METAGENAPP_REFS if you used a custom destination:
    export METAGENAPP_REFS=/my/path
"""

import argparse
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Model registry
# Replace URLs and checksums after uploading to Zenodo.
# ---------------------------------------------------------------------------

ZENODO_DOI = "10.5281/zenodo.20076968"

MODELS = {
    "naive_model_silva_v1.pkl": {
        "url": "https://zenodo.org/record/20076968/files/naive_model_silva_v1.pkl",
        "dest": "16S/models/naive_model_silva_v1.pkl",
        "size_gb": 2.31,
        "sha256": None,
    },
    "silva_reference.fasta": {
        "url": "https://zenodo.org/record/20076968/files/silva_reference.fasta",
        "dest": "16S/silva_reference.fasta",
        "size_gb": 0.022,
        "sha256": None,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST,
                        help=f"Directory to download models into (default: {DEFAULT_DEST})")
    parser.add_argument("--only", choices=list(MODELS.keys()), metavar="FILE",
                        help="Download only this file")
    args = parser.parse_args()

    dest_root = args.dest.resolve()
    print(f"MetagenApp model downloader")
    print(f"Zenodo: https://doi.org/{ZENODO_DOI}")
    print(f"Destination: {dest_root}")
    print()

    targets = {args.only: MODELS[args.only]} if args.only else MODELS

    for filename, info in targets.items():
        dest_path = dest_root / info["dest"]
        print(f"[{filename}]")
        download_file(info["url"], dest_path, filename, info["size_gb"], info["sha256"])
        print()

    print("Done.")
    if str(dest_root) != str(DEFAULT_DEST):
        print(f"\nAdd to your shell config:")
        print(f"  export METAGENAPP_REFS={dest_root}")


if __name__ == "__main__":
    main()
