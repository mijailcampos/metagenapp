"""
Upload MetagenApp reference models to Zenodo.

Usage:
    python upload_zenodo.py

You will be prompted for your Zenodo API token.
The token is never saved to disk.
"""

import getpass
import hashlib
import json
import subprocess
import sys
import urllib.request
import urllib.error
from pathlib import Path

# ---------------------------------------------------------------------------
# Files to upload
# ---------------------------------------------------------------------------

REFS_ROOT = Path("/data/databases/metagenapp_refs")

FILES = [
    REFS_ROOT / "16S/models/naive_model_silva_v1.pkl",
    REFS_ROOT / "16S/silva_reference.fasta",
]

# ---------------------------------------------------------------------------
# Zenodo metadata
# ---------------------------------------------------------------------------

METADATA = {
    "metadata": {
        "title": "MetagenApp SILVA v1 — 16S classifier model and reference (v1.0)",
        "upload_type": "dataset",
        "description": (
            "Reference models for MetagenApp v1.0, a fast 16S/18S metabarcoding pipeline. "
            "Includes:<br>"
            "<ul>"
            "<li><b>naive_model_silva_v1.pkl</b> — naive-v2 classifier trained on SILVA 138 NR99 "
            "(83,000 taxa, Wang-style bootstrap k-mer confidence, threshold 0.60). "
            "Achieves 91.1% genus resolution on pharyngeal microbiome benchmarks.</li>"
            "<li><b>silva_reference.fasta</b> — SILVA 138 NR99 reference sequences used for "
            "EDLib pairwise alignment fallback (ref mode).</li>"
            "</ul>"
            "See the MetagenApp repository for installation and usage instructions."
        ),
        "creators": [
            {
                "name": "Campos Compeán, José Mijail",
                "affiliation": "",
            }
        ],
        "keywords": [
            "16S rRNA",
            "metagenomics",
            "metabarcoding",
            "taxonomic classification",
            "SILVA",
            "microbiome",
            "bioinformatics",
        ],
        "license": "cc-by-4.0",
        "version": "1.0",
        "related_identifiers": [
            {
                "identifier": "https://github.com/mijailcampos/metagenapp",
                "relation": "isSupplementTo",
                "scheme": "url",
            }
        ],
    }
}

ZENODO_API = "https://zenodo.org/api"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _request(method: str, url: str, token: str, data=None, headers=None) -> dict:
    req_headers = {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}
    if headers:
        req_headers.update(headers)

    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers=req_headers, method=method)

    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        body = e.read().decode()
        raise RuntimeError(f"HTTP {e.code} {e.reason}: {body}") from e


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _upload_file(path: Path, bucket_url: str, token: str):
    """Upload using curl — handles multi-GB files reliably."""
    url = f"{bucket_url}/{path.name}"
    cmd = [
        "curl", "--progress-bar", "--fail-with-body",
        "-X", "PUT",
        "-H", f"Authorization: Bearer {token}",
        "-H", "Content-Type: application/octet-stream",
        "--upload-file", str(path),
        url,
    ]
    result = subprocess.run(cmd, capture_output=False)
    print()
    if result.returncode != 0:
        raise RuntimeError(f"curl upload failed (exit {result.returncode})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resume", metavar="DEP_ID", type=int,
                        help="Resume an existing draft deposition instead of creating a new one")
    args = parser.parse_args()

    print("MetagenApp — Zenodo uploader")
    print("=" * 50)

    # Validate files exist
    for f in FILES:
        if not f.exists():
            print(f"ERROR: File not found: {f}")
            sys.exit(1)
        size_gb = f.stat().st_size / 1e9
        print(f"  {f.name}  ({size_gb:.2f} GB)")

    print()
    token = getpass.getpass("Zenodo API token: ").strip()
    if not token:
        print("ERROR: Token required.")
        sys.exit(1)

    # 1. Create or resume deposition
    if args.resume:
        print(f"\n[1/4] Resuming deposition {args.resume}...")
        dep = _request("GET", f"{ZENODO_API}/deposit/depositions/{args.resume}", token)
    else:
        print("\n[1/4] Creating Zenodo deposition...")
        dep = _request("POST", f"{ZENODO_API}/deposit/depositions", token, data={})

    dep_id = dep["id"]
    bucket_url = dep["links"]["bucket"]
    deposit_url = dep["links"]["html"]
    print(f"  Deposition ID: {dep_id}")
    print(f"  Draft URL: {deposit_url}")

    # 2. Upload files
    print("\n[2/4] Uploading files...")
    for path in FILES:
        print(f"\n  {path.name}  ({path.stat().st_size/1e9:.2f} GB)")
        _upload_file(path, bucket_url, token)
        print(f"  Done: {path.name}")

    # 3. Add metadata
    print("\n[3/4] Setting metadata...")
    _request("PUT", f"{ZENODO_API}/deposit/depositions/{dep_id}", token, data=METADATA)
    print("  Metadata saved.")

    # 4. Done — DO NOT auto-publish, let the user review first
    print("\n[4/4] Upload complete.")
    print()
    print("=" * 50)
    print("NEXT STEPS:")
    print(f"  1. Review your deposit at: {deposit_url}")
    print(f"  2. Check title, description, and files look correct")
    print(f"  3. Click 'Publish' on the Zenodo web page")
    print(f"  4. Copy the DOI and update download_model.py")
    print("=" * 50)
    print()
    print("Deposit ID (save this):", dep_id)


if __name__ == "__main__":
    main()
