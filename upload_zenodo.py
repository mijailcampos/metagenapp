"""
Upload MetagenApp reference files to Zenodo.

Usage:
    python upload_zenodo.py                 # v1.1 reference files (default)
    python upload_zenodo.py --set v1.0      # legacy v1.0 model
    python upload_zenodo.py --resume DEP_ID # continue an existing draft

Files are read from $METAGENAPP_REFS (default: /data/databases/metagenapp_refs).
You will be prompted for your Zenodo API token (scope deposit:write).
The token is never saved to disk. The deposit is left as a draft: review it
on the Zenodo website and publish it there.
"""

import getpass
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.request
import urllib.error
from pathlib import Path

# ---------------------------------------------------------------------------
# Files to upload
# ---------------------------------------------------------------------------

REFS_ROOT = Path(os.getenv("METAGENAPP_REFS", "/data/databases/metagenapp_refs"))

SOFTWARE_DOI_V11 = "10.5281/zenodo.23325861"

SETS = {
    # v1.1: files used by `--mode ref --marker 16S --classifier wang`
    # (metagen_config.SILVA_REFMODE_* and WANG_MODEL_16S)
    "v1.1": {
        "files": [
            "16S_silva1382_full/wang_silva1382_full_species.npz",
            "16S_silva1382_full/silva1382_full_derep.fasta",
            "16S_silva1382_full/silva1382_full_derep.tax",
        ],
        "metadata": {
            "title": "MetagenApp v1.1 — 16S reference files (SILVA 138.2 NR99, Wang classifier)",
            "upload_type": "dataset",
            "description": (
                "Reference files for MetagenApp v1.1 "
                f"(software: https://doi.org/{SOFTWARE_DOI_V11}), the version used in the "
                "MetagenApp benchmark manuscript (ECCB 2026 Research Topic, Frontiers). "
                "They are the files read by <code>--mode ref --marker 16S --classifier wang</code>."
                "<ul>"
                "<li><b>wang_silva1382_full_species.npz</b> — Wang naive Bayes classifier "
                "(k = 8, 100 bootstraps of 1/8 of the k-mers; sparse format) trained on "
                "SILVA 138.2 NR99 Bacteria/Archaea (451,555 sequences) with classes down to "
                "species; MetagenApp reports the lineage to genus at confidence 0.80.</li>"
                "<li><b>silva1382_full_derep.fasta</b> / <b>silva1382_full_derep.tax</b> — "
                "full-length SILVA 138.2 NR99 Bacteria/Archaea, dereplicated "
                "(410,693 unique sequences, 5,128 genera; conflicting taxonomies resolved by "
                "≥90% consensus or LCA), used for reference mapping in ref mode.</li>"
                "</ul>"
                "<p>SHA-256 checksums:</p>{checksums}"
                "<p>Place the files in <code>$METAGENAPP_REFS/16S_silva1382_full/</code>, "
                "or run <code>python download_model.py</code> from the repository. "
                "Build scripts and usage: https://github.com/mijailcampos/metagenapp</p>"
                "<p>Derived from the SILVA SSU Ref NR99 138.2 database (CC BY 4.0). "
                "Please also cite: Quast C. et al. (2013) The SILVA ribosomal RNA gene "
                "database project: improved data processing and web-based tools. "
                "Nucleic Acids Research 41(D1):D590–D596. "
                "https://doi.org/10.1093/nar/gks1219</p>"
            ),
            "creators": [{"name": "Campos Compeán, José Mijail", "affiliation": ""}],
            "keywords": [
                "16S rRNA", "metabarcoding", "taxonomic classification", "SILVA",
                "naive Bayes classifier", "microbiome", "bioinformatics",
            ],
            "license": "cc-by-4.0",
            "version": "1.1",
            "related_identifiers": [
                {"identifier": SOFTWARE_DOI_V11, "relation": "isSupplementTo", "scheme": "doi"},
                {"identifier": "10.1093/nar/gks1219", "relation": "isDerivedFrom", "scheme": "doi"},
            ],
        },
    },
    # v1.0: legacy naive-v2 model (published as 10.5281/zenodo.20076968)
    "v1.0": {
        "files": [
            "16S/models/naive_model_silva_v1.pkl",
            "16S/silva_reference.fasta",
        ],
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
                "<p>SHA-256 checksums:</p>{checksums}"
                "See the MetagenApp repository for installation and usage instructions."
            ),
            "creators": [{"name": "Campos Compeán, José Mijail", "affiliation": ""}],
            "keywords": [
                "16S rRNA", "metagenomics", "metabarcoding", "taxonomic classification",
                "SILVA", "microbiome", "bioinformatics",
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
        },
    },
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
            raw = resp.read()
            return json.loads(raw) if raw else {}  # DELETE returns 204 with no body
    except urllib.error.HTTPError as e:
        body = e.read().decode()
        raise RuntimeError(f"HTTP {e.code} {e.reason}: {body}") from e


def _checksums(path: Path) -> tuple[str, str]:
    """Return (sha256, md5); Zenodo reports md5 for uploaded files."""
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest()


def _upload_file(path: Path, files_url: str, token: str, attempts: int = 6):
    """Upload one file through the InvenioRDM draft files API (init, content, commit).

    The legacy bucket PUT (/api/files/<bucket>/<name>) returned 502 from Zenodo's
    nginx for the 1.2 GB model; the RDM API is the one the Zenodo web uploader uses.
    Zenodo still answers 502 now and then in the middle of a long PUT, and after
    that the entry stays locked ("already being uploaded"), so curl's own retries
    are useless: each attempt deletes the entry and starts it again.
    """
    key = path.name
    for attempt in range(1, attempts + 1):
        if attempt > 1:
            wait = 60 * (attempt - 1)
            print(f"  Attempt {attempt}/{attempts} in {wait} s...")
            time.sleep(wait)
            try:
                _request("DELETE", f"{files_url}/{key}", token)
            except RuntimeError as e:
                print(f"  (delete before retry: {e})")
        _request("POST", files_url, token, data=[{"key": key}])
        cmd = [
            "curl", "--progress-bar", "--fail-with-body",
            # IPv4 + HTTP/1.1: long uploads over IPv6 were reset by peer (curl exit 56)
            "-4", "--http1.1",
            "-X", "PUT",
            # Token read from stdin so it does not show up in `ps`
            "-H", "@-",
            "-H", "Content-Type: application/octet-stream",
            "--upload-file", str(path),
            f"{files_url}/{key}/content",
        ]
        result = subprocess.run(cmd, input=f"Authorization: Bearer {token}\n", text=True)
        print()
        if result.returncode == 0:
            _request("POST", f"{files_url}/{key}/commit", token)
            return
        print(f"  curl failed (exit {result.returncode})")
    raise RuntimeError(f"upload of {key} failed after {attempts} attempts")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", dest="file_set", choices=list(SETS), default="v1.1",
                        help="Which reference set to upload (default: v1.1)")
    parser.add_argument("--resume", metavar="DEP_ID", type=int,
                        help="Resume an existing draft deposition instead of creating a new one")
    args = parser.parse_args()

    file_set = SETS[args.file_set]
    files = [REFS_ROOT / f for f in file_set["files"]]

    print(f"MetagenApp — Zenodo uploader (set {args.file_set})")
    print("=" * 50)

    # Validate files exist
    for f in files:
        if not f.exists():
            print(f"ERROR: File not found: {f}")
            sys.exit(1)
        print(f"  {f.name}  ({f.stat().st_size / 1e9:.2f} GB)")

    # Checksums go into the description and into SHA256SUMS.txt
    print("\nComputing SHA-256 checksums...")
    sums, md5s = {}, {}
    for f in files:
        sums[f.name], md5s[f.name] = _checksums(f)
        print(f"  {sums[f.name]}  {f.name}")

    tmpdir = Path(tempfile.mkdtemp(prefix="metagenapp_zenodo_"))
    sums_file = tmpdir / "SHA256SUMS.txt"
    sums_file.write_text("".join(f"{h}  {name}\n" for name, h in sums.items()))
    upload_paths = files + [sums_file]

    metadata = json.loads(json.dumps(file_set["metadata"]))
    checksums_html = "<ul>" + "".join(
        f"<li><code>{name}</code>: <code>{h}</code></li>" for name, h in sums.items()
    ) + "</ul>"
    metadata["description"] = metadata["description"].replace("{checksums}", checksums_html)

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
    deposit_url = dep["links"]["html"]
    print(f"  Deposition ID: {dep_id}")
    print(f"  Draft URL: {deposit_url}")

    # 2. Upload files
    # Files already complete in the draft (same size and md5) are skipped, so --resume
    # only re-sends what is missing or broken
    files_url = f"{ZENODO_API}/records/{dep_id}/draft/files"
    existing = {
        f["key"]: f
        for f in _request("GET", files_url, token).get("entries", [])
    }
    md5s[sums_file.name] = _checksums(sums_file)[1]

    print("\n[2/4] Uploading files...")
    for path in upload_paths:
        print(f"\n  {path.name}  ({path.stat().st_size/1e9:.2f} GB)")
        done = existing.get(path.name)
        if (done and done.get("status") == "completed"
                and done.get("size") == path.stat().st_size
                and (done.get("checksum") or "").removeprefix("md5:") == md5s[path.name]):
            print(f"  Already in the draft (size and md5 match), skipping")
            continue
        if done:
            # Pending or mismatched entry from an interrupted run: remove it first
            _request("DELETE", f"{files_url}/{path.name}", token)
        _upload_file(path, files_url, token)
        print(f"  Done: {path.name}")

    # 3. Add metadata
    print("\n[3/4] Setting metadata...")
    _request("PUT", f"{ZENODO_API}/deposit/depositions/{dep_id}", token, data={"metadata": metadata})
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
