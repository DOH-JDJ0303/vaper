#!/usr/bin/env python3

# vaper_sink.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import csv
import os
import argparse
import gzip as gzmod
from collections import defaultdict
import re
import subprocess
import zipfile
import uuid
import json
import shutil
import gzip

from vaper_utils import logging_config

# -------------------------------
#  CONFIG
# -------------------------------

LOGGER = logging_config()
ACC_RE = re.compile(r'(?:GCA_|GCF_)\d{9}\.\d')

# -------------------------------
#  FILE I/O
# -------------------------------

def _load_csv(f):
    """Load CSV into list of dicts."""
    data = []
    r = csv.reader(f, delimiter=",")
    header = None
    for i, rec in enumerate(r):
        if i == 0:
            header = rec
            continue
        row = {header[j]: k for j, k in enumerate(rec)}
        data.append(row)
    return data


def load_concatenated_json(path):
    """Load concatenated JSON objects from file."""
    dec = json.JSONDecoder()
    out = []
    with open(path, "r", encoding="utf-8") as f:
        s = f.read()
    i = 0
    n = len(s)
    while i < n:
        while i < n and s[i].isspace():
            i += 1
        if i >= n:
            break
        obj, j = dec.raw_decode(s, i)
        out.append(obj)
        i = j
    return out


def find_fasta_files(acc_dir):
    """Find all FASTA files recursively in directory."""
    exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
    out = []
    if not os.path.isdir(acc_dir):
        return out
    for root, _, files in os.walk(acc_dir):
        for fn in files:
            if fn.lower().endswith(exts):
                out.append(os.path.join(root, fn))
    return sorted(out)


def fasta_iter(path):
    """Yield (header, seq) tuples from FASTA file."""
    opener = gzmod.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        header = None
        chunks = []
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)

# -------------------------------
#  UTILITIES
# -------------------------------

def sanitize_id(s: str) -> str:
    """Clean contig ID: keep first token, replace special chars with underscore."""
    base = s.split()[0]
    base = re.sub(r"[^A-Za-z0-9._\-:+|=]", "_", base)
    return base


def ensure_unique_path(base_dir: str, filename: str) -> str:
    """Make filename unique by appending .2, .3, etc. if needed."""
    stem, ext = os.path.splitext(filename)
    candidate = os.path.join(base_dir, filename)
    i = 2
    while os.path.exists(candidate):
        candidate = os.path.join(base_dir, f"{stem}.{i}{ext}")
        i += 1
    return candidate

# -------------------------------
#  MAIN
# -------------------------------

def main():
    version = 1.0

    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+", help="Path(s) to Sourmash summary files (CSV).")
    parser.add_argument("--outdir", default=".", help="Base output directory (work dir will be created inside)")
    parser.add_argument("--version", action="version", version=str(version), help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    # Setup work directory
    run_id = uuid.uuid4().hex[:12]
    outdir_abs = os.path.abspath(args.outdir)
    workdir = os.path.join(outdir_abs, f"work-{run_id}")
    os.makedirs(workdir, exist_ok=False)
    LOGGER.info(f"Working directory: {workdir}")

    # Extract accessions from input files
    accessions = set()
    query_to_accs = defaultdict(set)

    for filepath in args.input:
        LOGGER.info(f"Processing: {filepath}")
        opener = gzmod.open if filepath.endswith(".gz") else open
        with opener(filepath, "rt", encoding="utf-8") as f:
            rows = _load_csv(f)

        for rec in rows:
            query = rec.get("query_name")
            name = rec.get("name")
            if not name:
                continue
            m = ACC_RE.search(name)
            if m:
                acc = m.group(0)
                accessions.add(acc)
                query_to_accs[query].add(acc)

    if not accessions:
        LOGGER.warning("No accessions found. Nothing to download.")
        return

    acc_list = sorted(accessions)
    LOGGER.info(f"Found {len(acc_list)} unique accessions.")

    old_cwd = os.getcwd()
    try:
        os.chdir(workdir)
        LOGGER.info(f"Changed working directory to {workdir}")

        # Write accession list
        with open("accessions.txt", "w", encoding="utf-8") as fh:
            fh.write("\n".join(acc_list) + "\n")
        LOGGER.info("Wrote accession list: accessions.txt")

        # Download NCBI dataset
        zip_path = "ncbi_dataset.zip"
        cmd = [
            "datasets", "download", "genome", "accession",
            "--assembly-level", "complete",
            "--inputfile", "accessions.txt",
            "--filename", zip_path
        ]
        LOGGER.info(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        if not os.path.isfile(zip_path) or os.path.getsize(zip_path) == 0:
            raise RuntimeError("Download appears to have failed: zip file missing or empty.")
        LOGGER.info(f"Downloaded dataset to {os.path.abspath(zip_path)}")

        # Extract dataset
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(".")
        LOGGER.info("Extracted dataset inside workdir")

        # Parse metadata report
        report_path = os.path.join("ncbi_dataset", "data", "assembly_data_report.jsonl")
        if not os.path.isfile(report_path):
            raise FileNotFoundError(f"Missing report file: {report_path}")

        report = load_concatenated_json(report_path)
        acc_to_species = {}
        for rec in report:
            acc = rec.get("currentAccession") or rec.get("accession")
            species = (rec.get("organism") or {}).get("organismName")
            if acc and species:
                acc_to_species[acc] = species

        # Process contigs
        tmp_contigs_dir = os.path.join(".", "contigs_tmp")
        os.makedirs(tmp_contigs_dir, exist_ok=True)

        contigs_dirname = "contigs"
        metadata_tmp = "metadata.jsonl.gz"
        manifest_tmp = "manifest.csv"
        acc_to_contigs = defaultdict(list)

        with gzip.open(metadata_tmp, "wt", encoding="utf-8") as meta_fh:
            for acc in acc_list:
                acc_dir = os.path.join("ncbi_dataset", "data", acc)
                fasta_files = find_fasta_files(acc_dir)
                if not fasta_files:
                    LOGGER.warning(f"No FASTA found for {acc} in {acc_dir}")
                    continue

                species = acc_to_species.get(acc)

                for fasta in fasta_files:
                    for hdr, seq in fasta_iter(fasta):
                        contig_id = sanitize_id(hdr)
                        tmp_name = f"{contig_id}.fa.gz"
                        tmp_path = ensure_unique_path(tmp_contigs_dir, tmp_name)
                        
                        with gzip.open(tmp_path, "wt", encoding="utf-8") as ofh:
                            ofh.write(f">{contig_id}\n{seq}\n")

                        rec_out = {
                            "name": contig_id,
                            "metadata": {
                                "accession": acc,
                                "species": species
                            }
                        }
                        meta_fh.write(json.dumps(rec_out) + "\n")
                        acc_to_contigs[acc].append(contig_id)

        # Write manifest
        rows = ["sample_id,ref_id"]
        for qname, accs in sorted(query_to_accs.items()):
            for acc in sorted(accs):
                for contig_id in acc_to_contigs.get(acc, []):
                    rows.append(f"{qname},{contig_id}")

        with open(manifest_tmp, "w", encoding="utf-8") as f:
            f.write("\n".join(rows) + "\n")

        # Move outputs to final location
        out_contigs_dir = os.path.join(outdir_abs, contigs_dirname)
        out_metadata = os.path.join(outdir_abs, "metadata.jsonl.gz")
        out_manifest = os.path.join(outdir_abs, "manifest.csv")

        if os.path.exists(out_contigs_dir):
            raise FileExistsError(f"Target contigs directory already exists: {out_contigs_dir}")
        
        shutil.move(tmp_contigs_dir, out_contigs_dir)
        shutil.move(metadata_tmp, out_metadata)
        shutil.move(manifest_tmp, out_manifest)

        LOGGER.info(f"Wrote metadata JSONL: {out_metadata}")
        LOGGER.info(f"Wrote manifest CSV:  {out_manifest}")
        LOGGER.info(f"Wrote contigs dir:   {out_contigs_dir}")

        # Cleanup
        extracted_root = "ncbi_dataset"
        if os.path.isdir(extracted_root):
            shutil.rmtree(extracted_root)
            LOGGER.info(f"Removed extracted directory: {extracted_root}")

        LOGGER.info("Done.")
    finally:
        try:
            os.chdir(old_cwd)
        except Exception:
            pass


if __name__ == "__main__":
    main()