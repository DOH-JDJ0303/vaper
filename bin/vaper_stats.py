#!/usr/bin/env python3
# vaper_stats.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from typing import Dict, Tuple, List
import screed
import os

from vaper_utils import logging_config, load_single_fasta_record

# -------------------------------
# GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

SINGLE_IUPAC = ["A","T","C","G"]
MIXED_IUPAC = ["R","Y","S","W","K","M","B","D","H","V"]

# -------------------------------
# FUNCTIONS
# -------------------------------

def read_aligned_records(path: str):
    """Return (ref_rec, query_rec) from an aligned FASTA with exactly two records."""
    recs = [r for r in screed.open(path)]
    if len(recs) != 2:
        raise ValueError(f"Expected exactly two aligned records in {path}, found {len(recs)}")
    return recs[0], recs[1]

def build_multifasta(ref_path: str, ref_contig: str, query_path: str, out_path: str = "multi.fa") -> str:
    """Create a two‑record multi‑FASTA (ref, query)."""
    LOGGER.info("Combining sequences into a multi‑FASTA")
    ref = load_single_fasta_record(ref_path, ref_contig)
    query = load_single_fasta_record(query_path)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f">{ref.name}\n{ref.sequence}\n>{query.name}\n{query.sequence}\n")
    return out_path


def run_mafft(in_fasta: str, out_fasta: str = "aligned.fa", log_path: str = "mafft.log", memsave: bool = False) -> str:
    """
    Run MAFFT to co‑align two sequences.
    Flags:
      --thread -1         use all cores
      --reorder           reorder sequences if needed
      --anysymbol         allow any IUPAC/ambiguous symbols
      --nomemsave         speed over memory
      --adjustdirectionaccurately   auto RC if needed
    """
    LOGGER.info("Aligning sequences with MAFFT")
    cmd = [
        "mafft",
        "--thread", "-1",
        "--reorder",
        "--anysymbol",
        "--adjustdirectionaccurately",
        in_fasta,
    ]

    if memsave:
        cmd = cmd + ["--nomemsave"]

    with open(out_fasta, "w", encoding="utf-8") as out, open(log_path, "w", encoding="utf-8") as log:
        subprocess.run(cmd, stdout=out, stderr=log, check=True)
    return out_fasta

def gather_metrics(aligned_path: str) -> Dict[str, int | float | str]:
    """
    Calculate alignment quality metrics from a pairwise MAFFT alignment.
    
    Logic:
    - Reads reference and query sequences from alignment file
    - Iterates through each aligned position to classify:
      1. Sequence lengths: counts non-gap characters in ref and query
      2. Query base quality: valid IUPAC bases, ambiguous codes, Ns, or invalid chars
      3. Alignment relationships: identical matches, insertions, deletions, or substitutions
    - Computes genome_fraction (non-N query coverage) and identity (match rate) relative to reference length
    
    Returns dictionary with sequence names, raw counts, and computed metrics.
    """
    ref_rec, qry_rec = read_aligned_records(aligned_path)
    ref_seq = ref_rec.sequence.upper()
    qry_seq = qry_rec.sequence.upper()

    if len(ref_seq) != len(qry_seq):
        raise ValueError("Aligned sequences are not the same length; check MAFFT output.")

    m = {
        "ref": ref_rec.name,
        "query": qry_rec.name,
        "rlen": 0,
        "qlen": 0,
        "identical": 0,
        "ins": 0,
        "del": 0,
        "sub": 0,
        "mixed": 0,
        "missing": 0
    }

    for rb, qb in zip(ref_seq, qry_seq):
        # skip double-gap positions
        if qb == '-' and rb == '-':
            LOGGER.warning("Query and reference are both dashes. Skipping!!")
            continue
            
        # track sequence lengths (non-gap characters)
        if rb != "-":
            m["rlen"] += 1
        if qb != "-":
            m["qlen"] += 1

        # count mixed bases
        if qb in MIXED_IUPAC:
            m["mixed"] += 1

        # classify alignment relationship
        if qb not in SINGLE_IUPAC + MIXED_IUPAC + ['-']:
            m["missing"] += 1
        elif qb == rb:
            m["identical"] += 1
        elif qb == "-":
            m["del"] += 1
        elif rb == "-":
            m["ins"] += 1
        else:
            m["sub"] += 1

    if m["rlen"] == 0:
        raise ValueError("Reference length is zero!")

    qlen_nomiss = m["qlen"] - m["missing"]
    m["genome_fraction"] = round(qlen_nomiss / m["rlen"], 3)
    m["identity"] = round(m["identical"] / qlen_nomiss, 3) if qlen_nomiss > 0 else 0

    return m


# -------------------------------
# MAIN
# -------------------------------

def main() -> None:
    version = "1.2.0"

    parser = argparse.ArgumentParser(description="Align two sequences with MAFFT and compute basic metrics.")
    parser.add_argument("--ref", required=True, help="Path to reference (single record).")
    parser.add_argument("--query", required=True, help="Path to query (single record).")
    parser.add_argument("--ref_name", help="Override reference display name.")
    parser.add_argument("--query_name", help="Override query display name.")
    parser.add_argument("--write_csv", action="store_true", help="Also write metrics.csv")
    parser.add_argument("--memsave", action="store_true", help="Use memory save mode with mafft")
    parser.add_argument("--version", action="version", version=f"%(prog)s {version}")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    multi = build_multifasta(args.ref, args.ref_name, args.query)
    aligned = run_mafft(multi, memsave = args.memsave)
    metrics = gather_metrics(aligned)

    if args.ref_name:
        metrics["ref"] = args.ref_name
    if args.query_name:
        metrics["query"] = args.query_name

    out_json = f"{metrics['query']}-{metrics['ref']}.assembly_stats.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
        f.write("\n")
    LOGGER.info(f"Wrote JSON: {out_json}")

    if args.write_csv:
        save_metrics_csv(metrics)


if __name__ == "__main__":
    main()
