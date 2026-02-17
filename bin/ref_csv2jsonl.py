#!/usr/bin/env python3

# ref_csv2jsonl.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import List
import screed

from vaper_utils import logging_config

LOGGER = logging_config()


def load_assembly_exact(basename: str) -> str:
    """
    Load an assembly file EXACTLY matching the basename given in the CSV.
    No fallbacks, no extension inference, no .gz checks.
    """

    p = Path(basename)

    if not p.is_file():
        LOGGER.error(f"Assembly file not found: {basename}")
        raise FileNotFoundError(f"Assembly not found: {basename}")

    LOGGER.info(f"Loading assembly: {p}")

    contigs = []
    with screed.open(str(p)) as db:
        for rec in db:
            if rec.sequence:
                contigs.append(rec.sequence)

    if not contigs:
        LOGGER.warning(f"No contigs found in assembly: {p}")
        return ""

    contigs.sort(key=len, reverse=True)
    return "".join(contigs)


def load_refsheet(ref_path: Path):
    ref_data = []

    with ref_path.open("r", encoding="utf-8", newline="") as csv_fh:
        reader = csv.DictReader(csv_fh)

        if "assembly" not in reader.fieldnames:
            raise ValueError(
                f"Refsheet {ref_path} is missing 'assembly' column. "
                f"Columns present: {reader.fieldnames}"
            )

        taxon_count = {}

        for row in reader:
            taxon = row.get("taxon", "unknown")
            segment = row.get("segment") or "wg"
            asm_basename = Path(row["assembly"]).name

            taxon_count.setdefault((taxon, segment), 0)
            taxon_count[(taxon, segment)] += 1
            variant = taxon_count[(taxon, segment)]

            metadata = {k: v for k, v in row.items() if k not in ["assembly", "taxon", "segment"]}

            rec = {
                "taxon": taxon,
                "segment": segment,
                "variant": variant,
                "basename": asm_basename,
                "metadata": metadata,
            }
            ref_data.append(rec)

    return ref_data


def main():
    version = "1.0"

    parser = argparse.ArgumentParser(description="Convert VAPER refsheet CSV → JSONL.GZ")
    parser.add_argument("refsheet", help="Path to refsheet CSV")
    parser.add_argument("--version", action="version", version=version)
    args = parser.parse_args()

    LOGGER.info(f"ref_csv2jsonl v{version}")
    LOGGER.info("Author: Jared Johnson")

    ref_path = Path(args.refsheet)
    if not ref_path.is_file():
        raise FileNotFoundError(f"Refsheet not found: {ref_path}")

    ref_data = load_refsheet(ref_path)

    # Load assemblies using exact match only
    for rec in ref_data:
        rec["sequence"] = load_assembly_exact(rec["basename"])

    # Write JSONL.GZ
    out_path = ref_path.with_suffix(".jsonl.gz")
    LOGGER.info(f"Writing JSONL.GZ → {out_path}")

    with gzip.open(out_path, "wt", encoding="utf-8") as gz:
        for rec in ref_data:
            gz.write(json.dumps(rec, ensure_ascii=False) + "\n")

    LOGGER.info("Done.")


if __name__ == "__main__":
    main()
