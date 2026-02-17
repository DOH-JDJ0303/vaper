#!/usr/bin/env python3
# vaper_tidy.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

from __future__ import annotations

import argparse
import gzip
import os
import re
from pathlib import Path

from vaper_utils import logging_config, load_single_fasta_record

# -------------------------------
# GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

NO_MIX_PATTERN = r'[^ATCGN]'

# -------------------------------
# FUNCTIONS
# -------------------------------

def get_file_stem(filepath: str) -> str:
    """Extract file stem, handling .gz and common FASTA extensions."""
    path = Path(filepath)
    name = path.name
    
    # Remove .gz if present
    if name.endswith('.gz'):
        name = name[:-3]
    
    # Remove FASTA extensions
    if name.endswith('.fa'):
        stem = name[:-3]
    elif name.endswith('.fasta'):
        stem = name[:-6]
    elif name.endswith('.fna'):
        stem = name[:-4]
    else:
        # Just remove the last extension
        stem = Path(name).stem
    
    return stem


def count_characters(s: str) -> dict[str, int]:
    """Count occurrences of each character in a sequence."""
    counts = {}
    for ch in s:
        counts[ch] = counts.get(ch, 0) + 1
    return counts


def remove_terminal_gaps(s: str) -> str:
    """Remove leading and trailing Ns from a sequence."""
    n = len(s)
    li, ti = 0, 0
    
    # ---- Leading gaps ----
    for i, ch in enumerate(s):
        if ch == 'N':
            li = i + 1     # +1 because replace *including* this N
        else:
            break

    # ---- Trailing gaps ----
    for i, ch in enumerate(reversed(s)):
        if ch == 'N':
            ti = i + 1     # +1 so slice removes these Ns too
        else:
            break

    # li = number of leading Ns
    # ti = number of trailing Ns

    LOGGER.info(f"Removing terminal gaps: {li} leading, {ti} trailing")

    # Remove them by slicing:
    if li == 0 and ti == 0:
        return s          # nothing to trim

    if ti == 0:
        return s[li:]      # only leading Ns

    return s[li : n - ti]  # remove both ends


def write_fasta_gzipped(output_path: str, name: str, sequence: str) -> None:
    """Write a FASTA record to a gzipped file."""
    LOGGER.info(f"Writing tidied sequence to: {output_path}")
    with gzip.open(output_path, 'wt') as f:
        f.write(f">{name}\n")
        # Write sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')
    LOGGER.info(f"Output file written successfully ({len(sequence)} bp)")


# -------------------------------
# MAIN
# -------------------------------

def main() -> None:
    version = "1.0.0"

    parser = argparse.ArgumentParser(description="Tidy FASTA sequences by removing mixed sites and/or Ns.")
    parser.add_argument("input", help="Path to assembly in FASTA format.")
    parser.add_argument("--no_mixed_sites", action="store_true", help="Replace non-ATCGN sites with N")
    parser.add_argument("--no_terminal_n", action="store_true", help="Remove terminal Ns")
    parser.add_argument("--no_n", action="store_true", help="Remove all Ns")
    parser.add_argument("--version", action="version", version=f"%(prog)s {version}")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")
    LOGGER.info(f"Input file: {args.input}")

    # Load input sequence
    LOGGER.info("Loading FASTA record...")
    data = load_single_fasta_record(args.input)
    seq = data['sequence']
    name = data['name']
    
    LOGGER.info(f"Original sequence length: {len(seq)} bp")
    LOGGER.info(f"name: {name}")

    # Apply transformations
    if args.no_mixed_sites:
        LOGGER.info("Replacing non-ATCGN sites with N...")
        # Find and count non-ATCGN characters
        mixed_chars = re.findall(NO_MIX_PATTERN, seq)
        if mixed_chars:
            char_counts = count_characters(''.join(mixed_chars))
            LOGGER.info(f"Non-ATCGN characters found: {dict(sorted(char_counts.items()))}")
            total_mixed = sum(char_counts.values())
            LOGGER.info(f"Total non-ATCGN sites replaced: {total_mixed}")
        else:
            LOGGER.info("No non-ATCGN characters found")
        seq = re.sub(NO_MIX_PATTERN, 'N', seq)

    if args.no_n:
        LOGGER.info("Removing all N characters...")
        n_count = seq.count('N')
        seq = re.sub("N", '', seq)
        LOGGER.info(f"Removed {n_count} Ns")
    elif args.no_terminal_n:
        seq = remove_terminal_gaps(seq)

    LOGGER.info(f"Final sequence length: {len(seq)} bp")

    # Generate output filename
    stem = get_file_stem(args.input)
    input_path = Path(args.input)
    output_path = input_path.parent / f"{stem}.tidy.fa.gz"
    
    # Write output
    write_fasta_gzipped(str(output_path), name, seq)
    
    LOGGER.info("Processing complete!")


if __name__ == "__main__":
    main()