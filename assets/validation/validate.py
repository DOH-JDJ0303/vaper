#!/usr/bin/env python3

# validate.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import logging
import subprocess
import sys
import textwrap
from collections import defaultdict
from typing import List, Dict, Any, Tuple, Optional
import os
import math
import itertools
import tempfile

import screed
import numpy as np

import boto3

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# IUPAC ambiguity codes for nucleotide matching
IUPAC_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'U': {'T'},  # Treat U as T
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'}
}

# --------------------------
# Core Functions
# --------------------------

def load_csv_data(file_path: str, delimiter: str = ',') -> List[Dict[str, str]]:
    """Load CSV data and return list of dictionaries."""
    with open(file_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def is_s3_uri(path: str) -> bool:
    return isinstance(path, str) and path.startswith("s3://")


def download_s3_uri(uri: str, tmpdir: str, s3_client) -> str:
    """
    Download s3://bucket/key to tmpdir and return local path.
    Filename is derived from the key basename.
    """
    rest = uri[5:]
    if "/" not in rest:
        raise ValueError(f"Invalid S3 URI (missing key): {uri}")
    bucket, key = rest.split("/", 1)
    local_path = os.path.join(tmpdir, os.path.basename(key) or "downloaded.fa.gz")
    s3_client.download_file(bucket, key, local_path)
    return local_path


def _safe_div(n: float, d: float) -> float:
    if d == 0 or (isinstance(d, float) and math.isnan(d)) or (isinstance(n, float) and math.isnan(n)):
        return float("nan")
    return n / d


def calculate_accuracy(rec: Dict[str, Any], ignore_termini: bool = False) -> None:
    """
    Accuracy + completeness, where one of the sequences is 'truth'.

    accuracy     = match / valid
    completeness = valid_test / valid_truth
    """
    if str(rec['seq1']).lower() == 'truth':
        tsfx, ssfx = "s1", "s2"
        termini_to_ignore = rec['del_terminal']
    else:
        tsfx, ssfx = "s2", "s1"
        termini_to_ignore = rec['ins_terminal']

    valid_truth = float(rec["valid_" + tsfx])
    valid_test  = float(rec["valid_" + ssfx])
    valid_both  = float(rec["valid"])
    match       = float(rec["match"])

    if ignore_termini:
        valid_truth = max(0.0, valid_truth - termini_to_ignore)
        valid_both  = max(0.0, valid_both  - termini_to_ignore)

    rec["accuracy"]     = _safe_div(match, valid_both)
    rec["completeness"] = _safe_div(valid_test, valid_truth)


def calculate_precision(rec: Dict[str, Any], ignore_termini: bool = False) -> None:
    """
    Precision-like score when there is no truth. Kept consistent with your current intent:
      precision = match / valid
    """
    match       = float(rec["match"])
    valid_s1    = float(rec["valid_s1"])
    valid_s2    = float(rec["valid_s2"])
    valid_both  = float(rec["valid"])
    valid_total = valid_both + (valid_s1 - valid_both) + (valid_s2 - valid_both)

    termini_to_ignore = float(rec['ins_terminal'] + rec['del_terminal'])

    if ignore_termini:
        valid_both  = valid_both  - termini_to_ignore
        valid_total = valid_total - termini_to_ignore
    

    rec["precision"]    = _safe_div(match, valid_both)
    rec["completeness"] = _safe_div(valid_both, valid_total)


def calculate_sequence_metrics(seq1, seq2, ignore_termini: bool = False) -> Dict[str, Any]:
    """
    Calculate similarity metrics between two aligned sequences using NumPy arrays.
    """
    logging.info(f'Calculating metrics (ref: {seq1.name}, query: {seq2.name})')

    s1 = np.array(list(seq1.sequence.upper()))
    s2 = np.array(list(seq2.sequence.upper()))

    if s1.shape != s2.shape:
        sys.exit("Error: Sequences are not the same length!")

    metrics: Dict[str, Any] = {
        'match': 0,
        'match_iupac': 0,
        'ins': 0,
        'del': 0,
        'sub': 0,
        'invalid': 0,
        'valid': 0,
        'valid_s1': 0,
        'invalid_s1': 0,
        'valid_s2': 0,
        'invalid_s2': 0,
        'del_terminal': 0,
        'del_internal': 0,
        'ins_terminal': 0,
        'ins_internal': 0,
    }
    metrics['seq1'] = seq1.name
    metrics['seq2'] = seq2.name

    # Validate sequences against IUPAC codes
    valid_bases = np.array(list(IUPAC_CODES.keys()) + ['-'])
    valid_s1 = np.isin(s1, valid_bases)
    valid_s2 = np.isin(s2, valid_bases)

    # Create mask for valid positions (excluding double gaps)
    mask = valid_s1 & valid_s2 & ~((s1 == '-') & (s2 == '-'))

    # Record validation statistics
    metrics['valid_s1'] = int(valid_s1.sum())
    metrics['invalid_s1'] = int((~valid_s1).sum())
    metrics['valid_s2'] = int(valid_s2.sum())
    metrics['invalid_s2'] = int((~valid_s2).sum())
    metrics['valid'] = int(mask.sum())
    metrics['invalid'] = int((~mask).sum())
    metrics['invalid_char_s1'] = list(set(list(s1[~valid_s1])))
    metrics['invalid_char_s2'] = list(set(list(s2[~valid_s2])))

    # Extract valid positions
    b1 = s1[mask]
    b2 = s2[mask]
    length = len(b1)
    position_indices = np.arange(length)

    # Find sequence boundaries (first and last non-gap positions)
    non_gap_b1 = np.flatnonzero(b1 != '-')
    first_b1 = int(non_gap_b1[0]) if non_gap_b1.size > 0 else length
    last_b1  = int(non_gap_b1[-1]) if non_gap_b1.size > 0 else -1

    non_gap_b2 = np.flatnonzero(b2 != '-')
    first_b2 = int(non_gap_b2[0]) if non_gap_b2.size > 0 else length
    last_b2  = int(non_gap_b2[-1]) if non_gap_b2.size > 0 else -1

    # Count matches
    metrics['match'] = int((b1 == b2).sum())

    # Calculate indels (relative to seq1)
    deletions_mask  = (b2 == '-') & (b1 != '-')
    insertions_mask = (b1 == '-') & (b2 != '-')

    metrics['del'] = int(deletions_mask.sum())
    metrics['ins'] = int(insertions_mask.sum())

    # Classify deletions as terminal vs internal (using seq2 boundaries)
    del_terminal_mask = deletions_mask & ((position_indices < first_b2) | (position_indices > last_b2))
    metrics['del_terminal']  = int(del_terminal_mask.sum())
    metrics['del_internal']  = int(deletions_mask.sum() - metrics['del_terminal'])

    # Classify insertions as terminal vs internal (using seq1 boundaries)
    ins_terminal_mask = insertions_mask & ((position_indices < first_b1) | (position_indices > last_b1))
    metrics['ins_terminal']  = int(ins_terminal_mask.sum())
    metrics['ins_internal']  = int(insertions_mask.sum() - metrics['ins_terminal'])

    # Calculate substitutions with IUPAC awareness
    difference_mask = (b1 != b2)

    single_bases = np.array(list('ATCG-'))
    single_base_mask = (np.isin(b1, single_bases) &
                        np.isin(b2, single_bases) &
                        difference_mask)

    simple_substitutions = (single_base_mask & (b1 != '-') & (b2 != '-'))
    metrics['sub'] = int(simple_substitutions.sum())

    iupac_differences = (difference_mask &
                         ~single_base_mask &
                         (b1 != '-') &
                         (b2 != '-'))

    if np.any(iupac_differences):
        for base1, base2 in zip(b1[iupac_differences], b2[iupac_differences]):
            if IUPAC_CODES[base1] & IUPAC_CODES[base2]:
                metrics['match_iupac'] += 1
            else:
                metrics['sub'] += 1

    # Determine whether truth is involved
    truth_involved = (str(seq1.name).lower() == 'truth') or (str(seq2.name).lower() == 'truth')
    if truth_involved:
        calculate_accuracy(metrics, ignore_termini)
    else:
        calculate_precision(metrics, ignore_termini)

    return metrics


def read_alignment_file(input_file: str) -> List:
    """Read alignment file and return list of sequence records."""
    logging.info('Reading alignment')
    return list(screed.open(input_file))


def run_mafft_alignment(sample_name: str, segment: str, input_file: str, outdir: str) -> str:
    logging.info('Aligning sequences with MAFFT')
    output_file = os.path.join(outdir, f'{sample_name}_{segment}.aligned.fa')

    with open(output_file, 'w') as out:
        subprocess.run(
            ["mafft", "--auto", input_file],
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True
        )
    return output_file


def create_multi_fasta(sample_name: str, segment: str, segment_data: List[Dict[str, str]], outdir: str) -> Tuple[str, Dict[str, int]]:
    logging.info(f'Combining sequences for {sample_name} {segment} into a multi-FASTA')
    output_file = os.path.join(outdir, f'{sample_name}_{segment}.multi.fa')
    lengths: Dict[str, int] = {}

    with tempfile.TemporaryDirectory(prefix=f"metrics_{sample_name}_{segment}_") as tmpdir:
        s3_client = boto3.client("s3")

        with open(output_file, 'w') as out:
            for row in segment_data:
                source = (row.get('source') or "").strip()
                assembly = (row.get('assembly') or "").strip()
                if not source or not assembly:
                    continue

                local_assembly = assembly
                if is_s3_uri(assembly):
                    local_assembly = download_s3_uri(assembly, tmpdir, s3_client)

                sequence_count = 0
                for record in screed.open(local_assembly):
                    sequence_count += 1
                    if sequence_count > 1:
                        sys.exit(f'Error: Each file must contain exactly one sequence ({assembly})')

                    out.write(f">{source}\n{record.sequence}\n")
                    lengths[source] = len(record.sequence)

    return output_file, lengths


def save_metrics_to_csv(metrics_data: List[Dict[str, Any]], filename: str) -> None:
    if not metrics_data:
        logging.warning(f"No data to save for {filename}")
        return

    # stable column order: union of keys, but keep common keys early
    preferred = [
        "sample", "segment", "seq1", "seq2",
        "length_s1", "length_s2",
        "match", "match_iupac", "sub", "ins", "del",
        "ins_terminal", "ins_internal", "del_terminal", "del_internal",
        "valid", "invalid", "valid_s1", "invalid_s1", "valid_s2", "invalid_s2",
        "completeness", "accuracy", "precision"
    ]
    all_keys = set()
    for r in metrics_data:
        all_keys.update(r.keys())
    fieldnames = [k for k in preferred if k in all_keys] + sorted([k for k in all_keys if k not in preferred])

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics_data)

    logging.info(f'Data saved to {filename}')


def _nan_stats(values: List[float]) -> Dict[str, float]:
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return dict(min=float("nan"), max=float("nan"), median=float("nan"), mean=float("nan"), stdev=float("nan"), n=0)
    return dict(
        min=float(np.min(arr)),
        max=float(np.max(arr)),
        median=float(np.median(arr)),
        mean=float(np.mean(arr)),
        stdev=float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
        n=int(arr.size),
    )


def _as_float(x) -> float:
    try:
        if x is None:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")

def summarize_and_report_metrics(
    all_metrics: List[Dict[str, Any]],
    args,
    summary_csv: str = "summary.csv",
    outcsv: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Summarize metrics focusing on accuracy and precision rows separately.
    
    Requires helper functions in your script:
      - _as_float(x) -> float (returns NaN on missing/bad)
      - _nan_stats(list_of_floats) -> dict with keys: n,min,max,median,mean,stdev
    """

    # Separate tracking for accuracy and precision rows
    accuracy_rows = {
        'completeness': [],
        'accuracy': [],
        'n_complete': 0,
        'n_pass': 0
    }
    
    precision_rows = {
        'completeness': [],
        'precision': [],
        'n_complete': 0,
        'n_pass': 0
    }

    failing_rows: List[Dict[str, Any]] = []

    min_comp = getattr(args, "min_completeness", None)
    min_acc  = getattr(args, "min_accuracy", None)
    min_prec = getattr(args, "min_precision", None)

    # Process each row
    for r in all_metrics:
        comp = _as_float(r.get("completeness"))
        acc  = _as_float(r.get("accuracy"))
        prec = _as_float(r.get("precision"))

        fails_threshold = False

        # Process accuracy rows
        if not math.isnan(acc):
            accuracy_rows['accuracy'].append(acc)
            
            if not math.isnan(comp):
                accuracy_rows['completeness'].append(comp)
                
                # Check if completeness meets threshold
                if min_comp is not None and comp >= min_comp:
                    accuracy_rows['n_complete'] += 1
                    
                    # Check accuracy threshold only if completeness passes
                    if min_acc is not None and acc >= min_acc:
                        accuracy_rows['n_pass'] += 1
                    elif min_acc is not None:
                        fails_threshold = True
                elif min_comp is not None:
                    fails_threshold = True

        # Process precision rows
        if not math.isnan(prec):
            precision_rows['precision'].append(prec)
            
            if not math.isnan(comp):
                precision_rows['completeness'].append(comp)
                
                # Check if completeness meets threshold
                if min_comp is not None and comp >= min_comp:
                    precision_rows['n_complete'] += 1
                    
                    # Check precision threshold only if completeness passes
                    if min_prec is not None and prec >= min_prec:
                        precision_rows['n_pass'] += 1
                    elif min_prec is not None:
                        fails_threshold = True
                elif min_comp is not None:
                    fails_threshold = True

        if fails_threshold:
            failing_rows.append(r)

    # Calculate statistics
    acc_comp_stats = _nan_stats(accuracy_rows['completeness'])
    acc_stats = _nan_stats(accuracy_rows['accuracy'])
    
    prec_comp_stats = _nan_stats(precision_rows['completeness'])
    prec_stats = _nan_stats(precision_rows['precision'])

    # Print summary
    print("\n" + "="*60)
    print("ACCURACY ROWS")
    print("="*60)
    print(f"Completeness: n={acc_comp_stats['n']}  "
          f"mean={acc_comp_stats['mean']:.4f} ± {acc_comp_stats['stdev']:.4f}  "
          f"median={acc_comp_stats['median']:.4f}  "
          f"range=[{acc_comp_stats['min']:.4f}, {acc_comp_stats['max']:.4f}]")
    
    print(f"Accuracy:     n={acc_stats['n']}  "
          f"mean={acc_stats['mean']:.4f} ± {acc_stats['stdev']:.4f}  "
          f"median={acc_stats['median']:.4f}  "
          f"range=[{acc_stats['min']:.4f}, {acc_stats['max']:.4f}]")

    if min_comp is not None:
        comp_pass_pct = (accuracy_rows['n_complete'] / acc_comp_stats['n'] * 100) if acc_comp_stats['n'] > 0 else 0.0
        print(f"\nCompleteness ≥{min_comp}: {accuracy_rows['n_complete']}/{acc_comp_stats['n']} ({comp_pass_pct:.1f}%)")
    
    if min_acc is not None and accuracy_rows['n_complete'] > 0:
        acc_pass_pct = (accuracy_rows['n_pass'] / accuracy_rows['n_complete'] * 100)
        print(f"Accuracy ≥{min_acc} (among complete): {accuracy_rows['n_pass']}/{accuracy_rows['n_complete']} ({acc_pass_pct:.1f}%)")

    print("\n" + "="*60)
    print("PRECISION ROWS")
    print("="*60)
    print(f"Completeness: n={prec_comp_stats['n']}  "
          f"mean={prec_comp_stats['mean']:.4f} ± {prec_comp_stats['stdev']:.4f}  "
          f"median={prec_comp_stats['median']:.4f}  "
          f"range=[{prec_comp_stats['min']:.4f}, {prec_comp_stats['max']:.4f}]")
    
    print(f"Precision:    n={prec_stats['n']}  "
          f"mean={prec_stats['mean']:.4f} ± {prec_stats['stdev']:.4f}  "
          f"median={prec_stats['median']:.4f}  "
          f"range=[{prec_stats['min']:.4f}, {prec_stats['max']:.4f}]")

    if min_comp is not None:
        comp_pass_pct = (precision_rows['n_complete'] / prec_comp_stats['n'] * 100) if prec_comp_stats['n'] > 0 else 0.0
        print(f"\nCompleteness ≥{min_comp}: {precision_rows['n_complete']}/{prec_comp_stats['n']} ({comp_pass_pct:.1f}%)")
    
    if min_prec is not None and precision_rows['n_complete'] > 0:
        prec_pass_pct = (precision_rows['n_pass'] / precision_rows['n_complete'] * 100)
        print(f"Precision ≥{min_prec} (among complete): {precision_rows['n_pass']}/{precision_rows['n_complete']} ({prec_pass_pct:.1f}%)")

    # Report failures
    if failing_rows:
        print("\n" + "="*60)
        print("BELOW-THRESHOLD COMPARISONS")
        print("="*60)
        for r in failing_rows:
            sample = r.get("sample", "")
            seg = r.get("segment", "")
            s1 = r.get("seq1", "")
            s2 = r.get("seq2", "")
            comp = _as_float(r.get("completeness"))
            acc  = _as_float(r.get("accuracy"))
            prec = _as_float(r.get("precision"))
            
            comp_str = f"{comp:.4f}" if not math.isnan(comp) else "N/A"
            acc_str = f"{acc:.4f}" if not math.isnan(acc) else "N/A"
            prec_str = f"{prec:.4f}" if not math.isnan(prec) else "N/A"
            
            print(f"{sample}\t{seg}\t{s1} vs {s2}\tcomp={comp_str}\tacc={acc_str}\tprec={prec_str}")
    else:
        print("\n✓ All comparisons meet thresholds")

    # Write summary CSV
    def _fmt(x: Any) -> str:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return ""
        return f"{x:.4f}" if isinstance(x, float) else str(x)

    rows = [
        {
            "row_type": "accuracy",
            "metric": "completeness",
            "n": acc_comp_stats['n'],
            "mean": _fmt(acc_comp_stats['mean']),
            "stdev": _fmt(acc_comp_stats['stdev']),
            "median": _fmt(acc_comp_stats['median']),
            "min": _fmt(acc_comp_stats['min']),
            "max": _fmt(acc_comp_stats['max']),
            "threshold": _fmt(min_comp),
            "n_pass": accuracy_rows['n_complete'] if min_comp else "",
            "pass_pct": f"{accuracy_rows['n_complete']/acc_comp_stats['n']*100:.1f}" if min_comp and acc_comp_stats['n'] > 0 else ""
        },
        {
            "row_type": "accuracy",
            "metric": "accuracy",
            "n": acc_stats['n'],
            "mean": _fmt(acc_stats['mean']),
            "stdev": _fmt(acc_stats['stdev']),
            "median": _fmt(acc_stats['median']),
            "min": _fmt(acc_stats['min']),
            "max": _fmt(acc_stats['max']),
            "threshold": _fmt(min_acc),
            "n_pass": accuracy_rows['n_pass'] if min_acc else "",
            "pass_pct": f"{accuracy_rows['n_pass']/accuracy_rows['n_complete']*100:.1f}" if min_acc and accuracy_rows['n_complete'] > 0 else ""
        },
        {
            "row_type": "precision",
            "metric": "completeness",
            "n": prec_comp_stats['n'],
            "mean": _fmt(prec_comp_stats['mean']),
            "stdev": _fmt(prec_comp_stats['stdev']),
            "median": _fmt(prec_comp_stats['median']),
            "min": _fmt(prec_comp_stats['min']),
            "max": _fmt(prec_comp_stats['max']),
            "threshold": _fmt(min_comp),
            "n_pass": precision_rows['n_complete'] if min_comp else "",
            "pass_pct": f"{precision_rows['n_complete']/prec_comp_stats['n']*100:.1f}" if min_comp and prec_comp_stats['n'] > 0 else ""
        },
        {
            "row_type": "precision",
            "metric": "precision",
            "n": prec_stats['n'],
            "mean": _fmt(prec_stats['mean']),
            "stdev": _fmt(prec_stats['stdev']),
            "median": _fmt(prec_stats['median']),
            "min": _fmt(prec_stats['min']),
            "max": _fmt(prec_stats['max']),
            "threshold": _fmt(min_prec),
            "n_pass": precision_rows['n_pass'] if min_prec else "",
            "pass_pct": f"{precision_rows['n_pass']/precision_rows['n_complete']*100:.1f}" if min_prec and precision_rows['n_complete'] > 0 else ""
        }
    ]

    with open(summary_csv, "w", newline="") as f:
        fieldnames = ["row_type", "metric", "n", "mean", "stdev", "median", "min", "max", "threshold", "n_pass", "pass_pct"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"\n✓ Summary written to: {summary_csv}")
    if outcsv:
        print(f"✓ Combined metrics written to: {outcsv}")

    return {
        "summary_csv": summary_csv,
        "failing_rows": failing_rows,
        "stats": {
            "accuracy_completeness": acc_comp_stats,
            "accuracy": acc_stats,
            "precision_completeness": prec_comp_stats,
            "precision": prec_stats,
        },
        "pass_rates": {
            "accuracy_complete": (accuracy_rows['n_complete'], acc_comp_stats['n']),
            "accuracy_pass": (accuracy_rows['n_pass'], accuracy_rows['n_complete']),
            "precision_complete": (precision_rows['n_complete'], prec_comp_stats['n']),
            "precision_pass": (precision_rows['n_pass'], precision_rows['n_complete']),
        },
    }



# --------------------------
# Main
# --------------------------

def main():
    """Main execution function."""
    version = "1.4"
    parser = argparse.ArgumentParser(
        description="Calculate sequence similarity metrics between aligned sequences."
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Path to samplesheet file (CSV with sample, segment, source, assembly columns)")
    parser.add_argument("--outdir", type=str, default='./', help="Output directory")
    parser.add_argument("--outcsv", type=str, default="all_metrics.csv",
                        help="Combined output CSV filename (written into --outdir unless absolute path).")

    parser.add_argument("--min-completeness", type=float, default=0.8,
                        help="Report to screen any comparisons with completeness below this value.")
    parser.add_argument("--min-accuracy", type=float, default=0.995,
                        help="Report to screen any truth comparisons with accuracy below this value.")
    parser.add_argument("--min-precision", type=float, default=0.995,
                        help="Report to screen any non-truth comparisons with precision below this value.")
    parser.add_argument("--ignore-termini", action="store_true",
                        help="Ignore indels on termini when calculating accuracy / precision.")
    parser.add_argument("--test", action="store_true",
                        help="Process only first sample+segment combination then stop.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {version}')
    args = parser.parse_args()

    print(textwrap.dedent(f"""\
        validate.py v{version}
        ----------------------
    """), flush=True)

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    outcsv = args.outcsv
    if not os.path.isabs(outcsv):
        outcsv = os.path.join(outdir, outcsv)

    # Load samplesheet and group by (sample, segment)
    samplesheet = load_csv_data(args.input)
    segment_groups = defaultdict(list)

    for record in samplesheet:
        sample = (record.get('sample') or '').strip()
        segment = (record.get('segment') or '').strip()
        if sample and segment:
            segment_groups[(sample, segment)].append(record)

    all_metrics: List[Dict[str, Any]] = []

    # Process each sample+segment group
    for count, ((sample_name, segment), segment_data) in enumerate(segment_groups.items(), 1):
        print(f"{sample_name} {segment} has {len(segment_data)} assemblies")

        # Create multi-FASTA and align sequences
        multi_fasta_file, lengths = create_multi_fasta(sample_name, segment, segment_data, outdir)
        aligned_fasta_file = run_mafft_alignment(sample_name, segment, multi_fasta_file, outdir)
        aligned_records = read_alignment_file(aligned_fasta_file)

        # Calculate all pairwise metrics
        pairwise_metrics = []
        for seq1, seq2 in itertools.combinations(aligned_records, 2):
            metrics = calculate_sequence_metrics(seq1, seq2, args.ignore_termini)
            metrics['sample'] = sample_name
            metrics['segment'] = segment
            metrics['length_s1'] = lengths.get(seq1.name)
            metrics['length_s2'] = lengths.get(seq2.name)
            pairwise_metrics.append(metrics)

        if pairwise_metrics:
            all_metrics.extend(pairwise_metrics)
        else:
            logging.warning(f"No pairwise comparisons computed for {sample_name} {segment}")

        if args.test and count == 1:
            break

    # Save metrics to file
    save_metrics_to_csv(all_metrics, outcsv)

    # Generate summary
    summarize_and_report_metrics(all_metrics, args)


if __name__ == "__main__":
    main()