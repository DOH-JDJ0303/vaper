#!/usr/bin/env python3

# vaper_condense.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import os
import re
import csv
import argparse
import datetime
import shutil
from typing import Dict, List, Set, Tuple
import gzip

import numpy as np
import screed
import sourmash
from sklearn.cluster import DBSCAN

from vaper_utils import logging_config, safe_filename

# -------------------------------
# GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
# HELPERS
# -------------------------------

def dict_to_csv(dicts, filename, cols):
    os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
    with open(filename, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, quoting=csv.QUOTE_ALL)
        w.writeheader()
        w.writerows(dicts)

def compute_distance_matrix(minhashes: Dict[str, Dict[str, object]]) -> Tuple[List[str], np.ndarray]:
    """
    Build a symmetric, precomputed distance matrix using sourmash containment_ani().dist.
    Returns (ids, NxN numpy array).
    """
    ids = sorted(minhashes.keys())
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        m1 = minhashes[ids[i]]['mh']
        for j in range(i + 1, n):
            m2 = minhashes[ids[j]]['mh']
            # containment_ani() returns an object with .dist (0..1 where 0 is identical)
            try:
                d = m1.containment_ani(m2).dist
                if d is None or np.isnan(d):
                    d = 1.0
            except Exception:
                # If comparison fails (e.g., empty sketches), push far apart
                d = 1.0
            mat[i, j] = d
            mat[j, i] = d
    return ids, mat

def labels_to_cluster_ids(labels: np.ndarray) -> Dict[str, int]:
    """
    Convert sklearn DBSCAN labels into contiguous non-negative cluster IDs.
    Noise (-1) becomes its own unique singleton cluster.
    """
    out: Dict[str, int] = {}
    return out

def assign_clusters(ids: List[str], labels: np.ndarray) -> Dict[str, int]:
    """
    Map DBSCAN labels to cluster integers. Noise (-1) gets unique IDs.
    """
    # Existing non-noise labels as-is, but we must ensure noise get unique ids
    next_id = (labels[labels >= 0].max() + 1) if np.any(labels >= 0) else 0
    mapped = {}
    for name, lab in zip(ids, labels):
        if lab == -1:
            mapped[name] = int(next_id)
            next_id += 1
        else:
            mapped[name] = int(lab)
    return mapped

def condense(data: Dict[str, Dict[str, object]]) -> Tuple[List[str], Dict[str, Dict[str, object]]]:
    clusters: Dict[int, List[str]] = {}
    for key, value in data.items():
        clusters.setdefault(value['cluster'], []).append(key)

    refs: List[str] = []
    for members in clusters.values():
        if len(members) > 1:
            max_covXdepth, select_ref = -1.0, None
            for id_ in members:
                try:
                    score = float(data[id_]['covXdepth'])
                except Exception:
                    score = -1.0
                if score > max_covXdepth:
                    select_ref, max_covXdepth = id_, score
        else:
            select_ref = members[0]
        refs.append(select_ref)
        for id_ in members:
            data[id_]['condensed'] = '' if id_ == select_ref else select_ref
    return refs, data

def save_selected_assemblies(refs: List[str], records: Dict[str, Dict[str, str]], outdir: str) -> int:
    os.makedirs(outdir, exist_ok=True)
    written = 0
    for rid in refs:
        rec = records.get(rid)
        if not rec:
            continue
        asm, seq = rec['assembly'], rec['sequence']
        fname = f"{safe_filename(asm)}.fa.gz"
        path = os.path.join(outdir, fname)
        if os.path.exists(path):
            LOGGER.info(f"Exists, not overwriting: {path}")
            continue
        with gzip.open(path, 'wt') as fh:
            fh.write(f">{asm}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
        written += 1
    return written

# -------------------------------
# MAIN
# -------------------------------

def main():
    version = "1.2"

    p = argparse.ArgumentParser()
    p.add_argument("--fasta", nargs="+", required=True,
                   help="One or more FASTA files.")
    p.add_argument("--stats", dest="read_stats", nargs="+", required=True,
                   help="One or more TSVs with #rname, coverage, meandepth. Same header across files.")
    p.add_argument("--dist_threshold", type=float, default=0.05,
                   help="DBSCAN eps on MinHash distance (0..1). Default 0.05")
    p.add_argument("--min_samples", type=int, default=1,
                   help="DBSCAN min_samples. Default 1 (allow singletons).")
    p.add_argument("--ksize", type=int, default=31, help="k-mer size for sourmash (default: 31)")
    p.add_argument("--scaled", type=int, default=10, help="scaled for sourmash MinHash (default: 10)")
    p.add_argument("--outdir", default="./", help="Output directory (CSV + selected FASTAs).")
    p.add_argument("--prefix", default='', help="Prefix for output names.")
    p.add_argument('--version', action='version', version="%(prog)s " + version)
    args = p.parse_args()

    LOGGER.info(f"vaper_condense v{version}")
    os.makedirs(args.outdir, exist_ok=True)

    # Short-circuit: if only one FASTA is provided, just copy/rename it and exit
    if len(args.fasta) == 1:
        fasta = args.fasta[0]
        base = os.path.basename(fasta)
        outpath = os.path.join(args.outdir, base)
        shutil.copy(fasta, outpath)
        LOGGER.info(f"Single FASTA supplied: copied to {outpath}")
        return

    # ---------------- Merge multiple stats ----------------
    read_stats: Dict[str, Dict[str, str]] = {}
    expected_header = None
    for path in args.read_stats:
        with open(path, newline='') as csvfile:
            rdr = csv.reader(csvfile, delimiter='\t')
            header = next(rdr)
            if expected_header is None:
                expected_header = header
                try:
                    idx_name = header.index('#rname')
                    idx_cov = header.index('coverage')
                    idx_depth = header.index('meandepth')
                except ValueError as e:
                    raise SystemExit(f"Stats file {path} missing required columns: {e}")
            else:
                if header != expected_header:
                    raise SystemExit(f"Stats header mismatch in {path}. "
                                     f"Expected: {expected_header}, got: {header}")

            for row in rdr:
                if not row:
                    continue
                name = row[idx_name]
                read_stats[name] = {'cov': row[idx_cov], 'depth': row[idx_depth]}
        LOGGER.info(f"Merged stats from {path}")

    LOGGER.info(f"Total stats entries: {len(read_stats)}")

    # ---------------- Load multiple FASTAs & build MinHashes ----------------
    prefix_re = re.compile(rf"^{re.escape(args.prefix)}_") if args.prefix else None
    minhashes: Dict[str, Dict[str, object]] = {}
    records: Dict[str, Dict[str, str]] = {}  # name -> {'assembly', 'sequence'}

    loaded = 0
    for fasta in args.fasta:
        with screed.open(fasta) as seqs:
            for rec in seqs:
                name = rec.name
                if prefix_re:
                    name = prefix_re.sub('', name)

                mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
                mh.add_sequence(rec.sequence, force=True)

                if name in read_stats:
                    try:
                        cov = float(read_stats[name]['cov']) / 100.0
                        depth = float(read_stats[name]['depth'])
                        covXdepth = cov * depth
                    except Exception:
                        covXdepth = 0.0
                else:
                    LOGGER.warning(f"No stats for {name}; covXdepth=0")
                    covXdepth = 0.0

                minhashes[name] = {
                    'assembly': rec.name,
                    'mh': mh,
                    'covXdepth': covXdepth
                }
                records[name] = {
                    'assembly': rec.name,
                    'sequence': rec.sequence
                }
                loaded += 1

    LOGGER.info(f'{datetime.datetime.now()}: Loaded {loaded} sequences from {len(args.fasta)} FASTA file(s)')

    # ---------------- DBSCAN clustering on MinHash distances ----------------
    LOGGER.info(f'{datetime.datetime.now()}: Building distance matrix for DBSCAN...')
    ids, dist_mat = compute_distance_matrix(minhashes)
    LOGGER.info(f'Distance matrix shape: {dist_mat.shape}')

    LOGGER.info(f'Clustering with DBSCAN (eps={args.dist_threshold}, min_samples={args.min_samples})')
    db = DBSCAN(eps=args.dist_threshold, min_samples=args.min_samples, metric='precomputed')
    labels = db.fit_predict(dist_mat)

    # Map labels to cluster IDs; give noise its own singleton ids
    cluster_map = assign_clusters(ids, labels)

    for n in minhashes:
        minhashes[n]['cluster'] = cluster_map.get(n, -1)

    n_clusters = len(set(cluster_map.values()))
    LOGGER.info(f'{datetime.datetime.now()}: {n_clusters} clusters formed.')

    # ---------------- Condense ----------------
    if n_clusters == len(minhashes):
        LOGGER.info(f'{datetime.datetime.now()}: Nothing to condense (all singletons).')
        for v in minhashes.values():
            v['condensed'] = ''
        refs = list(minhashes.keys())
    else:
        LOGGER.info(f'{datetime.datetime.now()}: Condensing references by cov×depth within clusters.')
        refs, minhashes = condense(minhashes)

    # ---------------- Write CSV (named by prefix) ----------------
    out_csv = os.path.join(args.outdir, f"{args.prefix}-condensed.csv" if args.prefix else "condensed.csv")
    rows = []
    for _, v in minhashes.items():
        out = dict(v)
        out.pop('mh', None)
        out.pop('cluster', None)
        rows.append(out)

    dict_to_csv(rows, out_csv, ["assembly", "condensed", "covXdepth"])
    LOGGER.info(f"Wrote {len(rows)} rows -> {out_csv}")

    # ---------------- Save selected assemblies (no overwrite) ----------------
    sel_dir = os.path.join(args.outdir)  # save directly to outdir as requested
    written = save_selected_assemblies(refs, records, sel_dir)
    LOGGER.info(f"Selected assemblies written: {written} (skipped existing files remain untouched)")

if __name__ == "__main__":
    main()
