#!/usr/bin/env python3

# vaper_metagenome.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import sys
import os
import argparse
import csv
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import json

from vaper_utils import logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

OKABE_ITO_BASE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
]

# -------------------------------
#  FUNCTIONS
# -------------------------------

def okabe_ito_ramp(n: int):
    """Interpolate the Okabe–Ito palette to n colors (reversed)."""
    if n <= 0:
        return []
    cmap = LinearSegmentedColormap.from_list("okabe_ito", OKABE_ITO_BASE[::-1])
    if n == 1:
        return [cmap(0.5)]
    return [cmap(i / (n - 1)) for i in range(n)]

def load_rows(path: str):
    """Return list of dict rows from CSV (or empty list if file missing)."""
    if not os.path.isfile(path):
        LOGGER.warning(f"File not found: {path}")
        return []
    LOGGER.info(f"Loading CSV: {path}")
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        rows = [dict(row) for row in r]
    LOGGER.info(f"Loaded {len(rows)} rows from {path}")
    return rows

def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def lineage_to_species(lineage: str) -> str:
    if lineage is None:
        return ""
    s = str(lineage)
    parts = [p for p in s.split(";") if p != ""]
    return parts[-1] if parts else s

def write_summary_line(prefix: str, data: str):
    out_summary = f"{prefix}.taxa-summary.json"
    LOGGER.info(f"Writing summary to {out_summary}")
    with open(out_summary, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=True)

def make_outputs(rows, prefix: str):
    out_plot = f"{prefix}.taxa-plot.jpg"

    LOGGER.info(f"Processing {len(rows)} total rows")
    species_rows = [r for r in rows if (r.get("rank") == "species")]
    LOGGER.info(f"Found {len(species_rows)} species-level rows")
    
    summaryline =  {"species": [], "fraction": []}
    if not species_rows:
        LOGGER.warning("No species rows detected → writing 'No Viruses Detected'")
        write_summary_line(prefix, summaryline)
        return

    unclass = [r for r in species_rows if r.get("lineage") == "unclassified"]
    unclass_frac = safe_float(unclass[0].get("fraction")) if unclass else 0.0
    denom = 1.0 - unclass_frac
    if denom <= 0:
        denom = 1e-12
        LOGGER.warning("Denominator for classified_fraction <= 0, forcing small epsilon")

    recs = []
    for r in species_rows:
        fraction = safe_float(r.get("fraction"))
        cf = fraction / denom
        label = lineage_to_species(r.get("lineage")) if cf >= 0.01 else "Other"
        recs.append({
            "Species": label,
            "fraction": fraction,
            "classified_fraction": cf
        })

    agg = {}
    for r in recs:
        sp = r["Species"]
        if sp not in agg:
            agg[sp] = {"fraction": 0.0, "classified_fraction": 0.0}
        agg[sp]["fraction"] += r["fraction"]
        agg[sp]["classified_fraction"] += r["classified_fraction"]

    items = sorted(agg.items(), key=lambda kv: kv[1]["fraction"], reverse=True)
    plot_items = [(sp, vals) for sp, vals in items if sp != "unclassified"]

    if not plot_items:
        LOGGER.warning("Only unclassified taxa found → writing 'No Viruses Detected'")
        write_summary_line(prefix, "No Viruses Detected")
        return

    labels = [sp for sp, _ in plot_items]
    sizes = [vals["classified_fraction"] * 100.0 for _, vals in plot_items]
    colors = okabe_ito_ramp(len(labels))

    LOGGER.info(f"Generating donut chart with {len(labels)} labels")
    fig, ax = plt.subplots(figsize=(15, 10))
    wedges, _ = ax.pie(
        sizes,
        labels=None,
        startangle=90,
        colors=colors,
        wedgeprops=dict(width=0.6)
    )
    ax.set(aspect="equal")
    ax.legend(wedges, labels, title="Species", loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    fig.savefig(out_plot, dpi=300)
    plt.close(fig)
    LOGGER.info(f"Saved plot: {out_plot}")

    sum_items = sorted(plot_items, key=lambda kv: kv[0])
    for sp, vals in sum_items:
        summaryline['species'].append(sp)
        summaryline['fraction'].append(round(vals['classified_fraction'], 3))
    write_summary_line(prefix, summaryline)

# -------------------------------
#  MAIN
# -------------------------------

def main():
    version = 1.0

    parser = argparse.ArgumentParser()
    parser.add_argument("sm_taxa", help="Path to Sourmash taxa CSV")
    parser.add_argument("prefix", help="Output prefix")
    parser.add_argument("--version", action="version", version=str(version), help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    rows = load_rows(args.sm_taxa)
    make_outputs(rows, args.prefix)
    LOGGER.info("Completed successfully")

if __name__ == "__main__":
    main()
