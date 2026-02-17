#!/usr/bin/env python3
# vaper_summary.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, TextIO, Union

from vaper_utils import logging_config

LOGGER = logging_config()
Number = Union[int, float]

# -------------------------------
#  FILE I/O
# -------------------------------

def _open_maybe_gzip(path: str) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

def load_json(path: str) -> Any:
    """Load JSON from file, handling gzip compression."""
    with _open_maybe_gzip(path) as f:
        return json.load(f)

def iter_jsonl(path: str) -> Iterable[Mapping[str, object]]:
    """Iterate over JSONL file, skipping invalid lines."""
    with _open_maybe_gzip(path) as f:
        for line in f:
            if line := line.strip():
                try:
                    yield json.loads(line)
                except Exception as e:
                    LOGGER.warning(f"Skipping bad JSON line in {path}: {e}")

# -------------------------------
#  PARSING
# -------------------------------

def samtoolstats_to_dict(stats_file: str) -> Dict[str, Number]:
    """Extract key metrics from samtools stats output."""
    targets = {
        "reads mapped:": ("reads_mapped", int),
        "bases mapped:": ("bases_mapped", int),
        "average length:": ("mean_mapped_read_length", float),
        "average quality:": ("mean_mapped_read_quality", float),
    }
    out = {}
    needed = {key for key, _ in targets.values()}
    
    with _open_maybe_gzip(stats_file) as fh:
        for line in fh:
            if not line.startswith("SN"):
                continue
            for label, (key, caster) in targets.items():
                if label in line:
                    val_str = line.split(":", 1)[1].split("#", 1)[0].strip()
                    try:
                        out[key] = caster(val_str)
                        needed.discard(key)
                        if not needed:
                            return out
                    except Exception:
                        LOGGER.warning(f"Failed to cast '{val_str}' for {key} in {stats_file}")
    return out

def load_assembly_stats(path: str) -> Optional[Dict[str, object]]:
    """Load assembly stats JSON, returning flattened structure."""
    try:
        rec = load_json(path)
    except Exception as e:
        LOGGER.error(f"Failed to parse assembly stats JSON {path}: {e}")
        return None
    
    if not isinstance(rec, dict) or "ref" not in rec:
        LOGGER.warning(f"Invalid assembly stats format: {path}")
        return None
    
    return {
        "ref": rec["ref"],
        "assembly": {k: v for k, v in rec.items() if k not in ("ref", "query")}
    }

def extract_ref_from_filename(filename: str, sample: str) -> str:
    """Extract reference name from BAM stats filename."""
    pattern = f"{sample}-(.+?)\\.stats\\.txt(?:\\.gz)?$"
    if match := re.search(pattern, filename):
        return match.group(1)
    raise ValueError(f"Could not infer reference from filename: {filename}")

# -------------------------------
#  DATA STRUCTURE HELPERS
# -------------------------------

def get_or_create_assembly(data: Dict[str, Any], ref: str) -> Dict[str, Any]:
    """Find or create an assemblies entry for the given reference."""
    assemblies = data.setdefault("assemblies", [])
    for entry in assemblies:
        if isinstance(entry, dict):
            if entry.get("reference") == ref:
                return entry
    
    new_entry = {"reference": ref}
    assemblies.append(new_entry)
    return new_entry

def get_ref_length(entry: Dict[str, Any]) -> Optional[float]:
    """Extract reference length from assembly or reference metadata."""
    length_keys = ("ref_len", "rlen", "length", "genome_length", "reference_length", "size", "bp")
    
    for section in ("assembly", "reference"):
        if block := entry.get(section):
            if isinstance(block, dict):
                for key in length_keys:
                    if isinstance(val := block.get(key), (int, float)):
                        return float(val)
    return None

def calculate_qc_status(assembly: Dict[str, Any], gf_min: Optional[float], depth_min: Optional[float]) -> List[Dict]:
    """Calculate QC status for genome fraction and read depth."""
    gf = assembly.get('assembly', {}).get('genome_fraction') or assembly.get('assembly', {}).get('gf')
    depth = assembly.get('alignment', {}).get('read_depth')
    
    return [
        {
            'key': 'genome_fraction',
            'value': gf,
            'threshold': gf_min,
            'status': gf >= gf_min if (gf is not None and gf_min is not None) else None
        },
        {
            'key': 'read_depth',
            'value': depth,
            'threshold': depth_min,
            'status': depth >= depth_min if (depth is not None and depth_min is not None) else None
        }
    ]

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    version = "1.7.0"
    parser = argparse.ArgumentParser(description="Summarize VAPER run outputs to a single JSON.")
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--fastp", help="Path to fastp summary (.json[.gz]).")
    parser.add_argument("--metagenome", help="Path to Sourmash metagenome summary (JSON).")
    parser.add_argument("--ref_info", nargs="+", help="Path(s) to reference metadata file(s) (JSONL[.gz]).")
    parser.add_argument("--bam_stats", nargs="+", help="Path(s) to samtools stats file(s) named '<sample>-<ref>.stats.txt[.gz]'.")
    parser.add_argument("--assembly_stats", nargs="+", help="Path(s) to assembly stats JSON file(s).")
    parser.add_argument("--max_depth", type=int, help="Floor for mapped reads used in depth calc.")
    parser.add_argument("--qc_gf_min", type=float, help="Minimum genome fraction to pass QC.")
    parser.add_argument("--qc_depth_min", type=float, help="Minimum estimated depth to pass QC.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
    parser.add_argument("--version", action="version", version=version)
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info("Author: Jared Johnson")

    data = {"reads": {"raw": {}, "filtered": {}}, "metagenome": [], "assemblies": []}
    filtered_total_bases = None

    # Load fastp summary
    if args.fastp:
        LOGGER.info(f"Loading fastp summary: {args.fastp}")
        summary = load_json(args.fastp).get("summary", {})
        data["reads"] = {
            "raw": summary.get("before_filtering", {}),
            "filtered": summary.get("after_filtering", {})
        }
        filtered_total_bases = data["reads"]["filtered"].get("total_bases")

    # Load metagenome table
    if args.metagenome:
        LOGGER.info(f"Loading metagenome table: {args.metagenome}")
        try:
            data["metagenome"] = load_json(args.metagenome)
        except Exception as e:
            LOGGER.error(f"Failed to load metagenome table: {e}")

    # Load assembly stats
    if args.assembly_stats:
        for filepath in args.assembly_stats:
            LOGGER.info(f"Loading assembly stats: {filepath}")
            if rec := load_assembly_stats(filepath):
                entry = get_or_create_assembly(data, str(rec["ref"]))
                entry["assembly"] = {**entry.get("assembly", {}), **rec["assembly"]}

    # Process BAM stats
    if args.bam_stats:
        for filepath in args.bam_stats:
            ref = extract_ref_from_filename(os.path.basename(filepath), args.sample)
            LOGGER.info(f"Parsing samtools stats for ref={ref}: {filepath}")
            
            stats = samtoolstats_to_dict(filepath)
            entry = get_or_create_assembly(data, ref)
            
            bases_mapped = int(stats.get("bases_mapped", 0))
            mean_len = float(stats.get("mean_mapped_read_length", 0.0))
            effective_bases = max(bases_mapped, args.max_depth or 0)
            
            # Calculate metrics
            pct_mapped = round(100.0 * bases_mapped / filtered_total_bases, 1) if filtered_total_bases else None
            
            read_depth = None
            if (ref_len := get_ref_length(entry)) and ref_len > 0 and mean_len > 0:
                depth = effective_bases / ref_len
                read_depth = min(int(depth), args.max_depth) if args.max_depth else int(depth)
            
            entry["alignment"] = {
                **entry.get("alignment", {}),
                **stats,
                "max_depth": args.max_depth,
                "pct_bases_mapped": pct_mapped,
                "read_depth": read_depth
            }

    # Apply QC thresholds
    for assembly in data.get('assemblies', []):
        if isinstance(assembly, dict):
            assembly['qc'] = calculate_qc_status(assembly, args.qc_gf_min, args.qc_depth_min)

    # Create top-level entry for sample
    data_out = {str(args.sample): data}

    # Append reference metadata
    if args.ref_info:
        ref_list = [a.get('reference') for a in data.get('assemblies', [{}]) ]

        ref_block = {}
        for filepath in args.ref_info:
            LOGGER.info(f"Loading reference metadata JSONL: {filepath}")
            for rec in iter_jsonl(filepath):
                if (ref_name := rec.get("name")) and (ref_name in ref_list):
                    ref_block[ref_name] = { k: v for k, v in rec.items() if k != 'name'}

        data_out['references'] = ref_block

    # Write output
    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, f"{args.sample}.json.gz")

    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        json.dump(
            data_out,
            out,
            separators=(",", ":"),
            sort_keys=True
        )
        out.write("\n")  # POSIX newline

    LOGGER.info(f"Wrote summary: {out_path}")

if __name__ == "__main__":
    main()