#!/usr/bin/env python3

# combine_summary.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import gzip
import json
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, List, Any, Tuple

from vaper_utils import logging_config

LOGGER = logging_config()
VERSION = "1.0"


def open_maybe_gzip(path: str):
    """Open file, handling gzip compression automatically."""
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def load_json(path: str) -> Dict[str, Any]:
    """Load JSON file with error handling and logging."""
    LOGGER.debug(f"Loading JSON file: {path}")
    try:
        with open_maybe_gzip(path) as f:
            data = json.load(f)
        LOGGER.debug(f"Successfully loaded {path}")
        return data
    except json.JSONDecodeError as e:
        LOGGER.error(f"Invalid JSON in file {path}: {e}")
        raise
    except FileNotFoundError:
        LOGGER.error(f"File not found: {path}")
        raise
    except Exception as e:
        LOGGER.error(f"Error loading {path}: {e}")
        raise


def extract_reference_info(ref_name: str, ref_data: Dict[str, Any]) -> Tuple[str, str]:
    """Extract species and segment information from reference data."""
    if ref_name not in ref_data:
        LOGGER.debug(f"Reference '{ref_name}' not found in reference data")
        return "", ""
    
    ref_rec = ref_data[ref_name]
    species = ref_rec.get('species') or ref_rec.get('metadata', {}).get('species') or []
    segment = ref_rec.get('segment') or ref_rec.get('metadata', {}).get('segment') or []
    
    species_str = '; '.join(species) if isinstance(species, list) else species
    segment_str = '; '.join(segment) if isinstance(segment, list) else segment
    
    return species_str, segment_str


def process_qc_tests(qc_tests: List[Dict[str, Any]]) -> Tuple[str, str]:
    """Process QC tests and return status and reasons."""
    if not qc_tests:
        return "FAIL", "no QC tests"
    
    qc_status = all(test.get('status') for test in qc_tests)
    reasons = [
        f"{test.get('key')} < {test.get('threshold')}" 
        for test in qc_tests 
        if not test.get('status')
    ]
    
    status = "PASS" if qc_status else "FAIL"
    reason = '; '.join(reasons) if reasons else ""
    
    return status, reason


def process_metagenome_summary(metagenome: Dict[str, Any]) -> str:
    """Generate species summary from metagenome data."""
    species_list = list(metagenome.get("species") or [])
    frac_list = list(metagenome.get("fraction") or [])
    
    if not species_list or not frac_list:
        return ""
    
    pairs = sorted(
        zip(species_list, frac_list), 
        key=lambda x: (-float(x[1] or 0), str(x[0] or ""))
    )
    
    return '; '.join(
        f"{round(100 * float(frac or 0), 1)}% {sp}" 
        for sp, frac in pairs
    )


def to_rows(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Convert a result object to one or more flattened row dictionaries."""
    sample_id, values, ref_data = None, {}, {}
    
    for k, v in result.items():
        if k == 'references':
            ref_data = v
        else:
            sample_id, values = k, v
    
    if not sample_id:
        LOGGER.warning("Result object missing sample ID")
        return []
    
    LOGGER.debug(f"Processing sample: {sample_id}")
    
    row = {'id': sample_id}

    # Process read metrics
    reads_data = values.get("reads") or {}
    for stage, metrics in reads_data.items():
        if isinstance(metrics, dict):
            for k, v in metrics.items():
                row[f"{stage}_{k}"] = v
    
    if reads_data:
        LOGGER.debug(f"  Processed {len(reads_data)} read stages")

    # Process metagenome species summary
    metagenome = values.get("metagenome") or {}
    row['species_summary'] = process_metagenome_summary(metagenome)

    # Handle missing assemblies
    assemblies = values.get("assemblies")
    if not assemblies:
        LOGGER.info(f"  No assemblies found for {sample_id}")
        row['reference'] = "no reference"
        return [row]

    # Process each assembly
    LOGGER.debug(f"  Processing {len(assemblies)} assemblies")
    data = []
    
    for idx, assembly in enumerate(assemblies, 1):
        a_row = row.copy()
        
        # Process assembly and alignment metrics
        for metric_type in ['assembly', 'alignment']:
            metrics = assembly.get(metric_type) or {}
            for k, v in metrics.items():
                if k == 'covered_non_N_on_ref':
                    continue
                
                # Rename specific keys
                key_map = {
                    'qlen': 'assembly_length', 
                    'rlen': 'reference_length'
                }
                new_key = key_map.get(k, f"assembly_{k}")
                a_row[new_key] = v
        
        # Process QC status
        qc_tests = assembly.get('qc')
        if qc_tests is not None:
            status, reason = process_qc_tests(qc_tests)
            a_row['assembly_qc'] = status
            a_row['assembly_qc_reason'] = reason
            LOGGER.debug(f"    Assembly {idx}: QC {status}")
        
        # Process reference information
        ref_name = assembly.get("reference")
        if ref_name:
            a_row['reference'] = ref_name
            species, segment = extract_reference_info(ref_name, ref_data)
            a_row['species'] = species
            a_row['segment'] = segment
        
        data.append(a_row)
    
    return data


def add_assembly_variant_numbering(rows: List[Dict[str, Any]]) -> None:
    """Add assembly variant numbering to rows (modifies in place)."""
    groups = defaultdict(list)

    for idx, r in enumerate(rows):
        sample_id = r.get("id")
        species = r.get("species")
        segment = r.get("segment")

        # Skip rows without species/segment info
        if not species or not segment:
            continue

        key = (sample_id, species, segment)
        groups[key].append(idx)

    variant_count = 0
    for idxs in groups.values():
        total = len(idxs)
        if total > 1:
            variant_count += total
        for i, idx in enumerate(idxs, start=1):
            rows[idx]["assembly_variant"] = f"{i} of {total}"
    
    if variant_count > 0:
        LOGGER.info(f"Added variant numbering to {variant_count} assemblies across {len(groups)} groups")


def write_combined_json(params_obj: Dict[str, Any], sample_data: Dict[str, Any], 
                       ref_data: Dict[str, Any], output_path: str = "VAPER-summary.json.gz") -> None:
    """Write combined JSON output."""
    LOGGER.info(f"Writing combined JSON to {output_path}")
    
    combined = {
        "params": params_obj, 
        "results": sample_data, 
        "references": ref_data
    }
    
    try:
        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            json.dump(combined, f, ensure_ascii=False, separators=(",", ":"))
            f.write("\n")
        
        file_size = Path(output_path).stat().st_size
        LOGGER.info(f"Successfully wrote {output_path} ({file_size:,} bytes)")
    except Exception as e:
        LOGGER.error(f"Failed to write {output_path}: {e}")
        raise


def write_csv(rows: List[Dict[str, Any]], output_path: str = "VAPER-summary.csv") -> None:
    """Write flattened CSV output."""
    LOGGER.info(f"Writing CSV to {output_path}")
    
    fieldnames = [
        "id", "species", "segment", "reference", "assembly_qc", "assembly_qc_reason",
        "assembly_variant", "assembly_length", "assembly_read_depth", "assembly_genome_fraction",
        "assembly_identity", "assembly_ins", "assembly_del", "assembly_sub", "assembly_missing",
        "assembly_mixed", "assembly_pct_bases_mapped", "assembly_reads_mapped",
        "assembly_mean_mapped_read_length", "assembly_mean_mapped_read_quality",
        "filtered_q30_rate", "filtered_read1_mean_length", "filtered_read2_mean_length",
        "filtered_total_reads", "raw_q30_rate", "raw_read1_mean_length",
        "raw_read2_mean_length", "raw_total_reads", "species_summary"
    ]

    # Create ordered rows
    ordered_rows = [
        OrderedDict((fn, r.get(fn, "")) for fn in fieldnames) 
        for r in rows
    ]
    
    try:
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore", restval="")
            writer.writeheader()
            writer.writerows(ordered_rows)
        
        LOGGER.info(f"Successfully wrote {len(ordered_rows)} rows to {output_path}")
    except Exception as e:
        LOGGER.error(f"Failed to write {output_path}: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Combine params and results JSON; emit flattened results CSV."
    )
    parser.add_argument("--params", required=True, help="Parameters JSON (optionally .gz)")
    parser.add_argument("--results", nargs="+", required=True, help="Result JSON files (optionally .gz)")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    LOGGER.info(f"VAPER combine_summary v{VERSION}")
    LOGGER.info(f"Author: Jared Johnson")
    LOGGER.info(f"Parameters file: {args.params}")
    LOGGER.info(f"Result files: {len(args.results)}")

    # Load parameters
    params_obj = load_json(args.params)
    LOGGER.info(f"Loaded parameters")

    # Load results
    LOGGER.info(f"Loading {len(args.results)} result files...")
    results_list = []
    for path in args.results:
        try:
            results_list.append(load_json(path))
        except Exception as e:
            LOGGER.error(f"Skipping {path} due to error: {e}")
            continue
    
    LOGGER.info(f"Successfully loaded {len(results_list)} result files")

    # Merge sample and reference data
    sample_data, ref_data = {}, {}
    for i, obj in enumerate(results_list):
        if not isinstance(obj, dict):
            LOGGER.warning(f"Skipping non-dict result at index {i}")
            continue

        for k1, v1 in obj.items():
            if k1 == 'references':
                if not isinstance(v1, dict):
                    continue
                for k2, v2 in v1.items():
                    if k2 not in ref_data:
                        ref_data[k2] = v2
            else:
                sample_data[k1] = v1

    LOGGER.info(f"Merged data: {len(sample_data)} samples, {len(ref_data)} references")

    # Write combined JSON
    write_combined_json(params_obj, sample_data, ref_data)

    # Flatten to rows
    LOGGER.info("Flattening results to rows...")
    rows = []
    for res in results_list:
        rows.extend(to_rows(res))
    
    LOGGER.info(f"Generated {len(rows)} total rows")

    # Add assembly variant numbering
    add_assembly_variant_numbering(rows)

    # Write CSV
    write_csv(rows)

    LOGGER.info("Processing complete")
    print("\nOutput files:")
    print("  - VAPER-summary.json.gz")
    print("  - VAPER-summary.csv")


if __name__ == "__main__":
    main()