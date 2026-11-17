"""Readers for legacy and Shiba 1.0 bulk output tables."""
import csv
import math
from pathlib import Path


def load_experiment_table(experiment_table):
    """Read named columns, retaining additional Shiba metadata such as technology."""
    with open(experiment_table, newline="") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        header_line = next(lines, "")
        # Retain support for the whitespace-delimited tables accepted by v0.1.
        split = (lambda line: next(csv.reader([line], delimiter="\t"))) if "\t" in header_line else str.split
        header = split(header_line.rstrip("\r\n"))
        if not {"sample", "bam", "group"}.issubset(header):
            raise ValueError("Experiment table must contain sample, bam, and group columns (bulk RNA-seq)")
        result = {}
        for number, line in enumerate(lines, 2):
            values = split(line.rstrip("\r\n"))
            if len(values) != len(header):
                raise ValueError(f"Experiment table row {number}: expected {len(header)} columns")
            row = dict(zip(header, values))
            sample = row.pop("sample")
            if not sample or sample in result:
                raise ValueError(f"Missing or duplicate sample name: {sample!r}")
            if not row["bam"] or not row["group"]:
                raise ValueError(f"Missing bam or group for sample {sample}")
            result[sample] = row
    if not result:
        raise ValueError("Experiment table has no samples")
    return result


def get_psi_values(positional_id, shiba_path):
    """Read sample PSI as percentages; preserve empty/trailing cells as NA."""
    path = Path(shiba_path) / "results/splicing/PSI_matrix_sample.txt"
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"event_id", "pos_id"}.issubset(reader.fieldnames):
            raise ValueError(f"PSI matrix must contain event_id and pos_id columns: {path}")
        samples = [key for key in reader.fieldnames if key not in ("event_id", "pos_id")]
        result = dict.fromkeys(samples, "NA")
        for row in reader:
            if row["pos_id"] == positional_id:
                for sample in samples:
                    try:
                        value = float(row[sample])
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(value) and 0 <= value <= 1:
                        result[sample] = value * 100
                break
    return result
