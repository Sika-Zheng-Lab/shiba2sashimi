"""Read Shiba's merged counts (including the Tosa-based v1.0 format)."""
import csv
from pathlib import Path
from .utils import coord2int, same_chrom


def extract_junctions_in_region(shiba_path, chrom, start, end, junction_list=None):
    """Read splice counts and explicitly requested one-base RI boundary counts."""
    path = Path(shiba_path) / "junctions/junctions.bed"
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chr", "start", "end", "ID"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid Shiba junction table header: {path}")
        samples = [key for key in reader.fieldnames if key not in required]
        result = {sample: {} for sample in samples}
        # Normalize only for matching; retain the event's ID for arc direction lookup.
        targets = None if junction_list is None else {
            (coord2int(j)[0].removeprefix("chr"), *coord2int(j)[1:]): j for j in junction_list
        }
        for row in reader:
            jchrom, left, right = coord2int(row["ID"])
            if not same_chrom(jchrom, chrom):
                continue
            boundary = right - left == 1
            # Never turn arbitrary one-base counts into splice-junction arcs.
            if boundary and targets is None:
                continue
            if not (start <= left <= end if boundary else start <= left < right - 1 <= end):
                continue
            key = (jchrom.removeprefix("chr"), left, right)
            if targets is not None and key not in targets:
                continue
            junction_id = targets[key] if targets is not None else row["ID"]
            for sample in samples:
                count = int(row[sample])
                if count < 0:
                    raise ValueError(f"Negative junction count for {sample}: {junction_id}")
                result[sample][junction_id] = count
    return result
