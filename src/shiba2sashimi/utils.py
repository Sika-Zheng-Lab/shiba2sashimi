"""Event coordinates. Plotting uses zero-based, half-open intervals."""
import csv
from pathlib import Path

EVENT_TYPES = ("SE", "FIVE", "THREE", "MXE", "RI", "AFE", "ALE", "MSE")


def normalize_chrom(chrom):
    # Match Shiba's chromosome normalization, including short contig names.
    return f"chr{chrom}" if chrom.isdecimal() or len(chrom) <= 2 else chrom


def same_chrom(left, right):
    return normalize_chrom(left) == normalize_chrom(right)


def coord2int(coordinate):
    """Parse an explicit zero-based, half-open plotting region."""
    try:
        chrom, positions = coordinate.rsplit(":", 1)
        start, end = map(int, positions.split("-"))
        if not chrom or start < 0 or end <= start:
            raise ValueError
    except (ValueError, AttributeError):
        raise ValueError(f"Invalid coordinate {coordinate!r}; expected chr:start-end with 0 <= start < end") from None
    return chrom, start, end


def read_event(positional_id, shiba_path):
    event_type = positional_id.split("@")[0]
    if event_type not in EVENT_TYPES:
        raise ValueError(f"Unsupported event type: {event_type}")
    path = Path(shiba_path) / "events" / f"EVENT_{event_type}.txt"
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"pos_id", "strand"}.issubset(reader.fieldnames):
            raise ValueError(f"Missing pos_id or strand in event table: {path}")
        for row in reader:
            if row["pos_id"] == positional_id:
                return row
    raise ValueError(f"Positional ID not found in EVENT file: {positional_id}")


def event_junctions(event):
    """Directions follow event definitions, never genomic ordering of exon A/B."""
    kind = event["pos_id"].split("@")[0]
    keys = {
        "SE": {"intron_c": "up", "intron_a": "down", "intron_b": "down"},
        "FIVE": {"intron_b": "up", "intron_a": "down"},
        "THREE": {"intron_b": "up", "intron_a": "down"},
        "MXE": {"intron_b1": "up", "intron_b2": "up", "intron_a1": "down", "intron_a2": "down"},
        "RI": {"intron_a": "up"},
        "AFE": {"intron_a": "up", "intron_b": "down"},
        "ALE": {"intron_a": "up", "intron_b": "down"},
        "MSE": {"intron": "down"},
    }[kind]
    directions = {}
    for key, direction in keys.items():
        for junction in event[key].split(";"):
            coord2int(junction)
            directions[junction] = direction
    if kind == "MSE":
        directions[event["intron"].split(";")[-1]] = "up"
    return directions


def event_region(event, extend_up, extend_down):
    if extend_up < 0 or extend_down < 0:
        raise ValueError("Event extensions must be nonnegative")
    directions = event_junctions(event)
    coords = [coord2int(junction) for junction in directions]
    chrom = coords[0][0]
    if any(not same_chrom(chrom, item[0]) for item in coords):
        raise ValueError("Event junctions span multiple chromosomes")
    # Shiba IDs contain the 1-based flanking exon bases, including in v1.0.
    # Extrema work for both old and new minus-strand MXE definitions.
    left = min(item[1] for item in coords) - 1
    right = max(item[2] for item in coords)
    upstream, downstream = (extend_down, extend_up) if event["strand"] == "-" else (extend_up, extend_down)
    return chrom, max(0, left - upstream), right + downstream


def event_boundaries(event):
    """RI exon–intron count IDs, distinct from spliced junction IDs.

    Shiba stores adjacent 1-based bases on either side of each boundary.
    Their left coordinate is the zero-based boundary position on either strand.
    """
    if event is None or not event["pos_id"].startswith("RI@"):
        return []
    chrom, left, right = coord2int(event["intron_a"])
    return [f"{chrom}:{left}-{left+1}", f"{chrom}:{right-1}-{right}"]


def posid2int(positional_id, shiba_path, extend_up, extend_down):
    event = read_event(positional_id, shiba_path)
    chrom, start, end = event_region(event, extend_up, extend_down)
    directions = event_junctions(event)
    return chrom, start, end, event["strand"], event.get("gene_name") or event.get("gene_id", ""), list(directions), directions
