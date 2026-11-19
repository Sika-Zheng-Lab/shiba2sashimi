"""Annotation tracks and event-local structures without inferred exon boundaries."""
import gzip
import logging
from collections import Counter, defaultdict
from pathlib import Path
import re

from .utils import coord2int, event_junctions, same_chrom

logger = logging.getLogger(__name__)


def exon_intervals(value):
    """Convert 1-based inclusive Shiba/GTF exons to half-open intervals."""
    intervals = []
    for item in value.split(";"):
        _, positions = item.rsplit(":", 1)
        start, end = map(int, positions.split("-"))
        if start < 1 or end < start:
            raise ValueError(f"Invalid exon coordinates: {item}")
        intervals.append((start - 1, end))
    return intervals


def event_features(event):
    if event is None:
        return []
    if event["pos_id"].startswith("RI@"):
        _, left, right = coord2int(event["intron_a"])
        return [(left, right - 1)]
    return sorted(set(interval for key, value in event.items()
                      if key == "exon" or key.startswith("exon_")
                      for interval in exon_intervals(value)))


def local_models(event):
    """Show only recorded exons and junctions; unknown flanking exons stay absent."""
    kind = event["pos_id"].split("@")[0]
    paths = {
        "SE": [("Inclusion", ["exon"], ["intron_a", "intron_b"]), ("Skipping", [], ["intron_c"])],
        "FIVE": [("Long exon", ["exon_a"], ["intron_a"]), ("Short exon", ["exon_b"], ["intron_b"])],
        "THREE": [("Long exon", ["exon_a"], ["intron_a"]), ("Short exon", ["exon_b"], ["intron_b"])],
        "MXE": [("Exon A", ["exon_a"], ["intron_a1", "intron_a2"]), ("Exon B", ["exon_b"], ["intron_b1", "intron_b2"])],
        "RI": [("Retention", ["exon_c"], []), ("Spliced", ["exon_a", "exon_b"], ["intron_a"])],
        "AFE": [("Distal", ["exon_a"], ["intron_a"]), ("Proximal", ["exon_b"], ["intron_b"])],
        "ALE": [("Distal", ["exon_a"], ["intron_a"]), ("Proximal", ["exon_b"], ["intron_b"])],
        "MSE": [("Inclusion", ["exon"], []), ("Skipping", [], [])],
    }[kind]
    result = []
    for i, (name, exon_keys, intron_keys) in enumerate(paths):
        exons = sorted(interval for key in exon_keys for interval in exon_intervals(event[key]))
        junctions = [j for key in intron_keys for j in event[key].split(";")]
        if kind == "MSE":
            junctions = event["intron"].split(";")[:-1] if i == 0 else event["intron"].split(";")[-1:]
        introns = [(coord2int(j)[1], coord2int(j)[2] - 1) for j in junctions]
        result.append({"name": f"{name} (local)", "exons": exons, "introns": introns, "strand": event["strand"]})
    return result


def read_gtf(path, chrom, start, end, gene_id=None):
    """Read exon structures from plain/gzipped GTF; keep transcript IDs intact."""
    opener = gzip.open if str(path).endswith(".gz") else open
    models = {}
    reference_counts = defaultdict(Counter)
    model_references = defaultdict(set)
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Invalid GTF row {number}: expected nine columns")
            if columns[2] != "exon" or not same_chrom(columns[0], chrom):
                continue
            attributes = dict(re.findall(r'(\S+)\s+"([^"]*)"', columns[8]))
            transcript = attributes.get("transcript_id")
            if not transcript:
                raise ValueError(f"GTF exon at line {number} is missing transcript_id")
            left, right = int(columns[3]) - 1, int(columns[4])
            if left < 0 or right <= left:
                raise ValueError(f"Invalid GTF exon coordinates at line {number}")
            key = (attributes.get("gene_id", ""), transcript)
            reference = attributes.get("ref_gene_id")
            if reference:
                reference_counts[key[0]][reference] += 1
                model_references[key].add(reference)
            model = models.setdefault(key, {"name": transcript, "strand": columns[6], "exons": []})
            model["exons"].append((left, right))
    result = []
    for key, model in models.items():
        if gene_id and key[0] != gene_id and gene_id not in model_references[key]:
            # StringTie novel transcripts lack ref_gene_id. Match Shiba's
            # most-common reference assignment within the assembled gene,
            # without replacing an explicit reference on another transcript.
            inferred = reference_counts[key[0]].most_common(1)
            if model_references[key] or not inferred or inferred[0][0] != gene_id:
                continue
        model["exons"] = sorted(set(model["exons"]))
        if model["exons"][0][0] >= end or model["exons"][-1][1] <= start:
            continue
        model["introns"] = [(a[1], b[0]) for a, b in zip(model["exons"], model["exons"][1:]) if a[1] < b[0]]
        result.append(model)
    return result


def load_models(shiba_path, chrom, start, end, event=None, gtf=None, max_transcripts=20):
    path = Path(gtf) if gtf else Path(shiba_path) / "annotation/assembled_annotation.gtf"
    if not gtf and not path.exists():
        path = Path(shiba_path) / "annotation/assembled.gtf"
    models = []
    if gtf or path.exists():
        models = read_gtf(path, chrom, start, end, event.get("gene_id") if event else None)
        if event:
            junctions = {(coord2int(j)[1], coord2int(j)[2] - 1) for j in event_junctions(event)}
            models.sort(key=lambda model: (-len(set(model["introns"]) & junctions), model["name"]))
        if len(models) > max_transcripts:
            logger.warning("Showing %d of %d transcripts; use --max-transcripts to show more", max_transcripts, len(models))
            models = models[:max_transcripts]
        if not models:
            logger.warning("No matching transcripts in %s for this region/gene", path)
    if not models and event:
        logger.info("Using event-local structures; unrecorded flanking exons are not inferred. Use --gtf for full transcript models.")
        models = local_models(event)
    return models
