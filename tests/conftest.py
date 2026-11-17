"""Small synthetic BAMs and Shiba output, generated without external tools."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import pysam
import pytest


SE_ID = "SE@chr1@201-250@150-301"


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_dataset(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    samples = ["sample1", "sample2", "sample3", "sample4"]
    for i, sample in enumerate(samples):
        path = root / f"{sample}.bam"
        with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": 2000}]}) as bam:
            for number in range(i+1):
                read = pysam.AlignedSegment()
                read.query_name = f"read{number}"
                read.query_sequence = "A" * 150
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 100
                read.mapping_quality = 60
                read.cigartuples = [(0, 50), (3, 50), (0, 50), (3, 50), (0, 50)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 150)
                bam.write(read)
        pysam.index(str(path))
    write_table(root / "experiment.tsv", ["sample", "bam", "group", "technology"], [
        {"sample": s, "bam": str(root/f"{s}.bam"), "group": "control" if i < 2 else "treatment", "technology": "short"}
        for i, s in enumerate(samples)])
    event = {"event_id": "SE_1", "pos_id": SE_ID, "exon": "chr1:201-250",
             "intron_a": "chr1:150-201", "intron_b": "chr1:250-301", "intron_c": "chr1:150-301",
             "strand": "+", "gene_id": "gene1", "gene_name": "Gene1", "label": "annotated"}
    write_table(root / "events/EVENT_SE.txt", list(event), [event])
    write_table(root / "results/splicing/PSI_matrix_sample.txt", ["event_id", "pos_id", *samples], [
        {"event_id": "SE_1", "pos_id": SE_ID, "sample1": "0.2", "sample2": "0.4", "sample3": "0.8", "sample4": ""}])
    rows = []
    for left, right, counts in [(150, 201, [1, 2, 3, 4]), (250, 301, [1, 2, 3, 4]), (150, 301, [4, 2, 1, 0]), (150, 151, [2, 4, 6, 8])]:
        rows.append(dict(zip(["chr", "start", "end", "ID", *samples], ["chr1", left, right, f"chr1:{left}-{right}", *counts])))
    write_table(root / "junctions/junctions.bed", ["chr", "start", "end", "ID", *samples], rows)
    annotation = root / "annotation/assembled_annotation.gtf"
    annotation.parent.mkdir(exist_ok=True)
    with annotation.open("w") as handle:
        for transcript, exons in [("tx_inclusion", [(101,150), (201,250), (301,350)]), ("tx_skipping", [(101,150), (301,350)])]:
            for start, end in exons:
                handle.write(f'1\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "gene1"; transcript_id "{transcript}";\n')
    return root


@pytest.fixture
def dataset(tmp_path):
    return make_dataset(tmp_path)
