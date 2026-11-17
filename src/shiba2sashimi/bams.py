"""Coverage from indexed BAM files."""
import numpy as np
import pysam
from .utils import same_chrom


def median_filter(data, window_size):
    if window_size < 1 or window_size % 2 == 0:
        raise ValueError("Smoothing window must be a positive odd integer")
    data = np.asarray(data, dtype=float)
    if window_size == 1:
        return data.copy()
    half = window_size // 2
    return np.array([np.median(data[max(0, i-half):i+half+1]) for i in range(len(data))])


def get_coverage(bam_path, chrom, start, end, window_size=1):
    """Raw A/C/G/T depth with pysam's default quality/read filters.

    Resolve exact chromosome names first, then an unambiguous chr alias. Let
    htslib locate BAI/CSI indexes (including sample.bai for sample.bam).
    """
    if start < 0 or end <= start:
        raise ValueError("Coverage region requires 0 <= start < end")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            raise FileNotFoundError(f"BAM index not found for {bam_path}; run samtools index")
        if chrom not in bam.references:
            candidates = [ref for ref in bam.references if same_chrom(ref, chrom)]
            if len(candidates) != 1:
                raise ValueError(f"Chromosome {chrom!r} not found unambiguously in {bam_path}")
            chrom = candidates[0]
        if end > bam.get_reference_length(chrom):
            raise ValueError(f"Region end {end} exceeds chromosome length for {chrom} in {bam_path}")
        counts = bam.count_coverage(chrom, start, end)
        coverage = np.sum(np.asarray(counts), axis=0, dtype=float)
    return median_filter(coverage, window_size)
