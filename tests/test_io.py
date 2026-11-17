import gzip
import numpy as np
import pytest
from shiba2sashimi import bams, junc, transcripts, utils
from conftest import SE_ID


def test_bam_depth_and_chr_alias(dataset):
    coverage = bams.get_coverage(str(dataset / "sample2.bam"), "chr1", 100, 350)
    assert coverage.shape == (250,)
    np.testing.assert_array_equal(coverage[:50], 2)
    np.testing.assert_array_equal(coverage[50:100], 0)
    np.testing.assert_array_equal(coverage[100:150], 2)
    np.testing.assert_array_equal(coverage[200:], 2)


def test_standard_basename_bai(dataset):
    (dataset / "sample1.bam.bai").rename(dataset / "sample1.bai")
    assert bams.get_coverage(str(dataset / "sample1.bam"), "1", 100, 101).tolist() == [1]


def test_csi_index(dataset):
    import pysam
    (dataset / "sample1.bam.bai").unlink()
    pysam.index("-c", str(dataset / "sample1.bam"))
    assert bams.get_coverage(str(dataset / "sample1.bam"), "1", 100, 101).tolist() == [1]


def test_missing_index_and_bad_chromosome(dataset):
    with pytest.raises(ValueError, match="Chromosome"):
        bams.get_coverage(str(dataset / "sample1.bam"), "chr9", 100, 200)
    with pytest.raises(ValueError, match="length"):
        bams.get_coverage(str(dataset / "sample1.bam"), "1", 100, 2001)
    (dataset / "sample1.bam.bai").unlink()
    with pytest.raises(FileNotFoundError, match="index"):
        bams.get_coverage(str(dataset / "sample1.bam"), "1", 100, 200)


def test_median_filter_preserves_fractional_means():
    np.testing.assert_array_equal(bams.median_filter([0., 1.5, 3.], 1), [0., 1.5, 3.])
    assert bams.median_filter([0, 1, 3], 3)[0] == .5
    with pytest.raises(ValueError):
        bams.median_filter([0, 1, 3], 0)


def test_tosa_counts_exclude_ri_boundaries(dataset):
    counts = junc.extract_junctions_in_region(dataset, "1", 100, 350)
    assert counts["sample1"] == {"chr1:150-201": 1, "chr1:250-301": 1, "chr1:150-301": 4}
    assert junc.extract_junctions_in_region(dataset, "1", 600, 700) == {s: {} for s in counts}


def test_junction_filter_retains_event_ids(dataset):
    counts = junc.extract_junctions_in_region(dataset, "chr1", 100, 350, ["1:150-201"])
    assert counts["sample2"] == {"1:150-201": 2}
    assert junc.extract_junctions_in_region(dataset, "chr1", 100, 350, [])["sample2"] == {}


def test_explicit_ri_boundaries_and_window_filter(dataset):
    ids = ["1:150-151", "1:300-301"]
    counts = junc.extract_junctions_in_region(dataset, "1", 100, 350, ids)
    assert counts["sample2"] == {"1:150-151": 4}  # Missing boundary is not fabricated.
    assert junc.extract_junctions_in_region(dataset, "1", 151, 350, ids)["sample2"] == {}


def test_gtf_conversion_and_event_highlight(dataset):
    event = utils.read_event(SE_ID, dataset)
    models = transcripts.load_models(dataset, "chr1", 100, 350, event)
    assert models[0]["name"] == "tx_inclusion"
    assert models[0]["exons"] == [(100, 150), (200, 250), (300, 350)]
    assert models[0]["introns"] == [(150, 200), (250, 300)]
    assert transcripts.event_features(event) == [(200, 250)]


def test_gz_gtf_and_missing_annotation_fallback(dataset):
    path = dataset / "annotation/assembled_annotation.gtf"
    zipped = dataset / "reference.gtf.gz"
    with gzip.open(zipped, "wt") as handle:
        handle.write(path.read_text())
    assert len(transcripts.read_gtf(zipped, "chr1", 100, 350)) == 2
    path.unlink()
    models = transcripts.load_models(dataset, "chr1", 100, 350, utils.read_event(SE_ID, dataset))
    assert models[0]["name"] == "Inclusion (local)"
    assert models[0]["exons"] == [(200, 250)]  # Never invent flanking exon extents.
    assert models[1]["exons"] == []
    with pytest.raises(FileNotFoundError):
        transcripts.load_models(dataset, "chr1", 100, 350, gtf="missing.gtf")


def test_gtf_gene_filter_and_limit(dataset, caplog):
    assert transcripts.read_gtf(dataset / "annotation/assembled_annotation.gtf", "chr1", 100, 350, "other_gene") == []
    result = transcripts.load_models(dataset, "chr1", 100, 350, utils.read_event(SE_ID, dataset), max_transcripts=1)
    assert len(result) == 1
    assert "Showing 1 of 2" in caplog.text


def test_header_only_junction_table(dataset):
    path = dataset / "junctions/junctions.bed"
    path.write_text(path.read_text().splitlines()[0] + "\n")
    result = junc.extract_junctions_in_region(dataset, "chr1", 0, 500)
    assert len(result) == 4
    assert all(value == {} for value in result.values())


def test_junctions_on_other_chromosomes_are_ignored(dataset):
    path = dataset / "junctions/junctions.bed"
    with path.open("a") as handle:
        handle.write("chr2\t160\t210\tchr2:160-210\t99\t99\t99\t99\n")
    assert all("chr2:160-210" not in junctions for junctions in junc.extract_junctions_in_region(dataset, "chr1", 0, 500).values())


def test_legacy_assembled_gtf_with_stringtie_reference_gene_ids(dataset):
    original = dataset / "annotation/assembled_annotation.gtf"
    legacy = dataset / "annotation/assembled.gtf"
    legacy.write_text(original.read_text().replace('gene_id "gene1";', 'gene_id "MSTRG.1"; ref_gene_id "gene1";'))
    original.unlink()
    models = transcripts.load_models(dataset, "chr1", 100, 350, utils.read_event(SE_ID, dataset))
    assert [m['name'] for m in models] == ['tx_inclusion', 'tx_skipping']
    assert models[0]['exons'] == [(100, 150), (200, 250), (300, 350)]
    assert transcripts.read_gtf(legacy, "chr1", 100, 350, "other_gene") == []


def test_stringtie_novel_transcript_inherits_gene_reference_independent_of_order(tmp_path):
    path = tmp_path / 'assembled.gtf'
    path.write_text(
        '1\ttest\texon\t201\t250\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "novel";\n'
        '1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "known"; ref_gene_id "gene1";\n'
        '1\ttest\texon\t301\t350\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "known"; ref_gene_id "gene1";\n'
        '1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "other"; ref_gene_id "gene2";\n')
    models = transcripts.read_gtf(path, 'chr1', 0, 400, 'gene1')
    assert [m['name'] for m in models] == ['novel', 'known']
    assert [m['name'] for m in transcripts.read_gtf(path, 'chr1', 0, 400, 'gene2')] == ['other']
