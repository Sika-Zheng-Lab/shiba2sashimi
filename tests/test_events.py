import pytest
from shiba2sashimi import utils, transcripts
from conftest import SE_ID, write_table


def event_row(kind, strand):
    row = {"event_id": f"{kind}_1", "pos_id": f"{kind}@chr1@test", "strand": strand, "gene_id": "gene1", "gene_name": "G"}
    if kind == "SE":
        row.update(exon="chr1:201-250", intron_a="chr1:150-201", intron_b="chr1:250-401", intron_c="chr1:150-401")
    elif kind in ("FIVE", "THREE"):
        if (kind == "FIVE" and strand == "+") or (kind == "THREE" and strand == "-"):
            row.update(exon_a="chr1:101-250", exon_b="chr1:101-200", intron_a="chr1:250-401", intron_b="chr1:200-401")
        else:
            row.update(exon_a="chr1:201-450", exon_b="chr1:251-450", intron_a="chr1:150-201", intron_b="chr1:150-251")
    elif kind == "MXE":
        row.update(exon_a="chr1:201-250", exon_b="chr1:301-350", intron_a1="chr1:150-201", intron_a2="chr1:250-401", intron_b1="chr1:150-301", intron_b2="chr1:350-401")
        if strand == "-":
            row.update(exon_a="chr1:301-350", exon_b="chr1:201-250", intron_a1="chr1:350-401", intron_a2="chr1:150-301", intron_b1="chr1:250-401", intron_b2="chr1:150-201")
    elif kind == "RI":
        row.update(exon_a="chr1:101-150", exon_b="chr1:401-450", exon_c="chr1:101-450", intron_a="chr1:150-401")
    elif kind in ("AFE", "ALE"):
        if (kind == "AFE" and strand == "+") or (kind == "ALE" and strand == "-"):
            row.update(exon_a="chr1:101-150;chr1:201-250", exon_b="chr1:301-350", intron_a="chr1:250-401;chr1:150-201", intron_b="chr1:350-401")
        else:
            row.update(exon_a="chr1:301-350;chr1:401-450", exon_b="chr1:201-250", intron_a="chr1:350-401;chr1:150-301", intron_b="chr1:150-201")
    else:
        row.update(exon="chr1:201-250;chr1:301-350", intron="chr1:150-201;chr1:250-301;chr1:350-401;chr1:150-401", mse_n="2")
    return row


@pytest.mark.parametrize("kind", utils.EVENT_TYPES)
@pytest.mark.parametrize("strand", ["+", "-"])
def test_all_event_types_and_strands(tmp_path, kind, strand):
    row = event_row(kind, strand)
    write_table(tmp_path / f"events/EVENT_{kind}.txt", list(row), [row])
    chrom, start, end, actual_strand, gene, junctions, directions = utils.posid2int(row["pos_id"], tmp_path, 20, 10)
    expected_left, expected_right = 149, 401
    if kind in ("FIVE", "THREE"):
        if (kind == "FIVE" and strand == "+") or (kind == "THREE" and strand == "-"):
            expected_left = 199
        else:
            expected_right = 251
    assert (chrom, start, end) == ("chr1", expected_left - (10 if strand == "-" else 20), expected_right + (20 if strand == "-" else 10))
    assert actual_strand == strand
    assert set(junctions) == set(directions)
    models = transcripts.local_models(row)
    assert len(models) == 2
    assert all(model["strand"] == strand for model in models)
    assert all(a < b for model in models for a, b in model["introns"] + model["exons"])


def test_minus_mxe_exon_a_is_upstream_in_transcript_orientation():
    event = event_row("MXE", "-")
    models = transcripts.local_models(event)
    assert models[0]["exons"] == [(300, 350)]
    assert set(models[0]["introns"]) == {(350, 400), (150, 300)}
    assert utils.event_junctions(event)["chr1:350-401"] == "down"


def test_legacy_minus_mxe_still_has_valid_window():
    event = event_row("MXE", "+")
    event["strand"] = "-"
    assert utils.event_region(event, 20, 10) == ("chr1", 139, 421)


def test_ri_highlight_is_retained_intron_only():
    assert transcripts.event_features(event_row("RI", "+")) == [(150, 400)]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_ri_boundaries_do_not_change_spliced_junctions(strand):
    event = event_row("RI", strand)
    assert utils.event_boundaries(event) == ["chr1:150-151", "chr1:400-401"]
    assert utils.event_junctions(event) == {"chr1:150-401": "up"}
    assert utils.event_boundaries(event_row("SE", strand)) == []
    assert utils.event_boundaries(None) == []


def test_missing_and_unknown_event(dataset):
    with pytest.raises(ValueError, match="not found"):
        utils.read_event("SE@chr1@unknown", dataset)
    with pytest.raises(ValueError, match="Unsupported"):
        utils.read_event("UNKNOWN@chr1", dataset)
    with pytest.raises(FileNotFoundError):
        utils.read_event("RI@chr1@100-200", dataset)


@pytest.mark.parametrize("coordinate", ["chr1:20-10", "chr1:-1-20", "chr1:10-10", "chr1", ":1-2"])
def test_invalid_coordinates(coordinate):
    with pytest.raises(ValueError):
        utils.coord2int(coordinate)


def test_explicit_chromosome_name_is_preserved():
    assert utils.coord2int("1:0-50") == ("1", 0, 50)
    assert utils.same_chrom("1", "chr1")
    assert utils.same_chrom("GL0001", "GL0001")


def test_single_base_exon_and_start_clipping():
    assert transcripts.exon_intervals("chr1:100-100") == [(99, 100)]
    assert utils.event_region(event_row("SE", "+"), 500, 0)[1] == 0
