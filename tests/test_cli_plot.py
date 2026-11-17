import subprocess
import sys
import numpy as np
import pytest
import matplotlib.pyplot as plt
from shiba2sashimi import main, plots
from conftest import SE_ID


def arguments(dataset, output="plot.svg"):
    return ["-e", str(dataset / "experiment.tsv"), "-s", str(dataset), "--id", SE_ID,
            "-c", "chr1:100-350", "-o", str(dataset / output), "--font_family", "DejaVu Sans"]


@pytest.mark.parametrize("group_average", [False, True])
@pytest.mark.parametrize("extension", ["svg", "pdf", "png"])
def test_end_to_end_bam_to_plot(dataset, group_average, extension):
    args = arguments(dataset, f"plot.{extension}")
    if group_average:
        args += ["--group-average", "--groups", "treatment,control"]
    assert main.main(args) == 0
    assert (dataset / f"plot.{extension}").stat().st_size > 1000
    assert not plt.get_fignums()
    if extension == "svg":
        text = (dataset / "plot.svg").read_text()
        assert "tx_inclusion" in text
        assert "Target event:" not in text
        assert "mean PSI = 80.00% (n = 1)" in text if group_average else "sample1; PSI = 20.00%" in text


def test_fractional_junctions_and_aligned_transcripts(dataset):
    figure = main.run(main.parse_args(arguments(dataset) + ["--group-average"]))
    texts = [t.get_text() for ax in figure.axes for t in ax.texts]
    assert "1.5" in texts
    assert "control (n = 2); mean PSI = 30.00% (n = 2)" in texts
    assert all(ax.get_xlim() == (100., 350.) for ax in figure.axes)
    # The junction ends match GTF exon boundaries exactly.
    points = figure.axes[0].patches[0].get_path().vertices
    assert points[0, 0] == 150
    assert points[-1, 0] == 200
    rectangles = figure.axes[-2].patches
    assert any(p.get_x() == 200 and p.get_width() == 50 for p in rectangles)
    # Highlight the exon itself without a separate target-event header or bar.
    assert len(figure.axes[-2].texts) == 2
    assert all(line.get_color() != plots.EVENT_COLOR for line in figure.axes[-2].lines)


def test_empty_coverage_and_junctions(tmp_path):
    figure = plots.sashimi({"a": np.zeros(10)}, {"a": {}}, {"a": {"group": "g"}},
                           None, None, None, 5, "chr1", 0, 10, tmp_path / "empty.svg", font_family="DejaVu Sans")
    assert figure.axes[0].get_ylim()[1] > 0
    assert not plt.get_fignums()


def test_sample_selection_and_no_annotations(dataset):
    figure = main.run(main.parse_args(arguments(dataset) + ["--samples", "sample3", "--no-transcripts", "--nolabel", "--nojunc"]))
    assert len(figure.axes) == 2
    assert not figure.axes[0].texts
    assert not figure.axes[0].patches


def test_coordinate_only_does_not_require_psi_or_events(dataset):
    args = arguments(dataset)
    index = args.index("--id")
    del args[index:index+2]
    (dataset / "results/splicing/PSI_matrix_sample.txt").unlink()
    (dataset / "events/EVENT_SE.txt").unlink()
    assert main.main(args) == 0


def test_missing_sample_column_is_error_not_zero(dataset, caplog):
    path = dataset / "junctions/junctions.bed"
    path.write_text(path.read_text().replace("sample1", "other_sample"))
    assert main.main(arguments(dataset)) == 1
    assert "Samples missing from junction table: sample1" in caplog.text
    assert not (dataset / "plot.svg").exists()


def test_unrelated_coordinate_is_rejected(dataset, caplog):
    args = arguments(dataset)
    args[args.index("-c")+1] = "chr2:100-350"
    assert main.main(args) == 1
    assert "same chromosome" in caplog.text


@pytest.mark.parametrize("option,value", [("--width", "0"), ("--dpi", "0"), ("--extend_up", "-1"), ("--smoothing_window_size", "0"), ("--minimum_junc_reads", "nan"), ("--max-transcripts", "0")])
def test_cli_invalid_options(dataset, option, value):
    with pytest.raises(SystemExit) as error:
        main.parse_args(arguments(dataset) + [option, value])
    assert error.value.code == 2


def test_cli_requires_region(dataset):
    with pytest.raises(SystemExit):
        main.parse_args(["-e", "experiment", "-s", "shiba", "-o", "plot"])


def test_entry_point_exit_status(dataset):
    result = subprocess.run([sys.executable, "-m", "shiba2sashimi.main", *arguments(dataset), "--samples", "unknown"], capture_output=True, text=True)
    assert result.returncode == 1
    assert "not found" in result.stderr


def test_requested_font_is_not_silently_substituted():
    with pytest.raises(ValueError, match="not installed"):
        with plots.figure_style("DefinitelyMissingFont123"):
            pass


def test_group_smoothing_occurs_after_mean(dataset):
    # median(mean(a,b)) != mean(median(a),median(b)) for these depths.
    a, b = [0, 0, 6], [0, 6, 0]
    figure = plots.sashimi({"a": a, "b": b}, {}, {"a": {"group": "g"}, "b": {"group": "g"}},
                          None, None, None, 5, "chr1", 0, 3, dataset / "smooth.svg",
                          group_average=True, smoothing_window_size=3, nojunc=True, font_family="DejaVu Sans")
    area = figure.axes[0].collections[0].get_paths()[0]
    assert area.contains_point((1.5, 2.9))  # Mean then median gives depth 3, not 0, at the middle base.


def test_group_threshold_uses_mean_not_sum(dataset):
    figure = main.run(main.parse_args(arguments(dataset) + ["--group-average", "--minimum_junc_reads", "2"]))
    assert len(figure.axes[0].patches) == 1  # 1.5-read inclusion means are hidden, not summed to 3.
    assert len(figure.axes[1].patches) == 2
    assert figure.axes[0].get_ylim() == figure.axes[1].get_ylim()


def test_three_groups_with_unequal_sizes(dataset):
    path = dataset / "experiment.tsv"
    path.write_text(path.read_text().replace("sample4.bam\ttreatment", "sample4.bam\tthird"))
    figure = main.run(main.parse_args(arguments(dataset) + ["--group-average", "--groups", "third,control,treatment"]))
    labels = [ax.texts[0].get_text() for ax in figure.axes[:3]]
    assert labels[0].startswith("third (n = 1); mean PSI = NA (n = 0)")
    assert labels[1].startswith("control (n = 2)")
    assert labels[2].startswith("treatment (n = 1)")


def test_event_default_window_and_even_smoothing(dataset):
    args = arguments(dataset)
    index = args.index("-c")
    del args[index:index+2]
    figure = main.run(main.parse_args(args + ["--smoothing_window_size", "2"]))
    assert figure.axes[0].get_xlim() == (0, 801)


def test_no_junction_file_needed_when_disabled(dataset):
    (dataset / "junctions/junctions.bed").unlink()
    assert main.main(arguments(dataset) + ["--nojunc"]) == 0


def test_group_average_with_only_sample_subset(dataset):
    figure = main.run(main.parse_args(arguments(dataset) + ["--samples", "sample2,sample3", "--group-average"]))
    assert figure.axes[0].texts[0].get_text().startswith("control (n = 1); mean PSI = 40.00%")


def test_missing_psi_sample_column_is_error(dataset, caplog):
    path = dataset / "results/splicing/PSI_matrix_sample.txt"
    path.write_text(path.read_text().replace("sample1", "other_sample"))
    assert main.main(arguments(dataset)) == 1
    assert "Samples missing from PSI matrix" in caplog.text


def test_large_fractional_group_count_uses_fixed_point(dataset):
    figure = plots.sashimi(
        {"a": np.ones(250)*100000, "b": np.ones(250)*100001},
        {"a": {"chr1:150-201": 3001}, "b": {"chr1:150-201": 3020}},
        {"a": {"group": "g"}, "b": {"group": "g"}},
        None, None, None, 5, "chr1", 100, 350, dataset / "large.svg",
        group_average=True, font_family="DejaVu Sans")
    assert "3010.5" in [t.get_text() for t in figure.axes[0].texts]
    assert [t.get_text() for t in figure.axes[0].get_yticklabels()] == ["0", "100000.5"]


def test_close_ri_boundary_labels_do_not_overlap(dataset):
    figure = plots.sashimi(
        {"a": np.ones(1000)}, {"a": {"chr1:500-501": 12345, "chr1:510-511": 67890}},
        {"a": {"group": "g"}}, None, None, None, 5, "chr1", 0, 1000,
        dataset / "short_ri.svg", font_family="DejaVu Sans")
    figure.canvas.draw()
    labels = [t for t in figure.axes[0].texts if t.get_text().startswith("EI:")]
    assert len(labels) == 2
    assert not labels[0].get_window_extent().overlaps(labels[1].get_window_extent())


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("grouped", [False, True])
def test_ri_boundary_counts_from_cli(dataset, strand, grouped):
    from conftest import write_table
    event = {"pos_id": "RI@chr1@150-301", "strand": strand,
             "intron_a": "1:150-301", "exon_a": "1:101-150",
             "exon_b": "1:301-350", "exon_c": "1:101-350", "gene_id": "gene1"}
    write_table(dataset / "events/EVENT_RI.txt", list(event), [event])
    psi_path = dataset / "results/splicing/PSI_matrix_sample.txt"
    psi_path.write_text(psi_path.read_text().replace(SE_ID, event["pos_id"]))
    with (dataset / "junctions/junctions.bed").open("a") as handle:
        handle.write("chr1\t300\t301\tchr1:300-301\t5\t0\t3\t0\n")
        handle.write("chr1\t250\t251\tchr1:250-251\t99\t99\t99\t99\n")
    args = arguments(dataset)
    args[args.index("--id")+1] = event["pos_id"]
    args += ["--minimum_junc_reads", "0"]
    if grouped:
        args += ["--group-average"]
    figure = main.run(main.parse_args(args))
    texts = [t.get_text() for t in figure.axes[0].texts]
    assert ("EI: 3" if grouped else "EI: 2") in texts
    assert ("EI: 2.5" if grouped else "EI: 5") in texts
    assert "EI: 99" not in texts
    assert len(figure.axes[0].patches) == 1  # Only the spliced junction is an arc.
    assert [line.get_xdata()[0] for line in figure.axes[0].lines[::2]] == [150, 300]
    # Boundaries remain visible when transcript models are disabled.
    hidden_models = main.run(main.parse_args(args + ["--no-transcripts"]))
    assert any(t.get_text().startswith("EI:") for t in hidden_models.axes[0].texts)
    # The mean threshold also applies to boundary counts.
    threshold = main.run(main.parse_args(args + ["--minimum_junc_reads", "4"]))
    boundary_labels = [t.get_text() for t in threshold.axes[0].texts if t.get_text().startswith("EI:")]
    assert boundary_labels == ([] if grouped else ["EI: 5"])
    disabled = main.run(main.parse_args(args + ["--nojunc"]))
    assert not disabled.axes[0].lines and not disabled.axes[0].patches
