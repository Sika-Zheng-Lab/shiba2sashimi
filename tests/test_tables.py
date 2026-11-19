import pytest
from shiba2sashimi import tables, tracks
from conftest import SE_ID


def test_shiba_1_experiment_metadata_and_sample_prefix(dataset):
    result = tables.load_experiment_table(dataset / "experiment.tsv")
    assert list(result) == ["sample1", "sample2", "sample3", "sample4"]
    assert result["sample1"]["group"] == "control"
    assert result["sample1"]["technology"] == "short"


def test_reordered_columns_and_legacy_whitespace(tmp_path):
    path = tmp_path / "experiment.tsv"
    path.write_text("group sample bam\ncontrol sample1 a.bam\n")
    assert tables.load_experiment_table(path) == {"sample1": {"group": "control", "bam": "a.bam"}}


@pytest.mark.parametrize("text", ["sample\tbam\tgroup\n", "sample\tbam\tgroup\ns\ta\tg\ns\tb\tg\n",
                                   "sample\talignment\tbarcode\ns\ta\tb\n", "sample\tbam\tgroup\ns\ta\n"])
def test_invalid_experiment(tmp_path, text):
    path = tmp_path / "experiment.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        tables.load_experiment_table(path)


def test_psi_empty_last_cell(dataset):
    assert tables.get_psi_values(SE_ID, dataset) == {"sample1": 20., "sample2": 40., "sample3": 80., "sample4": "NA"}
    assert set(tables.get_psi_values("unknown", dataset).values()) == {"NA"}


def test_psi_nonfinite_and_scientific_notation(dataset):
    path = dataset / "results/splicing/PSI_matrix_sample.txt"
    path.write_text(f"event_id\tpos_id\ta\tb\tc\td\te\tf\nSE_1\t{SE_ID}\t1e-1\tnan\tinf\tNA\t-0.1\t1.1\n")
    assert tables.get_psi_values(SE_ID, dataset) == {"a": 10, "b": "NA", "c": "NA", "d": "NA", "e": "NA", "f": "NA"}


def test_selection_order_and_groups_override(dataset):
    experiment = tables.load_experiment_table(dataset / "experiment.tsv")
    assert tracks.select_samples(experiment, "sample4,sample1") == ["sample4", "sample1"]
    assert tracks.select_samples(experiment, "missing", "treatment,control") == ["sample3", "sample4", "sample1", "sample2"]


@pytest.mark.parametrize("samples,groups", [("missing", None), (None, "missing"), ("sample1,sample1", None), (None, "control,")])
def test_invalid_selection(dataset, samples, groups):
    with pytest.raises(ValueError):
        tracks.select_samples(tables.load_experiment_table(dataset / "experiment.tsv"), samples, groups)
