import numpy as np
from shiba2sashimi.tracks import build_tracks


def test_group_means_include_zero_and_omit_missing_psi():
    experiment = {"a": {"group": "first"}, "b": {"group": "first"}, "c": {"group": "second"}, "d": {"group": "third"}}
    coverage = {"a": [2, 0], "b": [0, 4], "c": [3, 7], "d": [0, 0]}
    junctions = {"a": {"j": 3}, "b": {}, "c": {"j": 7}, "d": {}}
    psi = {"a": 20., "b": "NA", "c": 80., "d": float("nan")}
    result = build_tracks(coverage, junctions, experiment, ["d", "a", "b", "c"], psi, True)
    assert [t["label"] for t in result] == ["third", "first", "second"]
    assert result[0]["psi"] is None
    assert result[1]["n"] == 2
    np.testing.assert_array_equal(result[1]["coverage"], [1, 2])
    assert result[1]["junctions"]["j"] == 1.5
    assert result[1]["psi"] == 20
    assert result[1]["psi_n"] == 1
    assert result[2]["n"] == 1


def test_psi_is_sample_mean_not_read_weighted():
    result = build_tracks({"a": [10], "b": [90]}, {"a": {"j": 10}, "b": {"j": 90}},
                          {"a": {"group": "g"}, "b": {"group": "g"}}, ["a", "b"], {"a": 0., "b": 100.}, True)
    assert result[0]["psi"] == 50.
    assert result[0]["junctions"]["j"] == 50.


def test_sample_tracks_preserve_order_and_values():
    result = build_tracks({"a": [1], "b": [2]}, {}, {"a": {"group": "g"}, "b": {"group": "g"}}, ["b", "a"])
    assert [t["label"] for t in result] == ["b", "a"]
    assert result[0]["coverage"].tolist() == [2.]
    assert not result[0]["has_psi"]
