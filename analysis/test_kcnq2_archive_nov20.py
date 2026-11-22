"""Regression tests for substantive KCNQ2 and archive review corrections."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import run_deconvolution as moment
import run_rank_null as ranks
import run_survey_validation as archive

HERE = Path(__file__).resolve().parent


def test_additive_offsets_break_raw_moment_correction_but_centring_recovers_signal():
    signal = np.repeat([1.0, 2.0, 4.0], 4)
    ex = np.tile([-1.0, -1.0, 1.0, 1.0], 3)
    ey = np.tile([-2.0, 2.0, -2.0, 2.0], 3)
    cx, cy, target = 3.0, -2.0, 0.4
    xp, yp, xc, yc = signal + cx + ex, target * signal + cy + ey, cx + ex, cy + ey
    raw = moment.correct(xp, yp, xc, yc, "as-published")[0]
    expected_biased = (target * np.mean(signal ** 2) + (cy + cx * target) * signal.mean()) / (np.mean(signal ** 2) + 2 * cx * signal.mean())
    assert raw == pytest.approx(expected_biased)
    assert raw != pytest.approx(target)
    assert moment.correct(xp, yp, xc, yc, "control-centred")[0] == pytest.approx(target)
    assert moment.correct(xp, yp, xc, yc, "intercept")[0] == pytest.approx(target)


def test_se_comparator_recovers_unattenuated_signal_for_balanced_error_design():
    signal = np.repeat([1.0, 2.0, 4.0], 2)
    error = np.tile([-1.0, 1.0], 3)
    assert moment.se_correct(signal + error, 0.4 * signal, np.ones(6)) == pytest.approx(0.4)
    assert np.isnan(moment.se_correct([1, 1], [2, 2], [2, 2]))


def test_boltzmann_example_respects_conductance_and_voltage_dependence():
    assert ranks.normalized_current(-40, 0, 0.6) == pytest.approx(0.6)
    assert ranks.normalized_current(40, 0, 0.6) == pytest.approx(0.6)
    assert ranks.normalized_current(-40, 10) == pytest.approx(0.4158, abs=0.0001)
    assert ranks.normalized_current(40, 10) == pytest.approx(0.9974, abs=0.0001)
    output = json.loads((HERE / "rank_null.json").read_text())
    assert output["status"] == "independent-endpoint null withdrawn"
    assert all("null_with_control_correlations" not in row for row in output["pairs"].values())


def test_every_broad_screen_record_is_classified_and_recall_is_not_invented(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = archive.validate(decisions=HERE / "survey_decisions.csv")
    assert result["n_labelled_experiment_sets"] == 46
    assert result["n_unclassified_flagged"] == 0
    assert result["n_verified_designs"] == 5
    assert result["recall_estimate"] is None
    assert result["prevalence_estimate"] is None
    assert result["n_unflagged_reviewed"] == 0
    labels = pd.read_csv(archive.LABELS)
    labels = labels[labels.experiment_prefix != "urn:mavedb:00000107"]
    path = tmp_path / "missing_broad_label.csv"
    labels.to_csv(path, index=False)
    with pytest.raises(ValueError, match="flagged experiment sets without a validation label"):
        archive.validate(decisions=HERE / "survey_decisions.csv", labels=path)
