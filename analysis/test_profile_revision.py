"""Focused checks of the production profile estimator and its reporting."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import run_profile as profile


def test_production_rmse_pools_errors_before_taking_root(monkeypatch):
    seed, sigma, replicates = 92, 0.1, 20
    errors = np.random.default_rng(seed).normal(0, sigma, replicates) / 0.5
    monkeypatch.setattr(profile, "rng", np.random.default_rng(seed))
    actual = profile.rmse(2, [0.5], sigma, replicates,
                          mode="free", profile=np.array([0.3]))
    assert actual == pytest.approx(np.sqrt(np.mean(errors ** 2)))
    assert not np.isclose(actual, np.mean(np.abs(errors)))


def test_nonmonotone_fixed_activity_profile_has_exact_free_recovery():
    n, m = 3, 0.5
    g = np.array([0.0, 1.6, 0.2, 1.0])
    activities = 1 - (1 - m) * g
    ps = np.array([0.1, 0.9])
    measured = np.array([(1 - profile.bern(n, p) @ activities) / (1 - m)
                         for p in ps])
    X = profile.design(n, ps)
    b = measured - ps ** n
    assert np.all((activities >= 0) & (activities <= 1))
    assert np.allclose(profile.fit(X, b, "free"), g[1:-1])
    constrained = profile.fit(X, b, "monotone")
    assert np.sqrt(np.mean((constrained - g[1:-1]) ** 2)) > 0.4


def test_failed_profile_optimization_is_rejected(monkeypatch):
    monkeypatch.setattr(profile, "minimize", lambda *args, **kwargs:
                        SimpleNamespace(success=False, x=np.array([0.5]),
                                        message="iteration limit"))
    with pytest.raises(RuntimeError, match="iteration limit"):
        profile.fit(np.array([[0.5]]), np.array([0.2]), "monotone")


def test_reported_monotone_discrepancy_counts_include_both_tails():
    directory = Path(__file__).parent
    data = pd.read_csv(directory / "profile_scn5a.csv")
    report = json.loads((directory / "profile_identity.json").read_text())["scn5a"]
    below, above = data.g1_p50 < 0, data.g1_p50 > 1
    assert report["n_g1_below_zero"] == int(below.sum())
    assert report["n_g1_above_one"] == int(above.sum())
    assert report["n_violating_monotone_assumption"] == int((below | above).sum())
