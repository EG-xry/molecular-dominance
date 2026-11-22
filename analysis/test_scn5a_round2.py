"""Independent numerical checks of both controls and paired-scale revision."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from scn5a_scenarios import control_floor, test_floor as evaluate_floor

HERE = Path(__file__).resolve().parent


def test_inert_control_uses_both_control_variances_and_reported_cells():
    d = pd.read_csv(HERE / "scn5a_peak.csv")
    r = d.set_index("variant").loc["W822X"]
    assert (r.het, r.het_se, r.het_n) == (134.2, 5.2, 164)
    scenario = control_floor(218.4, 7.7, r.het, r.het_se, 199, int(r.het_n))
    expected_variance = (2 * 134.2 / 218.4 * 5.2) ** 2 + (134.2 ** 2 / 218.4 ** 2 * 7.7) ** 2
    assert scenario["p"] == pytest.approx(.3855311355)
    assert scenario["floor"] == pytest.approx(82.4617216117)
    assert scenario["floor_se"] ** 2 == pytest.approx(expected_variance)
    propagated = evaluate_floor(d, scenario["floor"], scenario["floor_se"], scenario["floor_df"])
    fixed = evaluate_floor(d, scenario["floor"], 0.)
    assert list(propagated["counts"].values()) == [20, 14, 15, 14, 10]
    assert list(fixed["counts"].values()) == [20, 20, 22, 16, 11]
    assert set(fixed["intersection"]) - set(propagated["intersection"]) == {"A735E"}


def test_joint_surface_reaches_corrected_t_boundary_and_exposes_nine_df_failure():
    rows = pd.read_csv(HERE / "scn5a_joint_uncertainty.csv")
    data = pd.read_csv(HERE / "scn5a_peak.csv").set_index("variant")
    floor, floor_se = 100 ** 2 / 218.4, 100 ** 2 * 7.7 / 218.4 ** 2
    sample = rows[(rows.scenario == "SB-normalized linear-dose interpretation") &
                  (rows.variant.isin(["D356N", "G1712C"])) & (rows.effective_df == 9)]
    boundaries = {}
    for row in sample.itertuples():
        r = data.loc[row.variant]
        inflation = row.max_all_se_inflation_log_ratio
        z = np.log(floor / r.het) / np.sqrt((inflation * r.het_se / r.het) ** 2 +
                                         (inflation * floor_se / floor) ** 2)
        assert 2 * 51 * stats.t.sf(z, 9) == pytest.approx(.05, rel=1e-12)
        boundaries[row.variant] = inflation
    assert boundaries["D356N"] < 1 < boundaries["G1712C"]


def test_rebased_robustness_is_tied_to_each_floor():
    rows = pd.read_csv(HERE / "scn5a_robustness.csv")
    r = rows[(rows.scenario == "SB-normalized linear-dose interpretation") &
             (rows.variant == "D356N")].iloc[0]
    assert r.z_log_ratio == pytest.approx(4.5509216005)
    assert r.tau_star_propagated == pytest.approx(7.0008551934)
    assert r.tau_star_log_ratio == pytest.approx(8.2748423716)
    assert r.effective_df_star_log_ratio == pytest.approx(10.2981787917)
    critical = stats.norm.isf(.05 / (2 * 51))
    z = np.log(r.floor / r.mixed_mean) / np.sqrt((r.mixed_se / r.mixed_mean) ** 2 +
            (r.floor_se / r.floor) ** 2 + (r.tau_star_log_ratio / r.floor) ** 2)
    assert z == pytest.approx(critical)


def test_power_law_identity_and_exact_intersection_stability_interval():
    references = np.array([180., 200., 218.4, 250.])
    assert np.allclose(100 * .5 ** np.log2(references / 100), 10000 / references)
    report = json.loads((HERE / "scn5a_floor_scenarios.json").read_text())
    interval = report["same_two_floor_intervals"]["observed_WT_WT_relative_SE"]
    d = pd.read_csv(HERE / "scn5a_peak.csv")
    lo, hi = interval["lower_exclusive"], interval["upper_inclusive"]
    for floor in (lo + 1e-4, (lo + hi) / 2, hi - 1e-4):
        assert evaluate_floor(d, floor, floor * 7.7 / 218.4)["intersection"] == ["D356N", "G1712C"]
    assert evaluate_floor(d, lo - 1e-4, (lo - 1e-4) * 7.7 / 218.4)["intersection"] == ["G1712C"]
    assert "D785N" in evaluate_floor(d, hi + 1e-4, (hi + 1e-4) * 7.7 / 218.4)["intersection"]


def test_dose_factors_are_paired_with_references_and_contrast_survives():
    report = json.loads((HERE / "bounds.json").read_text())
    raw = pd.read_csv(HERE / "scn5a_peak.csv")
    panel = raw[raw.variant != "W822X"]
    for name in ("observed WT+WT reference with matched scale", "additive reference with matched scale"):
        row = report["paired_reference_dose_scenarios"][name]
        assert row["scale"] == row["reference"] / 100
        x, y = row["reference"] - row["scale"] * panel.homo, row["reference"] - panel.het
        slope = np.dot(x, y) / np.dot(x, x)
        assert row["slope"] == pytest.approx(slope)
        assert row["implied_g1"] == pytest.approx((slope - row["p"] ** 2) / (2 * row["p"] * (1-row["p"])))
        assert row["D785N_normal_bonferroni"] < .05 and row["D785N_welch_bonferroni"] < .05
    assert report["scn5a_scale"]["scale"] == 2.184
    actual = pd.read_csv(HERE / "bounds.csv").set_index("variant")
    for variant, expected in report["paired_reference_dose_scenarios"]["observed WT+WT reference with matched scale"]["table_J"].items():
        for column, value in expected.items():
            assert actual.loc[variant, column] == pytest.approx(value)


def test_uncertainty_display_uses_paired_scale_and_legacy_profile_is_labeled():
    uncertainty = json.loads((HERE / "uncertainty.json").read_text())["slopes"]
    value = uncertainty["SCN5A matched-reference scaling slope (s=2.184, W=218.4)"]
    pairs = json.loads((HERE / "bounds.json").read_text())["paired_reference_dose_scenarios"]
    assert value["b"] == pytest.approx(pairs["observed WT+WT reference with matched scale"]["slope"])
    assert value["scale_assumptions"]["combined_scale"] == 2.184
    legacy = json.loads((HERE / "profile_identity.json").read_text())["scn5a"]
    assert "historical unmatched-scale sensitivity" in legacy["analysis_scope"]
    assert legacy["combined_scale"] == 2 and legacy["reference"] == 218.4


def test_inert_fraction_changes_profile_without_changing_the_paired_slope():
    pairs = json.loads((HERE / "bounds.json").read_text())["paired_reference_dose_scenarios"]
    inert = pairs["inert p.Trp822Ter fraction with matched scale"]
    wt = pairs["observed WT+WT reference with matched scale"]
    p = 1 - 134.2 / 218.4
    assert inert["p"] == pytest.approx(p)
    assert inert["slope"] == wt["slope"]
    assert inert["slope_ci"] == wt["slope_ci"]
    expected = [(x - p ** 2) / (2 * p * (1-p)) for x in wt["slope_ci"]]
    assert inert["implied_g1_ci"] == pytest.approx(expected)
    assert inert["implied_g1_ci"][0] > 1
    assert inert["R104W"]["g1"] == pytest.approx(1.131522160469416)
    displayed = pd.read_csv(HERE / "bounds.csv").set_index("variant")
    for v, values in inert["table_J"].items():
        assert displayed.loc[v, "g1_at_inert_control_p"] == pytest.approx(values["g1_at_reference_p"])


def test_higher_floor_counts_follow_the_full_family_and_uncertainty_scenarios():
    report = json.loads((HERE / "scn5a_floor_scenarios.json").read_text())
    sensitivity = report["higher_floor_robustness"]
    d = pd.read_csv(HERE / "scn5a_peak.csv")
    floor, se = sensitivity["floor"], sensitivity["floor_se"]
    base_var = (d.het_se / d.het) ** 2 + (se / floor) ** 2
    z = np.log(floor / d.het) / np.sqrt(base_var)
    hit = sorted(d.loc[2 * len(d) * stats.t.sf(z, 9) < .05, "variant"])
    assert hit == ["D356N", "D785N", "G1712C", "V1405M"]
    assert sensitivity["effective_df_scenarios"][0]["log_ratio_variants"] == hit
    assert [(row["propagated_normal_count"], row["log_ratio_count"])
            for row in sensitivity["discrepancy_scenarios"]] == [(14, 10), (10, 8), (3, 7), (0, 3)]
    for row in sensitivity["discrepancy_scenarios"]:
        tau = row["additional_sd"]
        z = np.log(floor / d.het) / np.sqrt(base_var + (tau / floor) ** 2)
        assert sorted(d.loc[2 * len(d) * stats.norm.sf(z) < .05, "variant"]) == row["log_ratio_variants"]
