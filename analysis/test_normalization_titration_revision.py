"""Independent checks for the source-normalization and titration corrections."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

HERE = Path(__file__).resolve().parent


def test_sb_reference_implies_revised_fraction_and_two_log_ratio_breaches():
    d = pd.read_csv(HERE / "scn5a_peak.csv")
    ctl = pd.read_csv(HERE / "scn5a_controls.csv").set_index("variant")
    sb, ww, se = ctl.loc["WT", "het"], ctl.loc["WT+WT", "het"], ctl.loc["WT+WT", "het_se"]
    fraction = (ww - sb) / ww
    floor, floor_se = sb ** 2 / ww, sb ** 2 * se / ww ** 2
    z = np.log(floor / d.het) / np.sqrt((d.het_se / d.het) ** 2 + (floor_se / floor) ** 2)
    survivors = sorted(d.loc[2 * len(d) * stats.norm.sf(z) < .05, "variant"])
    result = json.loads((HERE / "scn5a_specifications.json").read_text())["source_normalization"]
    assert fraction == pytest.approx(.5421245421)
    assert survivors == result["intersection_across_scenarios"] == ["D356N", "G1712C"]
    scenario = result["scenarios"]["SB-normalized linear-dose interpretation"]
    assert scenario["floor"] == pytest.approx(floor)
    assert scenario["floor_se"] == pytest.approx(floor_se)


def test_literal_ratio_sensitivity_preserves_documented_input_coordinates():
    from run_titration import coded_fractions
    raw = pd.DataFrame({"mutant_fraction_input": [0, .1, .25, .5, .75, .9, 1]})
    raw["wt_fraction_input"] = 1 - raw.mutant_fraction_input
    corrected = coded_fractions(raw, "Clatot", "stated_cDNA_ratios")
    assert corrected.mutant_fraction_input.tolist() == pytest.approx([0, 1/11, 1/5, 1/2, 4/5, 10/11, 1])
    assert raw.mutant_fraction_input.tolist() == [0, .1, .25, .5, .75, .9, 1]
    pd.testing.assert_frame_equal(coded_fractions(raw, "Clatot", "plotted_positions"), raw)


def test_candidate_multiplicity_and_failed_profile_fit_are_reported():
    from run_titration import CODINGS, holm_adjust
    out = json.loads((HERE / "titration.json").read_text())
    series = out["series"]
    for coding in CODINGS:
        fits = [s["codings"][coding] for s in series.values()]
        candidates = [c for fit in fits for c in fit["candidates"].values()]
        assert len(candidates) == out["candidate_multiplicity"]["family_size"] == 26
        for candidate in candidates:
            raw_p = stats.chi2.sf(candidate["delta_chi2"], candidate["dof"])
            assert candidate["p_bonferroni_26"] == pytest.approx(min(1, 26 * raw_p))
        assert [c["p_holm_26"] for c in candidates] == pytest.approx(holm_adjust([c["p_value"] for c in candidates]))
        for fit in fits:
            f = fit["one_profile_fit"]
            assert f["p_bonferroni_7"] == pytest.approx(min(1, 7 * f["p_value"]))
    literal, plotted = "stated_cDNA_ratios", "plotted_positions"
    for name in ("Clatot et al.: WT + L325R", "Clatot et al.: WT + L325R + difopein"):
        assert series[name]["codings"][literal]["one_profile_fit"]["lack_of_fit_after_bonferroni"]
        assert not series[name]["codings"][plotted]["one_profile_fit"]["lack_of_fit_after_bonferroni"]
    proportional = series["Hosy & Vivaudou: SUR2A wild type + DNN"]["codings"][literal]["candidates"]["proportional independent-subunits (g_k = k/n)"]
    assert proportional["p_bonferroni_26"] > .05 > proportional["p_holm_26"]


def test_fraction_and_endpoint_sensitivities_change_the_conclusion():
    out = json.loads((HERE / "titration.json").read_text())
    r104w = out["series"]["Clatot et al.: WT + R104W"]["sensitivity"]
    candidate = "heterodimer inactive (g1 = 1)"
    original = r104w["plotted_positions; fraction_shift=+0.00; mutant_endpoint_shift=0.00"]
    shifted = r104w["plotted_positions; fraction_shift=-0.01; mutant_endpoint_shift=0.00"]
    assert original["candidates"][candidate]["p_value"] < .05
    assert shifted["candidates"][candidate]["p_value"] > .05
    residual = r104w["plotted_positions; fraction_shift=+0.00; mutant_endpoint_shift=0.10"]
    assert residual["profile_estimate"] != original["profile_estimate"]
    assert out["sensitivity_grid"]["Clatot_mutant_endpoint_additions_relative_to_WT"][-1] == .10
