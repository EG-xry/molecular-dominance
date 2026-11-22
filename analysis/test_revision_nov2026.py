"""Checks for the analyses added in the revision of 8 November 2026.

Where possible each check recomputes a reported quantity independently from
the input files rather than re-reading the script's own output.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

HERE = Path(__file__).resolve().parent
ALPHA, FAMILY = 0.05, 51


# ---------------------------------------------------------------- SCN5A ----
@pytest.fixture(scope="module")
def scn5a():
    peak = pd.read_csv(HERE / "scn5a_peak.csv")
    ctl = pd.read_csv(HERE / "scn5a_controls.csv").set_index("variant")
    spec = json.loads((HERE / "scn5a_specifications.json").read_text())
    return peak, ctl, spec


def _log_ratio_set(peak, ww, ww_se, p):
    floor, fse = ww * (1 - p) ** 2, ww_se * (1 - p) ** 2
    z = np.log(floor / peak.het) / np.sqrt((peak.het_se / peak.het) ** 2 + (fse / floor) ** 2)
    padj = np.minimum(1, FAMILY * 2 * stats.norm.sf(z))
    return sorted(peak.loc[padj < ALPHA, "variant"])


def test_fraction_grid_intersection_is_the_log_ratio_set(scn5a):
    peak, ctl, spec = scn5a
    ww, ww_se = ctl.loc["WT+WT", "het"], ctl.loc["WT+WT", "het_se"]
    expected = {"0.40": 10, "0.50": 4, "0.52": 2, "0.60": 1}
    for p, n in expected.items():
        s = _log_ratio_set(peak, ww, ww_se, float(p))
        assert len(s) == n
        assert spec["fraction_sensitivity"]["sets"][p]["intersection"] == s
    for v in spec["fraction_sensitivity"]["sets"].values():
        assert v["intersection_equals_log_ratio"]
    assert all(spec["nesting_at_p_0_5"].values())


def test_primary_counts_unchanged_at_half(scn5a):
    _, _, spec = scn5a
    row = [r for r in spec["fraction_sensitivity"]["table"] if abs(r["p"] - 0.5) < 1e-12][0]
    assert [row[k] for k in (
        "fixed floor, two-sided normal",
        "floor uncertainty propagated, two-sided normal",
        "floor uncertainty propagated, one-sided normal",
        "Welch-Satterthwaite, reported cell counts",
        "log-ratio delta method, two-sided normal")] == [9, 8, 9, 7, 4]
    assert spec["intersection"] == ["D356N", "D785N", "G1712C", "V1405M"]


def test_log_ratio_significance_fraction_has_closed_form(scn5a):
    peak, ctl, spec = scn5a
    ww, ww_se = ctl.loc["WT+WT", "het"], ctl.loc["WT+WT", "het_se"]
    zc = stats.norm.isf(ALPHA / (2 * FAMILY))
    for v in ("D356N", "D785N", "G1712C", "V1405M"):
        r = peak.set_index("variant").loc[v]
        sigma = np.sqrt((r.het_se / r.het) ** 2 + (ww_se / ww) ** 2)
        # log(floor/het) = log(ww/het) + 2 log(1-p) = zc * sigma at the tipping point
        p_star = 1 - np.exp((zc * sigma - np.log(ww / r.het)) / 2)
        got = spec["tipping_points"][v]["p_significance_rescue"]["log-ratio delta method, two-sided normal"]
        assert got == pytest.approx(p_star, abs=1e-7)


def test_discrepancy_tipping_point_reaches_critical_value(scn5a):
    peak, ctl, spec = scn5a
    ww, ww_se = ctl.loc["WT+WT", "het"], ctl.loc["WT+WT", "het_se"]
    floor, fse = ww / 4, ww_se / 4
    zc = stats.norm.isf(ALPHA / (2 * FAMILY))
    for v in ("D356N", "D785N", "G1712C", "V1405M"):
        r = peak.set_index("variant").loc[v]
        tau = spec["tipping_points"][v]["tau_star_propagated_two_sided"]
        z = (floor - r.het) / np.sqrt(r.het_se ** 2 + fse ** 2 + tau ** 2)
        assert z == pytest.approx(zc, rel=1e-9)
    # Consistent with the shipped discrepancy scenarios.
    assert spec["survivors_by_hypothetical_discrepancy_sd"]["5"] == ["D356N", "D785N", "G1712C"]
    assert spec["survivors_by_hypothetical_discrepancy_sd"]["10"] == ["G1712C"]


# ----------------------------------------------------------------- TP53 ----
def test_hotspot_anchoring_reproduces_and_is_bracketed():
    import run_tp53_hotspot as H
    from run_tp53 import ALONE_STATE_SCREENS, SCREENS, cohorts, load
    out = json.loads((HERE / "tp53_hotspot.json").read_text())
    alone = load(H.ALONE_FILE)
    mixed = {n: load(f, c) for n, (f, c) in SCREENS.items()}
    common, _ = cohorts(alone, mixed, {n: load(f, c) for n, (f, c) in ALONE_STATE_SCREENS.items()})
    assert len(common) == out["common_missense_n"] == 7446
    slopes = H.anchored_slopes(alone, mixed, common, H.HOTSPOTS)
    for n, b in slopes.items():
        assert b == pytest.approx(out["anchored_slopes"][n], rel=1e-12)
        lo, hi = out["leave_one_out_range"][n]
        assert lo <= b <= hi
        clo, chi = out["anchored_bootstrap_ci"][n]
        assert clo < b < chi
    assert out["hotspot_set"] == ["R175H", "G245S", "R248Q", "R248W", "R249S", "R273C", "R273H", "R282W"]
    assert out["raw_ratio_second_over_first"] == pytest.approx(0.211 / 0.619, abs=0.01)


def test_anchoring_is_affine_invariant():
    """Rescaling and shifting a deposited axis must not change anchored slopes."""
    import run_tp53_hotspot as H
    from run_tp53 import ALONE_STATE_SCREENS, SCREENS, cohorts, load
    alone = load(H.ALONE_FILE)
    mixed = {n: load(f, c) for n, (f, c) in SCREENS.items()}
    common, _ = cohorts(alone, mixed, {n: load(f, c) for n, (f, c) in ALONE_STATE_SCREENS.items()})
    base = H.anchored_slopes(alone, mixed, common, H.HOTSPOTS)
    moved = {n: f.assign(score=3.7 * f.score - 1.2) for n, f in mixed.items()}
    alone2 = alone.assign(score=0.4 * alone.score + 2.0)
    again = H.anchored_slopes(alone2, moved, common, H.HOTSPOTS)
    for n in base:
        assert again[n] == pytest.approx(base[n], rel=1e-10)


# -------------------------------------------------------------- archive ----
def test_archive_validation_counts(tmp_path, monkeypatch):
    import run_survey_validation as V
    decisions = HERE / "survey_decisions.csv"
    monkeypatch.chdir(tmp_path)
    r = V.validate(decisions=decisions)
    assert (r["strict_screen"]["n_true"], r["strict_screen"]["experiment_sets"]) == (3, 25)
    assert (r["ratio_screen"]["n_stating_mutant_wt_ratio"], r["ratio_screen"]["experiment_sets"]) == (0, 13)
    assert r["n_verified_designs"] == 5
    assert r["detection_by_strict_text"] == [3, 5] and r["detection_by_keyword"] == [2, 5]
    assert r["designs_stating_ratio_or_fraction"] == [0, 5]


def test_bad_evidence_quote_is_rejected(tmp_path, monkeypatch):
    import run_survey_validation as V
    labels = pd.read_csv(V.LABELS)
    labels.loc[labels.experiment_prefix == "urn:mavedb:00000094", "evidence_quote"] = "not in the record"
    bad = tmp_path / "labels.csv"
    labels.to_csv(bad, index=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="evidence quote"):
        V.validate(decisions=HERE / "survey_decisions.csv", labels=bad)


def test_zheng_case_roles_and_link():
    case = json.loads((HERE / "survey_zheng_case.json").read_text())["archive_evidence"]
    assert case["n_coexpression_score_sets"] == 8
    assert case["n_variants_coexpressed"] == 1187
    assert case["score_set_role_counts"]["trial upload"] == 1
    assert case["score_set_role_counts"]["KCNQ1 variants alone"] == 1
    assert case["experiment_primary_publication_identifiers"] == ["35760561"]


# ------------------------------------------------------------ selection ----
def test_selection_uses_the_full_kcnq2_current_panel():
    import run_selection as S
    x, y, sx, sy = S.effects("a-1", "a-4", "Normalized current density at +40 mV", 1.0)
    assert len(x) == 301
    assert S.slope(x, y) == pytest.approx(0.534, abs=5e-4)
    out = json.loads((HERE / "selection.json").read_text())
    e = out["endpoints"]["current at +40 mV"]
    assert e["n_all"] == 301 and e["n_selected"] < 301
    assert abs(e["simulation_relative_bias_selected"]) < 0.05


# ------------------------------------------------------------ Kv7 source ----
def test_kv7_extracts_rederive_exactly_from_vendored_workbooks():
    r = json.loads((HERE / "kv7_provenance.json").read_text())
    for panel in ("KCNQ4", "KCNQ2_KCNQ3"):
        ex = r["exhaustiveness"][panel]
        assert ex["missing_from_extract"] == [] and not ex["extract_rows_not_produced_by_rules"]
    assert r["inputs_unchanged_during_run"] is True
    ctx = r["context_point_slopes"]
    assert ctx["KCNQ4 C+E (rules K1-K4)"]["n"] == 190
    assert ctx["KCNQ4 C+E (rules K1-K4)"]["b"] == pytest.approx(0.3645, abs=1e-4)
    assert ctx["KCNQ4 C+D+E+F (K1 dropped)"]["n"] == 419
    assert ctx["KCNQ4 C+D+E+F (K1 dropped)"]["b"] == pytest.approx(0.2945, abs=1e-4)
    import hashlib
    expected = {"zheng2022_supplemental_dataset.xlsx": "bdfe8de2bedbeea7a4e64d152ef97ffb77dabcfded58a744517a2f1a105d7942",
                "vanoye2022_sdt3.xlsx": "bfee4246fa5214db0e375ecf340adf051ee774922168021a7b97f3e1a300d41a",
                "vanoye2022_sdt4.xlsx": "70569b45a9efb6d17698a9c453079b3e6808722e8c9c6e31dfa786d68dd989b6"}
    for name, digest in expected.items():
        assert hashlib.sha256((HERE / "data/kv7/source" / name).read_bytes()).hexdigest() == digest


# --------------------------------------------------------- HGVS labels ----
def test_scn5a_hgvs_matches_clinvar_crosswalk_and_reference():
    from hgvs import scn5a_hgvs
    c = pd.read_csv(HERE / "clinvar_crosswalk.csv")
    for src, mane in zip(c.source_label, c.mane_hgvs_p):
        assert scn5a_hgvs(src, accession=True) == mane
    assert scn5a_hgvs("W822X") == "p.Trp822Ter"
    assert scn5a_hgvs("G1712C") == "p.Gly1711Cys" and scn5a_hgvs("V1405M") == "p.Val1404Met"
    with pytest.raises(ValueError):
        scn5a_hgvs("A1712C")          # wrong reference residue must fail


# ----------------------------------------------------- rank-correlation null ----
def test_rank_null_observed_values_match_rank_analysis():
    null = json.loads((HERE / "rank_null.json").read_text())["pairs"]
    rank = json.loads((HERE / "kv7_rank.json").read_text())["rows"]
    obs = {r["endpoints"].replace(" vs ", "|"): r["spearman"] for r in rank if r["set"] == "disease panel"}
    for key, row in null.items():
        p, q = key.split(" vs ")
        assert row["observed_rho"] == pytest.approx(obs[f"{p}|{q}"], abs=5e-4)
        assert "null_with_control_correlations" not in row
        assert "descriptive only" in row["interpretation"]


# ------------------------------------------------- alternative intervals ----
def test_residue_intervals_reproduce_published_intervals():
    u = json.loads((HERE / "uncertainty.json").read_text())["slopes"]
    published = {"KCNQ4 activation-voltage panel": (0.28, 0.44),
                 "KCNQ2/KCNQ3 activation-voltage panel": (0.36, 0.56),
                 "KCNQ2 synonymous controls, current at +40 mV": (0.802, 1.031),
                 "SCN5A matched-reference scaling slope (s=2.184, W=218.4)": (0.711, 0.797),
                 "TP53 Giacomelli, WT background + nutlin-3": (0.563, 0.676)}
    for name, (lo, hi) in published.items():
        r = u[name]["residue_ci"]
        assert r[0] == pytest.approx(lo, abs=0.012) and r[1] == pytest.approx(hi, abs=0.012)
    assert "two_stage_ci" not in u["TP53 Boettcher, K562 WT-activity reporter"]


# ------------------------------------------------- published titrations ----
def test_titration_inversion_reproduces_and_identifies():
    import run_titration as T
    out = json.loads((HERE / "titration.json").read_text())["series"]
    hosy = out["Hosy & Vivaudou: SUR2A wild type + DNN"]["codings"]["stated_cDNA_ratios"]
    assert hosy["n_subunits"] == 4 and len(hosy["mutant_fractions"]) == 6
    assert hosy["one_profile_fit"]["dof"] == 3 and hosy["one_profile_fit"]["p_value"] > 0.05
    for j in range(1, 5):
        assert hosy["candidates"][f"requires >= {j} wild-type subunits"]["p_value"] < 0.05
    for k in ("requires >= 2 wild-type, adjacent pairs only (random placement)",
              "requires >= 2 wild-type, opposite pairs only (random placement)"):
        assert hosy["candidates"][k]["p_value"] > 0.05
    l325r = out["Clatot et al.: WT + L325R"]["codings"]
    assert l325r["stated_cDNA_ratios"]["profile_estimate"][0] == pytest.approx(0.796, abs=0.01)
    assert l325r["stated_cDNA_ratios"]["one_profile_fit"]["lack_of_fit_after_bonferroni"]
    assert not l325r["plotted_positions"]["one_profile_fit"]["lack_of_fit_after_bonferroni"]
    X = T.design([0.5], 4)
    assert X.tolist()[0] == pytest.approx([0.25, 0.375, 0.25])


def test_clinvar_records_are_the_queried_substitution():
    """A ClinVar text search also returns other alleles at the codon (and large
    deletions spanning it); each accepted record must name the queried change."""
    c = pd.read_csv(HERE / "clinvar_crosswalk.csv")
    found = c[c.clinvar_variation_id.notna()]
    assert len(found) == 47
    for title, mane in zip(found.clinvar_title, found.mane_hgvs_p):
        assert title.startswith("NM_000335.5(SCN5A)") and title.endswith(f"({mane.split(':')[-1]})"), title
