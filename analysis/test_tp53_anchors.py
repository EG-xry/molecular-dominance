"""Checks for the TP53 candidate-anchor diagnostics, not a DN calibration."""
import numpy as np
import pandas as pd
import pytest

import run_tp53 as tp53


def test_anchor_parser_retains_stops_without_admitting_stop_loss_or_complex_calls():
    assert tp53.parse_anchor("NP_000537.3:p.Arg175Ter") == ("R175*", 175, "nonsense")
    assert tp53.parse_anchor("p.Arg175Arg") == ("R175=", 175, "synonymous")
    assert tp53.parse_anchor("p.Arg175=") == ("R175=", 175, "synonymous")
    assert tp53.parse_anchor("p.Arg175His") == ("R175H", 175, "missense")
    for hgvs in ("p.Ter394Tyr", "p.Ter394=", "p.Arg175HisfsTer3",
                 "p.Arg175His;Arg248Gln", "p.?"):
        assert tp53.parse_anchor(hgvs) is None
    # The regression parser continues to exclude truncations.
    assert tp53.parse("p.Arg175Ter") is None


def test_anchor_loader_uses_finite_protein_means_and_does_not_overweight_nt_alleles(tmp_path, monkeypatch):
    source = pd.DataFrame({
        "hgvs_pro": ["p.Arg175=", "p.Arg175Arg", "p.Arg175Ter", "p.Arg175Ter",
                     "p.Arg175His", "p.Arg248Ter", "p.Arg273Ter", "p.Ter394Tyr"],
        "score": [0.0, 2.0, 4.0, 8.0, 3.0, np.nan, np.inf, 100.0]})
    source.to_csv(tmp_path / "scores.csv", index=False)
    monkeypatch.setattr(tp53, "SRC", str(tmp_path))
    frame = tp53.load_anchor_variants("scores.csv").set_index("key")
    assert set(frame.index) == {"R175=", "R175*", "R175H"}
    assert frame.loc["R175=", "score"] == 1.0
    assert frame.loc["R175*", "score"] == 6.0
    assert frame.loc["R175=", "n_input"] == 2
    assert frame.loc["R175*", "n_input"] == 2


def test_anchor_bootstrap_preserves_shared_residue_offsets():
    # Resampling the two classes independently would introduce spurious gap
    # variation despite each residue having exactly the same class contrast.
    frame = pd.DataFrame({"pos": [1, 2, 3, 1, 2, 3],
                          "variant_class": ["synonymous"] * 3 + ["nonsense"] * 3,
                          "score": [-4.0, 0.0, 7.0, -2.0, 2.0, 9.0]})
    lo, hi, valid = tp53.anchor_gap_ci(frame, np.random.default_rng(19), n_boot=200)
    assert lo == hi == pytest.approx(2.0)
    assert valid == 200


def test_deposited_mixed_controls_are_not_equivalent_to_alone_lof_anchors():
    """Exercise exact screen/column selection and both coverage definitions."""
    rows = tp53.assess_anchors(n_boot=200, seed=918)
    full = {r["screen"]: r for r in rows if r["scope"] == "all deposited positions"}
    assert len(rows) == 8 and len(full) == 4
    for row in rows:
        assert row["n_valid_boot"] == 200
        if row["state"] == "mixed":
            assert row["nonsense_minus_synonymous"] < 0
            assert row["gap_ci_hi"] < 0
            assert row["prob_nonsense_above_synonymous"] < 0.5
        else:
            assert row["nonsense_minus_synonymous"] > 0
            assert row["gap_ci_lo"] > 0
            assert row["prob_nonsense_above_synonymous"] > 0.5
    giacomelli = full["Giacomelli, WT background + nutlin-3"]
    boettcher = full["Boettcher, K562 WT-activity reporter"]
    assert giacomelli["n_synonymous"] == 377
    assert boettcher["n_synonymous"] == 36
    assert giacomelli["nonsense_minus_synonymous"] == pytest.approx(-0.31332522999137047)
    assert boettcher["nonsense_minus_synonymous"] == pytest.approx(-0.25120081237269576)
