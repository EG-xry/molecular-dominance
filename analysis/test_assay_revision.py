"""Regressions for empirical defects identified in the PLOS ONE/G3 review."""
import numpy as np
import pandas as pd
import pytest

import run_deconvolution as deconvolution
import run_kv7_rank as rank
import run_tp53 as tp53


@pytest.mark.parametrize("module,function", [
    (rank, rank.paired_measurements), (deconvolution, deconvolution.pair),
])
def test_synonymous_pairs_keep_nucleotide_identity_and_reject_duplicates(monkeypatch, module, function):
    """Repeated p.= must not generate a Cartesian product within an axis."""
    alone = pd.DataFrame({"hgvs_nt": ["c.3A>G", "c.6C>T"],
                          "hgvs_pro": ["p.=", "p.="], "value": [1.0, 2.0]})
    mixed = alone.assign(value=[0.5, 0.8])
    frames = {"alone": alone, "mixed": mixed}
    monkeypatch.setattr(module, "load", lambda tag: frames[tag].copy())
    paired = function("alone", "mixed", "value", None)
    assert len(paired) == 2
    assert paired.x.tolist() == [1.0, 2.0]
    assert paired.y.tolist() == [0.5, 0.8]
    frames["mixed"] = pd.concat([mixed, mixed.iloc[[0]]], ignore_index=True)
    with pytest.raises(pd.errors.MergeError):
        function("alone", "mixed", "value", None)


def test_filtered_control_cohorts_have_only_six_five_five_unique_variants():
    """Exercise deposited data and joins directly, independent of saved outputs."""
    control = {name: rank.ratios(*spec) for name, spec in rank.CONTROL.items()}
    filtered = {}
    for name, frame in control.items():
        assert len(frame) == frame.id.nunique() == 57
        x = rank.paired_measurements(*rank.CONTROL[name]).x
        rms = np.sqrt(np.mean(x ** 2))
        filtered[name] = frame[frame.x.abs() >= rms]
    counts = []
    for first, second in [("dV1/2", "+40 mV"), ("dV1/2", "-40 mV"),
                          ("+40 mV", "-40 mV")]:
        pair = filtered[first].merge(filtered[second], on="id", validate="one_to_one")
        counts.append(len(pair))
    assert counts == [6, 5, 5]


def test_ratio_denominator_zero_is_excluded_but_retained_for_control_rms(monkeypatch):
    source = pd.DataFrame({"hgvs_nt": ["c.3A>G", "c.6C>T"],
                           "hgvs_pro": ["p.=", "p.="], "value": [0.0, 2.0]})
    monkeypatch.setattr(rank, "load", lambda tag: source.copy())
    measured = rank.paired_measurements("alone", "mixed", "value", None)
    ratios = rank.ratios("alone", "mixed", "value", None)
    assert len(measured) == 2
    assert np.sqrt(np.mean(measured.x ** 2)) == pytest.approx(np.sqrt(2))
    assert ratios.id.tolist() == ["c.6C>T"]


def test_moment_correction_distinguishes_raw_second_moments_and_variances():
    """A shifted panel makes the three estimands demonstrably different."""
    xp, yp = np.array([2.0, 4.0]), np.array([3.0, 7.0])
    xc, yc = np.array([1.0, 2.0]), np.array([2.0, 1.0])
    assert deconvolution.correct(xp, yp, xc, yc, "as-published")[0] == pytest.approx(2.0)
    assert deconvolution.correct(xp, yp, xc, yc, "control-centred")[0] == pytest.approx(2.5)
    assert deconvolution.correct(xp, yp, xc, yc, "intercept")[0] == pytest.approx(3.0)
    for origin in ("control-centred", "intercept"):
        expected = deconvolution.correct(xp, yp, xc, yc, origin)[0]
        shifted = deconvolution.correct(xp + 10, yp - 3, xc + 10, yc - 3, origin)[0]
        assert shifted == pytest.approx(expected)


def test_moment_correction_does_not_report_nonpositive_signal_moment():
    result = deconvolution.correct([1, 2], [2, 3], [2, 3], [1, 2])
    assert np.isnan(result[0])


def test_tp53_common_cohort_contains_only_two_mixed_screens():
    """HCT116 must not be assigned the 7,446-variant mixed-screen cohort."""
    alone = tp53.load("mavedb_00000068-b-1.csv")
    mixed = {name: tp53.load(*spec) for name, spec in tp53.SCREENS.items()}
    monoallelic = {name: tp53.load(*spec) for name, spec in tp53.ALONE_STATE_SCREENS.items()}
    common, cohorts = tp53.cohorts(alone, mixed, monoallelic)
    shared = [(name, panel) for scope, name, panel in cohorts if scope == "common variant set"]
    assert len(common) == 7446
    assert len(shared) == 2
    for name, panel in shared:
        assert name in mixed
        assert set(panel.key) == common
        assert not panel.syn_x.any()
        assert panel.pos_x.nunique() == 392
    hct = [(scope, panel) for scope, name, panel in cohorts if "HCT116" in name]
    assert len(hct) == 1 and hct[0][0] == "own coverage"
    assert (~hct[0][1].syn_x).sum() == 3423


def test_tp53_parser_does_not_accept_complex_or_truncating_as_missense():
    assert tp53.parse("p.Arg175His") == ("R175H", 175, False)
    assert tp53.parse("NP_000537.3:p.Arg175His") == ("R175H", 175, False)
    assert tp53.parse("p.Arg175Arg") == ("R175=", 175, True)
    assert tp53.parse("p.Arg175=") == ("R175=", 175, True)
    for hgvs in ("p.Arg175Ter", "p.Ter394Tyr", "p.Arg175HisfsTer3", "p.Arg175His;Arg248Gln"):
        assert tp53.parse(hgvs) is None
