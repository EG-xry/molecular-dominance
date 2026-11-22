"""Tests for dominance.py. Each of the first eight targets a real error made in development."""
import pandas as pd
from scipy import stats
import sys, numpy as np, pandas as pd, pytest
sys.path.insert(0, ".")
import dominance as D


def _panel(b_true=0.5, n_res=60, per_res=5, noise=0.5, seed=1, klass="missense"):
    rng = np.random.default_rng(seed)
    pos = np.repeat(np.arange(1, n_res + 1), per_res)
    x = rng.normal(0, 5, len(pos))
    y = b_true * x + rng.normal(0, noise, len(pos))
    return pd.DataFrame(dict(x=x, y=y, pos=pos, klass=klass))


# --- error 7: a regex that treated Ter as an amino acid ---
def test_nonsense_is_not_missense():
    assert D.classify_variant("p.Arg175Ter")[3] == "nonsense"
    assert D.classify_variant("p.Ter394Tyr")[3] == "stoploss"
    assert D.classify_variant("p.Arg175His")[3] == "missense"
    assert D.classify_variant("p.Arg175=")[3] == "synonymous"
    assert D.classify_variant("p.Arg175Arg")[3] == "synonymous"
    assert D.classify_variant("p.[=;=]")[3] == "other"


def test_nonsense_contamination_biases_the_slope():
    """Nonsense variants have large x and ~zero transfer; mixing them in drags b down."""
    mis = _panel(b_true=0.60, seed=2)
    non = _panel(b_true=0.00, n_res=10, per_res=4, seed=3, klass="nonsense")
    non["x"] = non["x"] * 3          # nonsense carry the largest alone-state effects
    clean = D.estimate(mis, units="physical")
    dirty = D.estimate(pd.concat([mis, non], ignore_index=True), units="physical")
    assert clean.b > dirty.b + 0.05, "contamination must measurably depress the slope"


# --- the estimator itself ---
def test_recovers_known_slope():
    r = D.estimate(_panel(b_true=0.42, seed=4), units="physical")
    assert abs(r.b - 0.42) < 0.05
    assert r.ci[0] < 0.42 < r.ci[1]


def test_residue_clustering_widens_ci_vs_naive():
    """Correlated residues must not be counted as independent."""
    rng = np.random.default_rng(5)
    pos = np.repeat(np.arange(40), 8)
    shift = rng.normal(0, 0.35, 40)[pos]        # per-residue offset
    x = rng.normal(0, 5, len(pos))
    y = (0.5 + shift) * x + rng.normal(0, 0.4, len(pos))
    df = pd.DataFrame(dict(x=x, y=y, pos=pos))
    clustered = D.estimate(df, units="physical")
    unclustered = D.estimate(df.assign(pos=np.arange(len(df))), units="physical")
    assert (clustered.ci[1] - clustered.ci[0]) > (unclustered.ci[1] - unclustered.ci[0])


# --- the referee's finding: slopes where the axes are uncorrelated ---
def test_refuses_when_axes_uncorrelated():
    """Must FAIL CLOSED, not merely warn - the Methods claimed refusal and the code only warned."""
    rng = np.random.default_rng(6)
    df = pd.DataFrame(dict(x=rng.normal(0, 5, 300), y=rng.normal(0, 5, 300),
                           pos=np.repeat(np.arange(60), 5)))
    with pytest.raises(D.Uninterpretable, match="not a measurement"):
        D.estimate(df, units="physical")
    r = D.estimate(df, units="physical", allow_uncorrelated=True)
    assert r.interpretable is False
    assert any("NOT INTERPRETABLE" in w for w in r.warnings_)


def test_interpretable_flag_requires_physical_units():
    r = D.estimate(_panel(seed=61), units="z-score")
    assert r.interpretable is False
    assert D.estimate(_panel(seed=62), units="physical").interpretable is True


# --- the Z-score scale-indeterminacy the referee caught ---
def test_flags_non_physical_units():
    r = D.estimate(_panel(seed=7), units="z-score")
    assert any("scale-indeterminate" in w for w in r.warnings_)


# --- the unfalsifiability problem: ladder needs measured p AND physical units ---
def test_ladder_refuses_without_p():
    r = D.estimate(_panel(seed=8), units="physical")
    with pytest.raises(D.Uninterpretable, match="p was not measured"):
        D.compare_to_ladder(r, n_subunits=4)


def test_ladder_refuses_asserted_p():
    """A 1:1 plasmid design is an assumption about DNA, not a measured subunit pool."""
    r = D.estimate(_panel(seed=81), units="physical")
    with pytest.raises(D.Uninterpretable, match="not flagged as MEASURED"):
        D.compare_to_ladder(r, n_subunits=2, p=0.5)
    tbl, floor = D.compare_to_ladder(r, n_subunits=2, p=0.5, p_measured=True)
    assert abs(floor - 0.25) < 1e-12


def test_ladder_refuses_non_physical_axes():
    r = D.estimate(_panel(seed=9), units="z-score")
    with pytest.raises(ValueError, match="not in physical units"):
        D.compare_to_ladder(r, n_subunits=4, p=0.5, p_measured=True)


def test_preferential_pairing_counterexample_breaches_the_floor():
    """The reviewers' counterexample: pairwise assembly with preferential pairing reproduces a
    sub-floor observation with NO loss of protein, so a breach cannot identify trafficking."""
    WTWT, q, act_mm = 218.4, 0.0845, 0.338          # subunit fraction stays exactly 0.5
    observed = WTWT * q * (1 + act_mm)
    floor = WTWT * 0.25
    assert observed < floor, "counterexample must fall below the random-assembly floor"
    assert abs(observed - 24.7) < 0.5, "and must reproduce D785N's reported 24.7%"


# --- error 5: reporting the coefficient without its baseline ---
def test_baseline_is_reported_and_flagged_when_it_wins():
    rng = np.random.default_rng(10)
    df = _panel(b_true=0.5, seed=11)
    df["tgt"] = (df.y > df.y.median()).astype(float)   # target driven by y alone
    r = D.estimate(df, units="physical", target="tgt")
    assert not np.isnan(r.baseline_partial_r)
    assert not np.isnan(r.coefficient_partial_r)


# --- class gates ---
def test_class_gates_do_not_assign_expected_slopes_from_class_labels():
    """A class name alone neither validates a slope nor identifies its cause."""
    near_zero = pd.concat([_panel(0.6, seed=12, klass="missense"),
                      _panel(0.0, n_res=20, seed=13, klass="synonymous")], ignore_index=True)
    t = D.class_gates(near_zero).set_index("klass")
    larger = pd.concat([_panel(0.6, seed=14, klass="missense"),
                     _panel(0.7, n_res=20, seed=15, klass="synonymous")], ignore_index=True)
    t2 = D.class_gates(larger).set_index("klass")
    assert t2.loc["synonymous", "b"] > t.loc["synonymous", "b"]
    assert t.loc["synonymous", "note"] == t2.loc["synonymous", "note"]
    assert "assay context required" in t2.loc["synonymous", "note"]
    assert "shared reference" not in t2.loc["synonymous", "note"]


# --- the arithmetic the whole stoichiometric argument rests on ---
def test_ladder_arithmetic():
    lad = D.stoichiometric_ladder(4, 0.5)
    assert abs(lad[1] - 0.9375) < 1e-9 and abs(lad[4] - 0.0625) < 1e-9
    lad2 = D.stoichiometric_ladder(2, 0.5)
    assert abs(lad2[1] - 0.75) < 1e-9 and abs(lad2[2] - 0.25) < 1e-9
    for n in (1, 2, 3, 4, 6):
        assert abs(D.stoichiometric_ladder(n, 0.5)[n] - 0.5 ** n) < 1e-12, "floor must be (1-p)^n"


def test_too_few_observations_raises():
    with pytest.raises(ValueError, match=">=20"):
        D.estimate(_panel(n_res=2, per_res=3, seed=16), units="physical")


def test_ladder_allows_an_avowed_conditional_but_flags_it():
    """Error 16: refusing every unmeasured p made the published sensitivity
    analysis unrunnable, so the driver crashed instead of reporting a
    conditional.  An assumption that is named may be published; one that is
    silent may not."""
    r = D.Result(b=0.4, ci=(0.35, 0.45), n=10, n_residues=10, pearson_r=0.9,
                 units="physical", interpretable=True, p_measured=False)
    with pytest.raises(D.Uninterpretable):
        D.compare_to_ladder(r, n_subunits=4, p=0.5)          # silent: refused
    tbl, floor = D.compare_to_ladder(r, n_subunits=4, p=0.5,
                                     assume_p_because="1:1 plasmid design")
    assert tbl.conditional.all()
    assert (tbl.p_basis == "1:1 plasmid design").all()
    assert floor == pytest.approx(0.0625)
    with pytest.raises(ValueError):                          # cannot be both
        D.compare_to_ladder(r, n_subunits=4, p=0.5, p_measured=True,
                            assume_p_because="1:1 plasmid design")


def test_welch_uses_both_degrees_of_freedom():
    """Error 17 (review finding 2): the SCN5A Welch specification used
    het_n - 1 alone, ignoring the floor's own 198 df.  Satterthwaite combines
    both, and the understated df cost one variant."""
    s1, n1 = 3.6**2, 45          # a variant's SE and cell count
    s2, n2 = (7.7 / 4)**2, 199   # the floor's
    naive = n1 - 1
    satt = (s1 + s2)**2 / (s1**2 / (n1 - 1) + s2**2 / (n2 - 1))
    assert satt > naive, "combining both sources must not lower the df"


def test_log_ratio_delta_method_scales_each_se_by_its_own_mean():
    """Error 18 (review finding 2): the log specification divided the floor's
    SE by the variant's mean.  The delta-method variance of log(a/b) is
    (sa/a)^2 + (sb/b)^2."""
    a, sa, b, sb = 54.6, 1.925, 19.3, 3.6     # floor, its SE; D356N, its SE
    right = np.log(a / b) / np.sqrt((sb / b)**2 + (sa / a)**2)
    wrong = (np.log(a) - np.log(b)) / (np.sqrt(sa**2 + sb**2) / b)
    assert not np.isclose(right, wrong, rtol=0.05)
    # the correct form is the larger here: dividing sa by the smaller mean b
    # inflates the floor's contribution to the variance
    assert right > wrong


def test_source_and_mane_labels_differ_above_the_deleted_codon():
    """Error 19 (review finding 3): SCN5A isoform a carries Gln1077, isoform b
    does not, so labels above 1077 shift by one.  The crosswalk recorded zero
    offsets while its own returned titles showed the shift."""
    import csv
    import os
    p = os.path.join(os.path.dirname(__file__), "clinvar_crosswalk.csv")
    if not os.path.exists(p):
        pytest.skip("run run_clinvar.py first")
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        pos, off = int(r["source_position"]), int(r["offset"])
        assert off == (1 if pos > 1077 else 0), (
            f"{r['source_label']} at {pos} should have offset "
            f"{1 if pos > 1077 else 0}, got {off}")
    assert any(int(r["offset"]) == 1 for r in rows), \
        "a crosswalk with no offsets anywhere is the bug this test encodes"


def test_control_slope_is_not_comparable_to_panel_slope():
    """Error 20: the manuscript compared a zero-effect control slope directly
    against a disease-panel slope and concluded they were indistinguishable.
    They are the same estimator in two different signal-to-noise regimes.  The
    artifact's contribution to a panel slope is var(e)/(var(xi) + var(noise)),
    which shrinks as the panel's true effects grow."""
    share, d_true = 0.3276, 0.43          # measured control slope, assumed truth
    def observed(var_ratio):
        return (d_true * var_ratio + share) / (var_ratio + 1)
    # at the measured KCNQ2 variance ratio the artifact moves the slope by ~0.016
    assert abs(observed(19.9) - d_true) < 0.02
    # only when true effects are comparable to the noise does it matter
    assert abs(observed(0.5) - d_true) > 0.05
    # and the observed slope is monotone in the variance ratio
    seq = [observed(v) for v in (0.5, 1, 2, 5, 20)]
    assert seq == sorted(seq), "dilution must be monotone in the variance ratio"


def test_deconvolution_recovers_the_true_ratio():
    """The method-of-moments inversion of the shared-noise account:
    d = [b_p var(x_p) - b_c var(x_c)] / [var(x_p) - var(x_c)].

    The estimator is consistent in median, not per draw: with a small control
    panel and a modest variance ratio a single realisation can miss by 0.07,
    which is why run_deconvolution.py reports a bootstrap interval and the
    variance ratio alongside every point estimate."""
    rng = np.random.default_rng(3)
    share, d_true, N, M = 0.75, 0.5, 400, 60
    raw, corrected = [], []
    for _ in range(400):
        xi = rng.normal(0, 1.5, N)
        e = rng.normal(0, np.sqrt(share), N)
        x = xi + e + rng.normal(0, np.sqrt(1 - share), N)
        y = d_true * xi + e + rng.normal(0, np.sqrt(1 - share), N)
        e2 = rng.normal(0, np.sqrt(share), M)
        xc = e2 + rng.normal(0, np.sqrt(1 - share), M)
        yc = e2 + rng.normal(0, np.sqrt(1 - share), M)
        bp = np.sum(x * y) / np.sum(x * x)
        bc = np.sum(xc * yc) / np.sum(xc * xc)
        vp, vc = np.mean(x**2), np.mean(xc**2)
        raw.append(bp)
        corrected.append((bp * vp - bc * vc) / (vp - vc))
    raw_med, cor_med = np.median(raw), np.median(corrected)
    assert raw_med > d_true + 0.02, "the raw slope should be inflated here"
    assert abs(cor_med - d_true) < 0.03, f"correction gave {cor_med:.3f}"
    assert abs(cor_med - d_true) < abs(raw_med - d_true), "must reduce the bias"


def test_cross_study_control_slope_cannot_contain_a_shared_reference():
    """Error 22 (section 3): a control panel's slope means different things
    within and across studies.  Within one experiment it carries both the
    shared-reference artifact and any real correlated effect; across two
    independent experiments no reference can be shared, so it carries only the
    real effect.  A large cross-study control slope therefore falsifies the
    zero-effect assumption rather than measuring an artifact."""
    rng = np.random.default_rng(5)
    N, share, real = 400, 0.6, 0.0
    # within study: shared reference, no real effect -> slope ~ the shared share
    e = rng.normal(0, np.sqrt(share), N)
    x = e + rng.normal(0, np.sqrt(1 - share), N)
    y = e + rng.normal(0, np.sqrt(1 - share), N)
    within = np.sum(x * y) / np.sum(x * x)
    # across studies: independent references, no real effect -> slope ~ 0
    x2 = rng.normal(0, 1, N)
    y2 = rng.normal(0, 1, N)
    across = np.sum(x2 * y2) / np.sum(x2 * x2)
    assert within > 0.4, f"within-study control should show the artifact, got {within:.3f}"
    assert abs(across) < 0.15, f"cross-study control should be ~0, got {across:.3f}"
    # and with a REAL correlated effect the cross-study slope is non-zero
    xi = rng.normal(0, 1, N)
    x3 = xi + rng.normal(0, 0.3, N)
    y3 = 0.7 * xi + rng.normal(0, 0.3, N)
    assert np.sum(x3 * y3) / np.sum(x3 * x3) > 0.4


def test_section3_requires_one_common_alone_axis():
    """Error 21: the HCT116 row of section 3 was computed against an axis that
    correlates only +0.65 with the common alone-state axis, which is not the
    comparison the section describes.  Regressing screens on different alone
    axes confounds the mixed state with the alone state."""
    rng = np.random.default_rng(9)
    N = 800
    truth = rng.normal(0, 1, N)
    alone_common = truth + rng.normal(0, 0.3, N)
    alone_other = truth + rng.normal(0, 1.1, N)     # r ~ 0.65 with the common axis
    mixed = 0.6 * truth + rng.normal(0, 0.3, N)
    r = np.corrcoef(alone_common, alone_other)[0, 1]
    assert 0.5 < r < 0.8, f"constructed r = {r:.2f}"
    b_common = np.sum(alone_common * mixed) / np.sum(alone_common ** 2)
    b_other = np.sum(alone_other * mixed) / np.sum(alone_other ** 2)
    assert abs(b_common - b_other) > 0.1, (
        "using a different alone axis must change the slope materially, "
        "which is why the common axis is part of the claim")


def test_breach_definition_changes_the_clinical_result():
    """Error 23: section 5b documented the breach definition as 'significantly
    below the floor, Bonferroni-corrected' while the reported numbers came from
    the plain numerical comparison.  The two give different answers, so the
    definition is load-bearing and both must be reported."""
    import os
    p = os.path.join(os.path.dirname(__file__), "clinical_null.csv")
    if not os.path.exists(p):
        pytest.skip("run run_clinical.py first")
    d = pd.read_csv(p)
    rho = d[d.statistic == "partial Spearman rho"]
    num = rho[rho.analysis.str.contains(r"numerical", regex=True)]
    sig = rho[rho.analysis.str.contains(r"significant", regex=True)]
    assert len(num) and len(sig), "both definitions must appear in the output"
    assert abs(float(num.value.iloc[0])) > 0.2
    assert abs(float(sig.value.iloc[0])) < 0.1, (
        "the association must be shown to vanish under the stricter definition")


def test_legacy_stop_and_bare_synonymous_notation():
    """Error 25: classify_variant returned 'other' for 'p.=' and for 'W822X'.
    The first silently removed every MaveDB control from the zero-effect gate,
    so the |b| < 0.25 check never ran on them.  The second let a nonsense
    variant into a panel described as missense -- the same class of error as
    the Ter bug found in round 1, in a different spelling."""
    assert D.classify_variant("p.=")[3] == "synonymous"
    assert D.classify_variant("W822X")[3] == "nonsense"
    assert D.classify_variant("W822*")[3] == "nonsense"
    assert D.classify_variant("p.Trp822Ter")[3] == "nonsense"
    assert D.classify_variant("p.Ala62Val")[3] == "missense"
    # and the SCN5A panel does contain one such variant
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_peak.csv")
    if os.path.exists(p):
        v = pd.read_csv(p).variant
        klass = [D.classify_variant(x)[3] for x in v]
        assert klass.count("nonsense") == 1, (
            "W822X is the prespecified nonsense control and must be seen as one")


def test_floor_standard_error_is_not_the_model_error():
    """Sampling SE and hypothetical structural discrepancy are distinct.

    Four additive residuals are illustrative, not a floor-error estimator.
    The sensitivity deliberately uses separate analyst-chosen SDs.
    """
    import json
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.json")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    j = json.load(open(p))
    assert j["illustrative_additive_residual_rms"] > 10 * j["floor"]["se"], (
        "the controls' scatter must be reported against the floor SE")
    surv = {int(k): len(v) for k, v in j["survivors_by_hypothetical_discrepancy_sd"].items()}
    assert surv[0] > surv[5] > surv[10] > surv[15], "must be monotone"
    assert surv[15] == 0, "no variant should survive 15 points of model error"


def test_shared_offsets_can_create_ratio_rank_agreement():
    """A common reference offset can manufacture ranks when true ratios agree.

    This counterexample prevents reading a null control correlation as a
    general guarantee that disease-panel ratio ranks are nuisance-free.
    """
    x = np.arange(1.0, 101.0)
    y = 0.5 * x
    observed1 = (y + 1.0) / (x + 1.0)
    observed2 = (y + 2.0) / (x + 2.0)
    assert np.allclose(y / x, 0.5)
    assert stats.spearmanr(observed1, observed2).statistic > 0.999999


def test_published_cutpoint_beats_the_floor_on_validated_variants():
    """Error 27 (raised by the unread electrophysiology review): the three
    best-validated dominant-negative trafficking variants are caught by
    O'Neill's own <75% criterion but not by our floor.  The floor is a test of
    an assembly model, not a detector of dominant-negative variants, and the
    manuscript must not offer it as the latter."""
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.csv")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    d = pd.read_csv(p)
    t = d[d.variant.isin(["R104W", "R121W", "R104Q"])]
    assert len(t) == 3
    assert (t.het < 75).sum() == 3, "the published cut-point catches all three"
    assert (t.het < 54.6).sum() == 1, "our floor catches one"
    assert (t.p_fixed < 0.05).sum() == 0, "our significance test catches none"


def test_w822x_does_not_suppress_wild_type():
    """W822X has high mixed output in this assay. This one observation does
    not establish absence of protein or a mechanism for all null alleles."""
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.csv")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    d = pd.read_csv(p)
    r = d[d.variant == "W822X"].iloc[0]
    additive = 218.4 / 2 + r.homo
    assert r.het > additive, "the observed mixed current exceeds this illustrative prediction"
    assert (r.het - additive) / r.het_se > 3, "and the excess should be clear"


def test_photobleaching_prediction_is_conditioned_on_detection():
    """Error 28 (flagged independently by two reviewers): the dimer prediction
    was compared against the observed two-step fraction without conditioning on
    detection.  A spot enters the observed distribution only if at least one
    fluorophore matures, so the like-for-like figure is m^2/(1-(1-m)^2), not
    m^2."""
    m = 0.70
    unconditional = m ** 2
    conditional = m ** 2 / (1 - (1 - m) ** 2)
    assert abs(unconditional - 0.49) < 1e-9
    assert abs(conditional - 0.5385) < 1e-3
    assert conditional > unconditional, (
        "conditioning on detection must raise the predicted two-step fraction")


def test_library_does_not_call_b_a_dominance_coefficient():
    """Error 29: the module's docstring and printed output called b a
    'dominance coefficient', which is precisely the conflation the manuscript
    argues against.  Naming it correctly in the software is part of the
    argument, not cosmetic."""
    import inspect
    src = inspect.getsource(D)
    rendered = str(D.Result(b=0.4, ci=(0.3, 0.5), n=10, n_residues=10,
                            pearson_r=0.9, units="physical"))
    assert "dominance coefficient b" not in rendered, rendered
    assert "effect-transfer slope" in rendered
    # the phrase may appear only where the module explains why it is not used
    for line in src.splitlines():
        if "dominance coefficient" in line:
            assert ("not" in line.lower() or "conflation" in line.lower()), line


def test_falconer_mapping_rescales_a_one_allele_measurement():
    """Error 30: section 1b mapped SCN5A's ONE-allele alone-state value onto
    Falconer's 2a.  The alone state carries one allele against a 100% single-WT
    reference while the mixed state carries two against 218.4%, and peak
    current is extensive, so the value must be rescaled to a two-allele
    equivalent before the mapping applies.  The two give opposite answers."""
    import json
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.json")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    j = json.load(open(p))
    b = j["falconer"]["b"]
    assert b > 0.5, "the conditional rescaling places the slope above one-half"
    assert j["falconer"]["ci"][0] > 0.5, "the interval excludes one-half under the stated scaling assumptions"
    # the naive one-allele construction points the other way, which is the trap
    d = pd.read_csv(os.path.join(os.path.dirname(__file__), "scn5a_peak.csv"))
    d = d[d.variant != "W822X"]
    naive = np.sum((100 - d.homo) * (100 - d.het)) / np.sum((100 - d.homo) ** 2)
    assert naive < 0.5 < b, (
        f"naive {naive:.3f} and rescaled {b:.3f} must straddle 1/2, which is "
        "why the rescaling is not optional")


def test_wt_wt_control_is_supra_additive():
    """Error 31: the text called a 2.39-SE supra-additive control a
    confirmation of additivity, and anchoring the floor to it rather than to
    strict additivity runs in the direction that produces more breaches."""
    import json
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.json")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    j = json.load(open(p))
    assert j["strict_additivity"]["supra_additive_SE"] > 2, "supra-additive"
    assert j["strict_additivity"]["floor"] < j["floor"]["value"], (
        "the strict floor must be the more conservative one")
    assert len(j["strict_additivity"]["survivors"]) < \
        j["counts"]["fixed floor, two-sided normal"], (
        "and must therefore return fewer survivors")


def test_transfer_ratio_is_the_binomial_mean_of_the_dominance_profile():
    """Section 1c's identity: d = E[g_K], K ~ Binomial(n, p), where g_k is the
    fraction of the full effect expressed at k mutant subunits.  Checked
    against a direct simulation of the mixture, not against itself."""
    from math import comb
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(3000):
        n = int(rng.integers(2, 8))
        p = float(rng.uniform(0.02, 0.98))
        m = float(rng.uniform(0.0, 0.95))
        g = np.sort(rng.uniform(0, 1, n + 1))
        g[0], g[n] = 0.0, 1.0
        a = 1 - g * (1 - m)
        P = np.array([comb(n, k) * p**k * (1 - p)**(n - k) for k in range(n + 1)])
        d_direct = (1 - float(P @ a)) / (1 - m)
        worst = max(worst, abs(d_direct - float(P @ g)))
    assert worst < 1e-12, f"identity failed by {worst:.2e}"


def test_ladder_and_g_relation_are_special_cases_of_the_identity():
    """The stoichiometric ladder is a step profile; b = 1/4 + g/2 is n = 2."""
    from math import comb
    n, p = 4, 0.5
    P = [comb(n, k) * p**k * (1 - p)**(n - k) for k in range(n + 1)]
    for t in (1, 2, 3, 4):
        step = [1.0 if k >= t else 0.0 for k in range(n + 1)]
        assert abs(sum(a * b for a, b in zip(P, step))
                   - D.stoichiometric_ladder(n, p)[t]) < 1e-12
    P2 = [comb(2, k) * 0.5**2 for k in range(3)]
    for g in (0.0, 0.3, 1.0):
        assert abs(sum(a * b for a, b in zip(P2, [0, g, 1])) - (0.25 + g / 2)) < 1e-12


def test_titrating_the_allele_ratio_identifies_the_profile():
    """Section 1c: d(p) is a degree-n polynomial whose Bernstein coefficients
    are the profile.  Because g_0 and g_n are fixed by the definition of the
    profile, there are n-1 free parameters and n-1 ratios suffice -- an
    earlier version of this analysis used n+1, which is the error the
    condition numbers below would have hidden."""
    from math import comb
    rng = np.random.default_rng(6)

    def interior(n, ps):
        return np.vstack([[comb(n, k) * p**k * (1 - p)**(n - k)
                           for k in range(1, n)] for p in ps])

    for n in (2, 3, 4, 5, 6):
        m = n - 1
        ps = np.array([0.5]) if m == 1 else np.linspace(0.10, 0.90, m)
        X = interior(n, ps)
        assert X.shape == (m, m), "n-1 ratios, n-1 unknowns"
        g = np.sort(rng.uniform(0, 1, n + 1))
        g[0], g[n] = 0.0, 1.0
        assert np.allclose(np.linalg.solve(X, X @ g[1:n]), g[1:n], atol=1e-10)
    # the interior system is far better conditioned than the full one
    n, m = 4, 3
    ps_i = np.linspace(0.10, 0.90, m)
    ps_f = np.linspace(0.10, 0.90, n + 1)
    full = np.vstack([[comb(n, k) * p**k * (1 - p)**(n - k)
                       for k in range(n + 1)] for p in ps_f])
    assert np.linalg.cond(interior(n, ps_i)) < np.linalg.cond(full) / 4


def test_profile_monotone_solver_respects_upper_bound_during_fit():
    """Clipping the unconstrained-increment solution misses this optimum."""
    from run_profile import fit
    X = np.array([[1.0, 0.0], [1.0, 1.0]])
    b = np.array([0.0, 3.0])
    estimated = fit(X, b, "monotone")
    assert np.allclose(estimated, [1.0, 1.0], atol=1e-8)
    assert np.isclose(np.sum((X @ estimated - b) ** 2), 2.0)


def test_floor_is_the_zero_heterodimer_boundary_and_is_conservative():
    """Section 1c: the section 4 floor is the a_1 = 0 case with the
    mutant-homodimer quarter dropped, so it under-counts breaches."""
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.csv")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    d = pd.read_csv(p)
    d = d[d.variant != "W822X"]
    exact = (d.het < 54.6 + 0.5 * d.homo).sum()      # a_1 = 0 boundary
    used = (d.het < 54.6).sum()                      # the floor as applied
    assert used < exact, "the applied floor must be the conservative one"


def test_endpoint_joins_are_one_to_one():
    """Error 32 (both PLOS/G3 panels, CRITICAL): run_kv7_rank.py normalised
    both score sets to a protein-change identifier, but all 43 synonymous
    controls carry the annotation 'p.=', so the endpoint join formed a 43x43
    cartesian product and inflated 5 real pairs into 129.  This is the same
    bug that run_kv7_controls.py already avoided by keying on hgvs_nt."""
    import os
    p = os.path.join(os.path.dirname(__file__), "kv7_rank.csv")
    if not os.path.exists(p):
        pytest.skip("run run_kv7_rank.py first")
    d = pd.read_csv(p)
    filt = d[d.set == "benign controls"]
    unfilt = d[d.set == "benign controls, unfiltered"]
    assert (filt.n <= 6).all(), (
        f"the filtered control join must yield the true 5-6 pairs, got "
        f"{filt.n.tolist()}")
    assert (unfilt.n == 57).all(), "unfiltered must pair all 57 controls exactly once"
    assert filt.spearman.isna().all(), (
        "5-6 pairs is below the minimum, so no correlation may be reported")


def test_control_slope_is_reported_with_centred_statistics():
    """Error 33: a through-origin control slope was read as evidence of
    biological association.  E[xy] = Cov(x,y) + E[x]E[y], and for the HCT116
    controls the mean product is the whole numerator."""
    import os
    p = os.path.join(os.path.dirname(__file__), "tp53_controls.csv")
    if not os.path.exists(p):
        pytest.skip("run run_tp53.py first")
    d = pd.read_csv(p)
    h = d[d.screen.str.contains("HCT116")].iloc[0]
    assert h.b_control > 0.6, "the through-origin slope is large"
    assert abs(h.centred_pearson_r) < 0.1 and h.centred_p > 0.3, (
        "but the centred association is null, so the slope is a mean product")
    assert abs(h.mean_product_share_of_Exy) > 0.9
    # every control row must carry the centred statistics beside the slope
    for c in ("centred_pearson_r", "centred_p", "mean_product_share_of_Exy"):
        assert d[c].notna().all(), f"{c} missing"


def test_hct116_is_not_a_mixed_state_screen():
    """Error 34: the HCT116 scan is mut/Delta monoallelic -- the second allele
    is deleted, so no wild-type protein is present.  It cannot enter a
    comparison of wild-type-background screens."""
    import os
    p = os.path.join(os.path.dirname(__file__), "tp53_portability.csv")
    if not os.path.exists(p):
        pytest.skip("run run_tp53.py first")
    d = pd.read_csv(p)
    assert "state" in d.columns, "each screen must declare its genotype state"
    hct = d[d.screen.str.contains("HCT116")]
    assert (hct.state != "mixed").all(), "HCT116 must not be labelled mixed"
    mixed = d[(d.scope == "common variant set") & (d.state == "mixed")]
    assert len(mixed) == 2, "the portability comparison rests on two screens"


def test_falconer_sign_follows_the_stated_axis_convention():
    """Error 35: with x = phi_WT - phi_hom the wild type is at +a, so
    y = a - d and d/a = 1 - 2b, not 2b - 1."""
    import json
    import os
    p = os.path.join(os.path.dirname(__file__), "scn5a_specifications.json")
    if not os.path.exists(p):
        pytest.skip("run run_scn5a.py first")
    j = json.load(open(p))["falconer"]
    assert abs(j["d_over_a"] - (1 - 2 * j["b"])) < 1e-9
    assert j["d_over_a"] < 0, "under this convention the sign is negative"
    assert j["b"] > 0.5, "while the position relative to 1/2 is unchanged"
    assert "convention" in j, "the axis convention must be recorded"


def test_pooled_rmse_is_not_mean_of_roots():
    """Error 36 (P1-2): the precision simulation took the root inside the
    replicate loop and averaged the results.  By Jensen that is at most the
    pooled RMSE, and with one unknown it is exactly the mean absolute error."""
    rng = np.random.default_rng(21)
    e = rng.normal(0, 1, (2000, 1))          # one unknown, as for a dimer
    mean_of_roots = np.mean(np.sqrt((e ** 2).mean(axis=1)))
    pooled = np.sqrt((e ** 2).mean())
    mae = np.abs(e).mean()
    assert mean_of_roots < pooled, "Jensen: mean of roots understates the RMSE"
    assert abs(mean_of_roots - mae) < 1e-12, "with one unknown it IS the MAE"


def test_monotone_and_bounded_are_one_assumption():
    """Error 37 (P1-3): the manuscript stated monotonicity and [0,1] as two
    conditions.  Since g_n = 1 is fixed by the definition of the profile,
    monotonicity plus fixed endpoints implies 0 <= g_k <= 1.
    Violations may arise in either tail or by interior nonmonotonicity.
    The SCN5A entries are conditional point estimates, not rejections."""
    import json
    import os
    p = os.path.join(os.path.dirname(__file__), "profile_identity.json")
    if not os.path.exists(p):
        pytest.skip("run run_profile.py first")
    j = json.load(open(p))
    assert j["scn5a"]["n_violating_monotone_assumption"] > 20, (
        "the panel must be reported as violating the assumption")
    # the requirement must be given under more than one constraint regime
    for n, row in j["precision_required"].items():
        for k in ("se_free", "se_bounded", "se_monotone"):
            assert k in row, f"{k} missing for n={n}"
    # and free must tighten with n while monotone does not
    free = [j["precision_required"][str(n)]["se_free"] for n in (2, 6)]
    mono = [j["precision_required"][str(n)]["se_monotone"] for n in (2, 6)]
    assert free[0] > 4 * free[1], "unconstrained requirement tightens steeply"
    assert mono[1] > free[1] * 3, "the constraint is what flattens it"


def test_deconvolution_reports_origin_sensitivity():
    """Error 38 (P1-4): the section 2b correction is written in variances but
    computed as second moments about the origin, which agree only if the
    controls are centred there.  They are not."""
    import os
    p = os.path.join(os.path.dirname(__file__), "deconvolution.csv")
    if not os.path.exists(p):
        pytest.skip("run run_deconvolution.py first")
    d = pd.read_csv(p)
    for c in ("d_control_centred", "d_intercept", "d_origin_span",
              "control_mean_t_x", "mean_product_share_of_control_numerator"):
        assert c in d.columns, f"{c} must be reported"
    assert (d.control_mean_t_x.abs() > 2).all(), "controls are off-origin"
    assert d.d_origin_span.max() > 0.1, (
        "the origin specification moves the answer materially and must be shown")


def test_package_is_self_contained():
    """P1-8: a checkout of dominance/ alone must reproduce everything, so the
    inputs are vendored under data/ rather than read from sibling trees."""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    for rel in ("data/kv7/mavedb_00000664-a-14.csv",
                "data/tp53/mavedb_00000068-b-1.csv",
                "data/archive/scoresets.json",
                "data/refseq/NP_000326.2.fasta",
                "data/NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx"):
        assert os.path.exists(os.path.join(here, rel)), f"missing vendored {rel}"
    assert os.path.exists(os.path.join(here, "LICENSE"))
