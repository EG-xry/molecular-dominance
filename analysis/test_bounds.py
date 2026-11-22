"""Checks of the partial-identification bounds in run_bounds.py.

The closed forms (B1)/(B2) are checked against the linear program they claim to
solve, the dimer special case against the general formula, and each asserted
property of g_1(p) -- the crossing at p*, where it is and is not monotone, and
the two p-free rejections -- against the function itself. The last tests read
bounds.csv/bounds.json, so run `python3 run_bounds.py` first.
"""
import json
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import linprog

import run_bounds as B

HERE = Path(__file__).resolve().parent
GRID_N = range(2, 7)
GRID_P = (0.05, 0.1, 0.25, 1 / 3, 0.4, 0.5, 0.55, 0.6, 0.75, 0.9)
GRID_D = (-0.2, 0.0, 0.03, 0.1, 0.25, 0.42, 0.5, 0.619, 0.75, 0.9, 0.999, 1.0, 1.05, 1.4)


def lp_both(n, p, d, k):
    """LP optimum of g_k over the same monotone set, or None if infeasible."""
    return B.lp_bound(n, p, d, k, "min"), B.lp_bound(n, p, d, k, "max")


# --- (1) the closed form is the linear program -----------------------------
def test_closed_form_equals_the_linear_program():
    worst, checked = 0.0, 0
    for n in GRID_N:
        for p in GRID_P:
            for d in GRID_D:
                gmin, gmax, feasible = B.bounds_at_p(n, p, d)
                if not feasible:
                    continue
                for k in range(1, n):
                    lo, hi = lp_both(n, p, d, k)
                    assert lo is not None and hi is not None, (n, p, d, k)
                    worst = max(worst, abs(gmin[k - 1] - lo), abs(gmax[k - 1] - hi))
                    checked += 1
    assert checked > 500
    assert worst < 1e-9, f"closed form departs from the LP by {worst:.2e}"


def witness(n, p, d, k, target):
    """The step profile that attains (B1) or (B2); the derivation's witness.

    max: everything below k is zero and everything from k up is the bound,
    unless the bound is capped at 1, when the slack moves below k instead.
    min: everything up to k is the bound and everything above it is 1, unless
    the bound is capped at 0, when the slack moves above k.
    """
    w = B.weights(n, p)
    g = np.zeros(n + 1)
    g[n] = 1.0
    if target == "max":
        raw = (d - w[n]) / w[k:n].sum()
        if raw <= 1.0:
            g[k:n] = raw
        else:
            g[k:n] = 1.0
            if k > 1:
                g[1:k] = (d - w[k:].sum()) / w[1:k].sum()
    else:
        raw = (d - w[k + 1:].sum()) / w[1:k + 1].sum()
        if raw >= 0.0:
            g[1:k + 1] = raw
            g[k + 1:n] = 1.0
        elif k + 1 < n:
            g[k + 1:n] = (d - w[n]) / w[k + 1:n].sum()
    return g


def test_closed_form_bounds_are_attained_by_a_feasible_profile():
    """Sharpness means the optimum is reached, not merely respected."""
    checked = 0
    for n, p, d in ((4, 0.5, 0.5), (5, 0.3, 0.42), (6, 0.6, 0.8), (3, 0.45, 0.25),
                    (4, 0.2, 0.5), (5, 0.7, 0.2), (3, 0.5, 0.8)):
        w = B.weights(n, p)
        gmin, gmax, feasible = B.bounds_at_p(n, p, d)
        assert feasible
        for k in range(1, n):
            for target, value in (("max", gmax[k - 1]), ("min", gmin[k - 1])):
                g = witness(n, p, d, k, target)
                assert np.all(np.diff(g) >= -1e-9), (n, p, d, k, target, g)
                assert np.all((g >= -1e-9) & (g <= 1 + 1e-9)), (n, p, d, k, target, g)
                assert float(w @ g) == pytest.approx(d, abs=1e-9)
                assert g[k] == pytest.approx(value, abs=1e-9)
                checked += 1
    assert checked > 40


def test_infeasible_cases_agree_with_the_lp_up_to_its_own_tolerance():
    """Where the closed form reports an empty set the LP must too, except in a
    band the width of HiGHS's feasibility tolerance around d = w_n."""
    disagreements = []
    for n in GRID_N:
        for p in GRID_P:
            for d in GRID_D:
                gmin, gmax, feasible = B.bounds_at_p(n, p, d)
                if feasible:
                    continue
                lo, hi = lp_both(n, p, d, 1)
                if lo is not None and hi is not None:
                    slack = max(B.weights(n, p)[n] - d, d - (1 - B.weights(n, p)[0]))
                    disagreements.append((n, p, d, slack))
    assert all(s < 1e-6 for *_, s in disagreements), disagreements


# --- (2) the dimer point case against the general formula -------------------
def test_dimer_point_equals_the_general_bounds_when_inside_the_unit_interval():
    for p in np.linspace(0.05, 0.95, 19):
        for d in np.linspace(-0.2, 1.2, 29):
            point = float(B.dimer_g1(p, d))
            gmin, gmax, feasible = B.bounds_at_p(2, p, d)
            if -1e-9 <= point <= 1 + 1e-9:
                assert feasible
                assert gmin[0] == pytest.approx(min(max(point, 0.0), 1.0), abs=1e-12)
                assert gmax[0] == pytest.approx(min(max(point, 0.0), 1.0), abs=1e-12)
            else:
                assert not feasible, (p, d, point)


def test_dimer_needs_no_monotonicity_but_n_above_two_does():
    """One interior unknown and one equation: g_1 is identified outright, which
    is why g_1 > 1 is a conclusion for a dimer and an assumption violation for
    larger n. At the same (p, d) a trimer leaves an interval half as wide as
    the whole parameter space."""
    p, d = 0.5, 0.50
    assert B.dimer_g1(p, d) == pytest.approx((0.50 - 0.25) / 0.5) == pytest.approx(0.5)
    assert B.bounds_at_p(2, p, d)[0][0] == pytest.approx(B.bounds_at_p(2, p, d)[1][0])
    gmin, gmax, _ = B.bounds_at_p(3, p, d)
    assert (gmin[0], gmax[0]) == pytest.approx((0.0, 0.5))


# --- (3) the allele-ratio tolerance ----------------------------------------
def test_p_star_is_the_crossing_of_g1_equals_one():
    assert float(B.allele_ratio_tolerance(0.0)) == 0.0   # p* leaves the open range
    for d in np.linspace(1e-3, 0.999, 60):
        ps = float(B.allele_ratio_tolerance(d))
        assert B.dimer_g1(ps, d) == pytest.approx(1.0, abs=1e-8)
        if ps > 0.02:
            assert B.dimer_g1(ps * 0.9, d) > 1.0          # below p*: a_1 < m required
        assert B.dimer_g1(min(0.999, ps + 0.01), d) < 1.0  # above it: not required
    for d in (-1e-9, -0.1, -0.5):                          # no admissible boundary
        assert np.isnan(B.allele_ratio_tolerance(d))
        assert np.all(B.dimer_g1(np.linspace(1e-3, 1 - 1e-3, 999), d) < 1.0)


def test_p_star_hand_worked():
    """d = 3/4 puts the crossing at p = 1/2: g_1(1/2) = (3/4 - 1/4)/(1/2) = 1."""
    assert float(B.allele_ratio_tolerance(0.75)) == pytest.approx(0.5)
    assert float(B.dimer_g1(0.5, 0.75)) == pytest.approx(1.0)
    assert float(B.allele_ratio_tolerance(0.0)) == pytest.approx(0.0)
    assert float(B.allele_ratio_tolerance(1.0)) == pytest.approx(1.0)
    assert np.isnan(B.allele_ratio_tolerance(1.3))             # model incompatible
    assert float(B.supra_wildtype_tolerance(0.25)) == pytest.approx(0.5)


def test_d_above_one_rejects_at_every_fraction_and_every_order():
    ps = np.linspace(1e-4, 1 - 1e-4, 20001)
    for d in (1.000001, 1.01, 1.0134, 1.1534, 1.3396, 2.0, 5.0):
        assert B.dimer_g1(ps, d).min() > 1.0
        turn_p, turn_g = B.dimer_turning_point(d)
        assert turn_g == pytest.approx(B.dimer_g1(turn_p, d), abs=1e-9)
        assert turn_g == pytest.approx(B.dimer_g1(ps, d).min(), abs=1e-4)
        for n in range(2, 9):
            for p in (0.05, 0.3, 0.5, 0.7, 0.95):
                assert not B.bounds_at_p(n, p, d)[2]
        assert B.feasible_p_interval(3, d) is None


def test_d_below_zero_rejects_in_the_other_direction():
    ps = np.linspace(1e-4, 1 - 1e-4, 20001)
    for d in (-1e-6, -0.1053, -0.5):
        assert B.dimer_g1(ps, d).max() < 0.0
        turn_p, turn_g = B.dimer_turning_point(d)
        assert turn_p < 0.5 and turn_g < 0.0
        assert turn_g == pytest.approx(B.dimer_g1(ps, d).max(), abs=1e-4)
        for n in range(2, 9):
            assert not B.bounds_at_p(n, 0.5, d)[2]


# --- (4) monotonicity in p, where it holds and where it does not ------------
def test_dimer_g1_is_strictly_decreasing_in_p_exactly_on_zero_to_one():
    ps = np.linspace(1e-3, 1 - 1e-3, 4000)
    for d in np.linspace(0.0, 1.0, 41):
        assert np.all(np.diff(B.dimer_g1(ps, d)) < 0), d
    for d in (-0.2, 1.2):                      # NOT monotone outside [0, 1]
        assert not np.all(np.diff(B.dimer_g1(ps, d)) < 0)
        assert not np.all(np.diff(B.dimer_g1(ps, d)) > 0)


def test_derivative_sign_matches_the_closed_form():
    h = 1e-7
    for d in np.linspace(-0.5, 2.0, 51):
        for p in np.linspace(0.05, 0.95, 19):
            numeric = (B.dimer_g1(p + h, d) - B.dimer_g1(p - h, d)) / (2 * h)
            if abs(numeric) > 1e-3:
                assert np.sign(numeric) == np.sign(d * (2 * p - 1) - p ** 2)


def test_each_bound_is_monotone_in_p_across_the_feasible_range():
    for n in (2, 3, 4, 5):
        for d in (0.2, 0.3645, 0.4227, 0.619, 0.8):
            lo_p, hi_p = B.feasible_p_interval(n, d)
            ps = np.linspace(lo_p + 1e-6, hi_p - 1e-6, 600)
            curves = np.array([np.concatenate(B.bounds_at_p(n, p, d)[:2]) for p in ps])
            for col in range(curves.shape[1]):
                diffs = np.diff(curves[:, col])
                assert np.all(diffs <= 1e-12), (n, d, col)


def test_hull_over_p_uses_the_endpoints_and_the_dimer_hull_adds_the_turning_point():
    lo, hi, share = B.bounds_over_p(4, 0.4227, 0.4, 0.6)
    assert share == 1.0
    assert hi[0] == pytest.approx(B.bounds_at_p(4, 0.4, 0.4227)[1][0], abs=1e-6)
    # a dimer whose interior minimum lies inside the range: endpoints are not enough
    d = 1.0 / (2 * 0.5 - 0.5 ** 2 / 0.5)       # any d > 1 with a turn near 0.5
    d = 2.0                                     # turning point 2 - sqrt(2) = 0.586
    turn_p, turn_g = B.dimer_turning_point(d)
    assert 0.4 < turn_p < 0.6
    lo2, hi2 = B.dimer_g1_hull(d, 0.4, 0.6)
    assert lo2 == pytest.approx(turn_g)
    assert lo2 < min(B.dimer_g1(0.4, d), B.dimer_g1(0.6, d))


# --- (5) hand-worked examples ----------------------------------------------
def test_hand_worked_tetramer_at_half():
    """n = 4, p = 1/2 gives w = (1,4,6,4,1)/16. With d = 1/2 the arithmetic is
        g_1 <= (1/2 - 1/16)/(14/16) = 1/2,   g_1 >= (1/2 - 11/16)/(4/16) < 0
        g_2 in [(1/2 - 5/16)/(10/16), (1/2 - 1/16)/(10/16)] = [3/10, 7/10]
        g_3 >= (1/2 - 1/16)/(14/16) = 1/2,   g_3 <= (7/16)/(4/16) > 1
    all of which are checked here against the exact fractions."""
    gmin, gmax, feasible = B.bounds_at_p(4, 0.5, 0.5)
    assert feasible
    assert gmax == pytest.approx([0.5, 0.7, 1.0])
    assert gmin == pytest.approx([0.0, 0.3, 0.5])
    assert B.weights(4, 0.5) == pytest.approx(np.array([1, 4, 6, 4, 1]) / 16)


def test_hand_worked_tetramer_panel_bound():
    """A tetramer panel slope of 0.4227 at p = 0.4: w_4 = 0.0256 and
    w_1 + w_2 + w_3 = 0.8448, so g_1 <= (0.4227 - 0.0256)/0.8448 = 0.4701."""
    _, gmax, _ = B.bounds_at_p(4, 0.4, 0.4227)
    assert gmax[0] == pytest.approx((0.4227 - 0.4 ** 4) / 0.8448, abs=1e-12)
    assert gmax[0] == pytest.approx(0.4701, abs=5e-4)


def test_feasible_p_interval_matches_the_direct_condition():
    for n in (2, 3, 4, 5):
        for d in np.linspace(0.02, 0.98, 25):
            lo, hi = B.feasible_p_interval(n, d)
            for p in np.linspace(0.02, 0.98, 97):
                inside = lo - 1e-9 <= p <= hi + 1e-9
                assert B.bounds_at_p(n, p, d)[2] == inside, (n, d, p)
            assert B.feasible_band(n, lo)[1] == pytest.approx(d, abs=1e-9)
            assert B.feasible_band(n, hi)[0] == pytest.approx(d, abs=1e-9)


# --- (6) the empty set is a rejection, not a wide interval ------------------
def test_empty_set_branch():
    gmin, gmax, feasible = B.bounds_at_p(3, 0.5, 1.2)
    assert not feasible
    assert gmin[0] > gmax[0]
    assert gmin[0] == pytest.approx((1.2 - (3 / 8 + 1 / 8)) / (3 / 8))
    assert B.lp_bound(3, 0.5, 1.2, 1, "max") is None
    lo, hi, share = B.bounds_over_p(3, 1.2, 0.4, 0.6)
    assert share == 0.0 and np.isnan(hi).all()


def test_partially_feasible_range_reports_its_share():
    """A transfer ratio only some of the range can produce must not be read as
    if every fraction in the range were admissible."""
    n, d = 2, 0.30                       # feasible only for p in [0.163, 0.548]
    lo_p, hi_p = B.feasible_p_interval(n, d)
    assert 0.4 < hi_p < 0.6
    _, _, share = B.bounds_over_p(n, d, 0.4, 0.6)
    assert 0.0 < share < 1.0
    assert share == pytest.approx((hi_p - 0.4) / 0.2, abs=1e-12)


# --- the panel-level reading, and the published numbers ---------------------
def test_panel_bounds_apply_to_the_effect_weighted_mean_profile():
    """b = sum_i w_i d_i, and the monotone set is convex with fixed endpoints,
    so the mean profile of any mixture of admissible profiles is admissible and
    the same bounds apply to it."""
    rng = np.random.default_rng(20261006)
    n, p = 4, 0.45
    w = B.weights(n, p)
    for _ in range(200):
        profiles = []
        for _ in range(5):
            g = np.sort(rng.uniform(0, 1, n + 1))
            g[0], g[n] = 0.0, 1.0
            profiles.append(g)
        share = rng.dirichlet(np.ones(5))
        mean = np.average(profiles, axis=0, weights=share)
        d = float(np.average([w @ g for g in profiles], weights=share))
        gmin, gmax, feasible = B.bounds_at_p(n, p, d)
        assert feasible
        assert np.all(mean[1:n] >= gmin - 1e-9)
        assert np.all(mean[1:n] <= gmax + 1e-9)


def test_W_cancels_from_the_linear_test_at_one_but_not_other_thresholds():
    """Algebraic cancellation is conditional on fixed dose scaling; no claim
    about frequentist conservatism or absence of shared normalization follows."""
    sc = pd.read_csv(HERE / "bounds.csv")
    base = B.threshold_test(sc, 218.4, 7.7, 199, 2., 1.)
    alternate = B.threshold_test(sc, 300., 35., 500, 2., 1.)
    for key in base:
        assert base[key] == pytest.approx(alternate[key])
    other_q = B.threshold_test(sc, 218.4, 7.7, 199, 2., .84)
    other_W = B.threshold_test(sc, 300., 35., 500, 2., .84)
    assert not np.allclose(other_q['normal_p'], other_W['normal_p'])


def test_the_p_free_rejection_is_exactly_mixed_below_the_paired_scale():
    sc = pd.read_csv(HERE / "bounds.csv")
    scale = json.loads((HERE / "bounds.json").read_text())["scn5a_scale"]["scale"]
    assert set(sc.loc[sc.d > 1, "variant"]) == set(sc.loc[sc.het < scale * sc.homo, "variant"])
    assert set(sc.loc[sc.d > 1, "variant"]) == {
        "D785N", "G1712C", "G1740R", "N1722D", "V1405M"}


def test_reported_counts_reproduce_from_the_table():
    sc = pd.read_csv(HERE / "bounds.csv")
    rep = json.loads((HERE / "bounds.json").read_text())
    counts = next(r for r in rep['dose_scaling_sensitivity'] if r['scale'] == rep['scn5a_scale']['scale'])
    assert counts["point_count"] == int((sc.d > 1).sum())
    assert counts['normal_nominal_set'] == sorted(sc.loc[sc.unit_normal_p < .05, 'variant'])
    assert counts['normal_bonferroni_set'] == sorted(sc.loc[sc.unit_normal_bonf < .05, 'variant'])
    assert rep["family_size"] == len(sc) == 51
    assert sum(rep["scenario_class_counts"].values()) == len(sc)
    # the headline must not outrun the correction
    assert counts['point_count'] >= len(counts['normal_nominal_set']) >= len(counts['normal_bonferroni_set'])


def test_p_star_is_consistent_with_the_stored_transfer_ratios():
    sc = pd.read_csv(HERE / "bounds.csv")
    expected = B.allele_ratio_tolerance(sc.d.to_numpy())
    assert np.allclose(sc.p_star, expected, atol=1e-12, equal_nan=True)
    inside = sc[(sc.d > 0) & (sc.d < 1)]
    assert np.allclose(B.dimer_g1(inside.p_star.to_numpy(), 0) * 0 + 1.0,
                       [B.dimer_g1(ps, dd) for ps, dd in zip(inside.p_star, inside.d)],
                       atol=1e-8)


def test_interval_bounds_keep_a_single_feasible_boundary_and_a_narrow_slice():
    """A finite grid can miss either a feasible singleton or a narrow slice."""
    lo, hi, share = B.bounds_over_p(2, .36, .6, .8)
    assert share == pytest.approx(0, abs=1e-12)
    assert lo == pytest.approx([0], abs=1e-12)
    assert hi == pytest.approx([0], abs=1e-12)
    upper = np.sqrt(.36 + 1e-8)
    lo, hi, share = B.bounds_over_p(2, .36 + 1e-8, .6, .8)
    assert share == pytest.approx((upper-.6)/.2, abs=1e-12)
    assert np.isfinite(lo).all() and np.isfinite(hi).all()


@pytest.mark.parametrize('args', [(1,.5,.4),(2,0,.4),(2,1,.4),(2,.5,np.nan)])
def test_invalid_bound_domains_fail_explicitly(args):
    with pytest.raises(ValueError):
        B.bounds_at_p(*args)


def test_ratio_confidence_set_inverts_the_normal_contrast_with_shared_W():
    sc, ctl = B.scn5a_table()
    z = ctl['critical_values']['95']
    W, Ws, s = ctl['WW'], ctl['WW_se'], ctl['scale']
    for _, row in sc.iterrows():
        for q in (row.d_95_lo, row.d_95_hi):
            contrast = (1-q)*W + q*s*row.homo-row.het
            variance = ((1-q)*Ws)**2 + (q*s*row.homo_se)**2 + row.het_se**2
            assert contrast**2/variance == pytest.approx(z*z, rel=1e-10)
        assert row.d_95_lo <= row.d <= row.d_95_hi
        assert row.g1_95_envelope_lo <= row.g1_hull_lo
        assert row.g1_95_envelope_hi >= row.g1_hull_hi


def test_uncertain_ratio_denominator_does_not_create_a_false_finite_interval():
    all_real = B.ratio_confidence_set(0, 0, 1, 1, 0, 1.96)
    assert all_real == [(-np.inf,np.inf)]
    split = B.ratio_confidence_set(0, 1, 1, .1, 0, 1.96)
    assert len(split) == 2 and split[0][0] == -np.inf and split[1][1] == np.inf
    assert split[0][1] < 0 < split[1][0]


def test_boundary_interpretation_requires_loss_of_function_sign():
    """For m>1, g1>1 implies a1>m rather than a1<m."""
    g1 = 1.2
    for m in (.5, 1.5):
        a1 = 1-g1*(1-m)
        assert (a1 < m) == (m < 1)


def test_corrected_claim_is_sensitive_to_dose_factor_and_not_an_uncalibrated_panel_bound():
    report=json.loads((HERE/'bounds.json').read_text())
    dose={r['scale']:r for r in report['dose_scaling_sensitivity']}
    assert dose[1.]['normal_bonferroni_set'] == []
    assert dose[2.]['normal_bonferroni_set'] == ['D785N']
    assert dose[2.]['welch_bonferroni_set'] == ['D785N']
    assert 'panels' not in report
    assert all(r['label']=='synthetic mathematical example' for r in report['synthetic_examples'])


def test_dose_threshold_matches_the_tabulated_scenarios():
    """S1 Fig B marks the dose factor at which the corrected one sided bound for
    p.Asp785Asn reaches zero. It must be a root of that bound, lie between the
    scenarios that do not reject (s = 1.5) and those that do (s = 2), and give
    the same decisions as the tabulated tests under both approximations."""
    from scipy import stats
    report = json.loads((HERE / "bounds.json").read_text())
    th = report["dose_threshold"]
    family = th["family_size"]
    assert family == report["family_size"] == 51
    row = pd.read_csv(HERE / "scn5a_peak.csv").set_index("variant").loc["D785N"]

    def se(s):
        return np.hypot(s * row.homo_se, row.het_se)

    def welch_df(s):
        parts = np.array([(s * row.homo_se) ** 2, row.het_se ** 2])
        return parts.sum() ** 2 / (parts ** 2 / (np.array([row.homo_n, row.het_n]) - 1)).sum()

    z = stats.norm.ppf(1 - 0.05 / family)
    assert th["critical_normal"] == pytest.approx(z, rel=1e-12)
    s_n, s_w = th["normal"], th["welch"]
    assert s_n * row.homo - row.het - z * se(s_n) == pytest.approx(0, abs=1e-8)
    t = stats.t.ppf(1 - 0.05 / family, welch_df(s_w))
    assert s_w * row.homo - row.het - t * se(s_w) == pytest.approx(0, abs=1e-8)
    assert 1.5 < s_n < s_w < 2
    dose = {r["scale"]: r for r in report["dose_scaling_sensitivity"]}
    for s, r in dose.items():
        assert ("D785N" in r["normal_bonferroni_set"]) == (s > s_n)
        assert ("D785N" in r["welch_bonferroni_set"]) == (s > s_w)
