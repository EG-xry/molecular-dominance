#!/usr/bin/env python3
"""The transfer ratio is a binomial expectation, and the design that inverts it.

    python3 run_profile.py [--quick]

Section 1c.

THE IDENTITY. Let a complex of n subunits carrying k mutant copies have
activity a_k, on a scale where the all-wild-type complex has activity 1, and
let m = a_n be the all-mutant activity. Provided m != 1, define

    g_k = (1 - a_k) / (1 - m)                                            (1)

the fraction of the variant's full effect expressed at k mutant subunits, so
g_0 = 0 and g_n = 1. Call (g_0 ... g_n) the DOMINANCE PROFILE. Equation (1) is
a reparameterisation, not an assumption: every activity vector with m != 1 has
exactly one such profile.

Under random assembly at mutant subunit fraction p, complexes with k mutant
subunits occur at P(K = k) = C(n,k) p^k (1-p)^(n-k), and the population
activity is sum_k P(K=k) a_k. With x = 1 - m the effect measured alone and
y = 1 - sum_k P(K=k) a_k the effect measured in the mixed state,

    y = 1 - sum_k P(K=k)[1 - g_k(1-m)] = (1-m) sum_k P(K=k) g_k

so that

    d = y/x = sum_k P(K=k) g_k = E[g_K],    K ~ Binomial(n, p).          (2)

The proof is two lines and the content is not depth but visibility: (2) shows
that a paired assay returns one linear functional of the profile.

CONSEQUENCES.
 (a) Non-identification. One measurement, and -- since g_0 and g_n are fixed
     by (1) -- exactly n-1 free parameters. Profiles sharing a binomial mean
     are indistinguishable.
 (b) Identification. d(p) is a polynomial of degree n in p whose Bernstein
     coefficients are the profile. Subtracting the known endpoints,

         d(p) - p^n = sum_{k=1}^{n-1} C(n,k) p^k (1-p)^(n-k) g_k         (3)

     is a linear system in n-1 unknowns, so n-1 distinct interior subunit
     fractions identify a fixed profile exactly, conditional on known p,
     random assembly, constant total complex abundance and a linear readout.
 (c) Conditional precision. Monotonicity and bounds are optional restrictions,
     not consequences of (1). This finite simulation estimates noise levels
     yielding a pooled profile RMSE near 0.1 under its stated assumptions.
     It does not prospectively establish experimental feasibility.
 (d) Dimers. n-1 = 1, so a single measurement identifies the profile if p is
     measured, and the cost of assuming p instead is a number.

Writes profile_identity.json and profile_scn5a.csv.
"""
import json
import sys
from datetime import date
from math import comb

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, minimize

QUICK = "--quick" in sys.argv
SEED = 20260929
rng = np.random.default_rng(SEED)
REP = 200 if QUICK else 600
LO, HI = 0.10, 0.90
NAMES = {2: "dimer", 3: "trimer", 4: "tetramer", 5: "pentamer", 6: "hexamer"}


def bern(n, p, lo=0, hi=None):
    hi = n if hi is None else hi
    return np.array([comb(n, k) * p ** k * (1 - p) ** (n - k)
                     for k in range(lo, hi + 1)])


def design(n, ps):
    """Rows of the interior system (3): columns k = 1 .. n-1."""
    return np.vstack([bern(n, p, 1, n - 1) for p in ps])


def fit(X, b, mode, gmax=None):
    """Least squares for the interior profile under one of three regimes.

    'free'      no constraint -- the general identifiability case.
    'bounded'   0 <= g_k <= gmax. With 0 <= m < 1, imposing 0 <= a_k <= 1
                gives 0 <= g_k <= 1/(1-m). Non-negativity alone supplies
                only the upper bound. The simulations choose gmax = 2,
                corresponding to m = 0.5, and do not restrict ordering.
    'monotone'  0 <= g_1 <= ... <= g_{n-1} <= 1, i.e. activity moves
                monotonically from 1 to m (declining when m < 1). With
                g_n = 1 fixed by definition, monotonicity already implies
                g_k <= 1, so this is ONE assumption, not two.

    The monotone problem uses explicit linear inequalities in g. An infeasible
    or unsuccessful numerical result is rejected rather than silently used.
    """
    m = X.shape[1]
    if mode == "free":
        return np.linalg.lstsq(X, b, rcond=None)[0]
    if mode == "bounded":
        hi = np.inf if gmax is None else gmax
        result = lsq_linear(X, b, bounds=(0.0, hi), max_iter=300)
        if not result.success:
            raise RuntimeError(f"Bounded profile fit failed: {result.message}")
        return result.x
    if mode != "monotone":
        raise ValueError(f"Unknown profile fit mode: {mode}")
    # Solve the actual constrained problem rather than solving it free and
    # clipping the answer: 0 <= g_1 <= ... <= g_{n-1} <= 1 as explicit linear
    # inequalities, minimising ||Xg - b||^2 under them.
    D = np.eye(m) - np.eye(m, k=-1)          # g_k - g_{k-1} >= 0, with g_0 = 0
    cons = [{"type": "ineq", "fun": lambda g, D=D: D @ g,
             "jac": lambda g, D=D: D},
            {"type": "ineq", "fun": lambda g: 1.0 - g[-1],
             "jac": lambda g: -np.eye(m)[-1]}]
    g0 = np.maximum.accumulate(
        np.clip(np.linalg.lstsq(X, b, rcond=None)[0], 0.0, 1.0))
    r = minimize(lambda g: float(((X @ g - b) ** 2).sum()), g0,
                 jac=lambda g: 2 * X.T @ (X @ g - b), constraints=cons,
                 bounds=[(0.0, 1.0)] * m, method="SLSQP",
                 options={"maxiter": 200, "ftol": 1e-12})
    if (not r.success or not np.isfinite(r.x).all()
            or np.min(D @ r.x) < -1e-8 or r.x[-1] > 1 + 1e-8):
        raise RuntimeError(f"Monotone profile fit failed: {r.message}")
    return r.x


def rmse(n, ps, sd, rep, mode="monotone", gmax=None, profile=None):
    """POOLED root mean squared error: squared errors are accumulated across
    every replicate and every profile entry before the square root is taken.

    An earlier version took the root inside the replicate loop and averaged the
    results.  By Jensen that is at most the pooled RMSE, and for a dimer -- one
    unknown -- it is exactly the mean absolute error, so it reported a
    different and more flattering statistic than the text claimed.
    """
    X = design(n, ps)
    sq = []
    for _ in range(rep):
        if profile is None:
            g = np.sort(rng.uniform(0, 1, n + 1))
            g[0], g[n] = 0.0, 1.0
            t = g[1:n]
        else:
            t = profile
        b = X @ t + rng.normal(0, sd, len(ps))
        sq.append((fit(X, b, mode, gmax) - t) ** 2)
    return float(np.sqrt(np.concatenate(sq).mean()))



def a1_at(p, mix, m):
    return (mix - (1 - p) ** 2 - p ** 2 * m) / (2 * p * (1 - p))



def main():
    global rng
    rng = np.random.default_rng(SEED)
    # --- (1) the identity, against a direct simulation of the mixture ----------
    worst = 0.0
    for _ in range(20000 if not QUICK else 2000):
        n = int(rng.integers(2, 8))
        p = float(rng.uniform(0.02, 0.98))
        m = float(rng.uniform(0.0, 0.95))
        g = np.sort(rng.uniform(0, 1, n + 1))
        g[0], g[n] = 0.0, 1.0
        a = 1 - g * (1 - m)
        P = bern(n, p)
        worst = max(worst, abs((1 - float(P @ a)) / (1 - m) - float(P @ g)))

    # --- (2) exact identification from n-1 ratios ------------------------------
    exact = {}
    for n in range(2, 7):
        m = n - 1
        ps = np.array([0.5]) if m == 1 else np.linspace(LO, HI, m)
        X = design(n, ps)
        err = 0.0
        for _ in range(200):
            g = np.sort(rng.uniform(0, 1, n + 1))
            g[0], g[n] = 0.0, 1.0
            b = X @ g[1:n]
            err = max(err, float(np.abs(np.linalg.solve(X, b) - g[1:n]).max()))
        exact[n] = {"n_unknowns": m, "n_ratios": m,
                    "design": [round(float(v), 3) for v in ps],
                    "max_recovery_error": err,
                    "condition_number": float(np.linalg.cond(X))}

    # --- (3) precision required, with the constraint ---------------------------


    precision = {}
    for n in range(2, 7):
        m = n - 1
        ps = np.array([0.5]) if m == 1 else np.linspace(LO, HI, m)
        row = {"complex": NAMES[n], "n_ratios": m}
        for mode, gmax in (("free", None), ("bounded", 2.0), ("monotone", None)):
            lo, hi = 1e-5, 0.60
            for _ in range(16 if QUICK else 22):
                mid = (lo + hi) / 2
                if rmse(n, ps, mid, REP, mode, gmax) > 0.1:
                    hi = mid
                else:
                    lo = mid
            row[f"se_{mode}"] = round((lo + hi) / 2, 4)
        # Monte Carlo spread of the monotone figure across seeds
        seeds = []
        for sd_seed in range(4 if QUICK else 8):
            g_rng = np.random.default_rng(1000 + sd_seed)
            saved = rng
            globals()["rng"] = g_rng
            lo, hi = 1e-5, 0.60
            for _ in range(14):
                mid = (lo + hi) / 2
                if rmse(n, ps, mid, max(120, REP // 3), "monotone") > 0.1:
                    hi = mid
                else:
                    lo = mid
            seeds.append((lo + hi) / 2)
            globals()["rng"] = saved
        row["se_monotone_seed_range"] = [round(min(seeds), 3), round(max(seeds), 3)]
        precision[n] = row

    # --- (4) additional monotone profiles, not an exhaustive robustness test ---
    # The finite set of step profiles is within the monotone model. Performance
    # here cannot establish robustness to nonmonotone or changing profiles.
    steps = []
    for n in (4, 5):
        m = n - 1
        ps = np.linspace(LO, HI, m)
        X = design(n, ps)
        for t in range(1, n + 1):
            prof = np.array([1.0 if k >= t else 0.0 for k in range(n + 1)])
            if prof[0] != 0 or prof[n] != 1:
                continue
            true = prof[1:n]
            su, sc = [], []
            for _ in range(REP):
                b = X @ true + rng.normal(0, 0.05, m)
                su.append((np.linalg.solve(X, b) - true) ** 2)
                sc.append((fit(X, b, "monotone") - true) ** 2)
            ru = float(np.sqrt(np.concatenate(su).mean()))
            rc = float(np.sqrt(np.concatenate(sc).mean()))
            steps.append({"n": n, "tolerance_t": t,
                          "rmse_unconstrained": round(ru, 4),
                          "rmse_constrained": round(rc, 4),
                          "gain": round(ru / rc, 2)})

    # A valid fixed-activity random-assembly system outside the monotone model.
    # Every activity is nonnegative and <= 1; monotonicity still fails. There is
    # no measurement noise, so the constrained error is model misspecification.
    counter_n, counter_m = 3, 0.5
    counter_g = np.array([0.0, 1.6, 0.2, 1.0])
    counter_a = 1 - (1 - counter_m) * counter_g
    counter_ps = np.array([0.1, 0.9])
    counter_X = design(counter_n, counter_ps)
    counter_b = counter_X @ counter_g[1:-1]
    counter_free = fit(counter_X, counter_b, "free")
    counter_monotone = fit(counter_X, counter_b, "monotone")
    counterexample = {
        "n": counter_n, "m": counter_m, "p": counter_ps.tolist(),
        "profile": counter_g.tolist(), "activities": counter_a.tolist(),
        "true_transfer_ratios": (counter_b + counter_ps ** counter_n).tolist(),
        "free_estimate": counter_free.tolist(),
        "monotone_estimate": counter_monotone.tolist(),
        "free_profile_rmse": float(np.sqrt(np.mean(
            (counter_free - counter_g[1:-1]) ** 2))),
        "monotone_profile_rmse": float(np.sqrt(np.mean(
            (counter_monotone - counter_g[1:-1]) ** 2))),
        "interpretation": "Noiseless misspecification error under a monotone restriction; "
                          "the unrestricted identity and inversion still hold."
    }

    # --- (5) historical SCN5A unmatched-scale sensitivity ---------------------
    # Retained as a labeled sensitivity; the displayed paired-scale results
    # now come from run_bounds.py (W=218.4, s=2.184, p=.542125).
    d = pd.read_csv("scn5a_specifications.csv")
    d = d[d.variant != "W822X"].copy()
    WW = float(pd.read_csv("scn5a_controls.csv").set_index("variant").loc["WT+WT", "het"])
    d["m"] = 2 * d.homo / WW
    d["mix"] = d.het / WW
    legacy_scope = "historical unmatched-scale sensitivity: s=2, W=218.4; displayed paired results are in bounds.csv"
    d["analysis_scope"] = legacy_scope
    d["combined_scale"] = 2.
    d["reference"] = WW




    for p in (0.40, 0.45, 0.50, 0.55, 0.60):
        tag = int(p * 100)
        d[f"a1_p{tag}"] = a1_at(p, d.mix, d.m)
        d[f"g1_p{tag}"] = (1 - d[f"a1_p{tag}"]) / (1 - d.m)
    d["requires_negative_a1"] = d.a1_p50 < 0
    d["worse_than_mutant_homodimer"] = d.a1_p50 < d.m
    d["g1_below_zero"] = d.g1_p50 < 0
    d["g1_above_one"] = d.g1_p50 > 1
    d["outside_monotone_profile"] = d.g1_below_zero | d.g1_above_one
    d[["variant", "homo", "het", "m", "mix", "in_intersection", "a1_p50",
       "g1_p40", "g1_p45", "g1_p50", "g1_p55", "g1_p60",
       "requires_negative_a1", "worse_than_mutant_homodimer",
       "g1_below_zero", "g1_above_one", "outside_monotone_profile",
       "analysis_scope", "combined_scale", "reference"]].to_csv(
        "profile_scn5a.csv", index=False)

    x = (WW - 2 * d.homo).values
    y = (WW - d.het).values
    d_panel = float(np.sum(x * y) / np.sum(x * x))
    g1_by_p = {p: float((d_panel - p ** 2) / (2 * p * (1 - p)))
               for p in (0.40, 0.45, 0.50, 0.55, 0.60)}
    cost = abs(g1_by_p[0.40] - g1_by_p[0.60])

    # what transfer-ratio precision does this assay actually deliver?
    se_d = np.sqrt((d.het_se / x) ** 2 + ((y / x ** 2) * 2 * d.homo_se) ** 2)
    dimer_req = precision[2]["se_monotone"]
    tetramer_req = precision[4]["se_monotone"]

    json.dump({"generated": date.today().isoformat(),
               "simulation_settings": {
                   "mode": "quick" if QUICK else "full",
                   "seed": SEED,
                   "identity_draws": 2000 if QUICK else 20000,
                   "identity_n_range_inclusive": [2, 7],
                   "identity_p_uniform_range": [0.02, 0.98],
                   "identity_m_uniform_range": [0.0, 0.95],
                   "identification_draws_per_n": 200,
                   "precision_n_range_inclusive": [2, 6],
                   "replicates_per_noise_level": REP,
                   "noise_search_iterations": 16 if QUICK else 22,
                   "noise_search_bounds": [1e-5, 0.60],
                   "target_pooled_profile_rmse": 0.1,
                   "profile_sampling": "Sort n+1 independent Uniform(0,1) draws, "
                                       "then replace endpoints by 0 and 1.",
                   "noise": "Independent homoscedastic Gaussian noise on d(p); "
                            "new profiles and errors at each bisection evaluation.",
                   "bounded_gmax": 2.0,
                   "bounded_interpretation": "0 <= a_k <= 1 with m = 0.5; "
                                             "not merely nonnegative activities.",
                   "seed_range_seeds": list(range(1000, 1000 + (4 if QUICK else 8))),
                   "seed_range_replicates_per_noise_level": max(120, REP // 3),
                   "seed_range_search_iterations": 14,
                   "step_profile_noise_sd": 0.05,
                   "step_profile_n": [4, 5],
                   "conditional_assumptions": [
                       "Known n and n-1 distinct p strictly between 0 and 1.",
                       "Random independent subunit assembly at measured p.",
                       "A single fixed activity/profile across all p.",
                       "Constant total complex abundance and linear activity readout.",
                       "Known endpoints g_0=0, g_n=1 and m != 1."
                   ],
                   "omitted_uncertainty": [
                       "No p measurement error or uncertainty in complex abundance.",
                       "No covariance from shared references or ratio denominators.",
                       "No profile drift, heteroscedasticity, or readout nonlinearity."
                   ],
                   "scope": "Finite conditional Monte Carlo noise thresholds, not "
                            "prospective experimental validation or confidence limits."
               },
               "identity_max_error": worst,
               "exact_identification": exact,
               "precision_required": precision,
               "step_profile_robustness": steps,
               "nonmonotone_counterexample": counterexample,
               "scn5a": {"analysis_scope": legacy_scope, "combined_scale": 2., "reference": WW,
                         "panel_transfer_ratio": d_panel, "g1_by_p": g1_by_p,
                         "cost_of_p_uncertainty_0.10": cost,
                         "assay_se_median": float(se_d.median()),
                         "assay_se_iqr": [float(se_d.quantile(.25)),
                                          float(se_d.quantile(.75))],
                         "assay_se_scope": "First-order per-variant ratio SE assuming "
                             "independent homo/het errors and a fixed WT+WT reference; "
                             "reference uncertainty and covariance omitted. Fractions "
                             "below simulation thresholds are descriptive comparisons, "
                             "not validation of a multi-ratio recovery experiment.",
                         "fraction_meeting_dimer_requirement":
                             float((se_d <= dimer_req).mean()),
                         "fraction_meeting_tetramer_requirement":
                             float((se_d <= tetramer_req).mean()),
                         "dimer_requirement": dimer_req,
                         "tetramer_requirement": tetramer_req,
                         "n_violating_monotone_assumption":
                             int(d.outside_monotone_profile.sum()),
                         "n_g1_above_one": int(d.g1_above_one.sum()),
                         "n_g1_below_zero": int(d.g1_below_zero.sum()),
                         "n_requiring_negative_a1": int(d.requires_negative_a1.sum()),
                         "n_worse_than_mutant_homodimer":
                             int(d.worse_than_mutant_homodimer.sum()),
                         "n_panel": len(d)}},
              open("profile_identity.json", "w"), indent=2)

    print(f"1. IDENTITY  d = E[g_K]: max error {worst:.2e} over random (n,p,m,profile)\n")
    print("2. EXACT IDENTIFICATION from n-1 ratios (g_0, g_n are fixed by definition)")
    print(f"   {'n':>2s} {'unknowns':>9s} {'ratios':>7s} {'recovery err':>13s} {'cond':>7s}")
    for n, v in exact.items():
        print(f"   {n:2d} {v['n_unknowns']:9d} {v['n_ratios']:7d} "
              f"{v['max_recovery_error']:13.2e} {v['condition_number']:7.1f}")
    print("\n3. CONDITIONAL NOISE THRESHOLDS for POOLED profile RMSE near 0.1")
    print(f"   {'complex':>10s} {'ratios':>7s} {'free':>8s} {'0<=g<=2':>8s} "
          f"{'monotone':>9s}  {'monotone, seed range':>21s}")
    for n, v in precision.items():
        sr = v["se_monotone_seed_range"]
        print(f"   {v['complex']:>10s} {v['n_ratios']:7d} {v['se_free']:8.3f} "
              f"{v['se_bounded']:8.3f} {v['se_monotone']:9.3f}  "
              f"{f'{sr[0]:.3f} - {sr[1]:.3f}':>21s}")
    print("\n4. FINITE STEP-PROFILE SENSITIVITY within the monotone model:")
    print(f"   gains {min(s['gain'] for s in steps):.1f}x to {max(s['gain'] for s in steps):.1f}x "
          "for these step profiles; no guarantee outside the model.")
    print("   Nonmonotone noiseless counterexample: constrained profile RMSE "
          f"{counterexample['monotone_profile_rmse']:.3f}; "
          f"free recovery RMSE {counterexample['free_profile_rmse']:.2e}.")
    print(f"\n5. SCN5A, a dimer: panel transfer ratio d = {d_panel:.4f}")
    for p, g in g1_by_p.items():
        print(f"     p = {p:.2f}  ->  g_1 = {g:6.3f}")
    print(f"   p from 0.40 to 0.60 spans {cost:.2f} in g_1; bounds are not universal")
    print(f"   assay SE on the transfer ratio: median {se_d.median():.3f} "
          f"(IQR {se_d.quantile(.25):.3f}-{se_d.quantile(.75):.3f})")
    print(f"   {(se_d <= dimer_req).mean():.0%} fall below the conditional DIMER threshold "
          f"({dimer_req:.3f}); {(se_d <= tetramer_req).mean():.0%} fall below the "
          f"TETRAMER threshold ({tetramer_req:.3f})")
    print("   These SEs condition on a fixed reference and independent measurement errors.")
    print(f"   {int(d.outside_monotone_profile.sum())}/{len(d)} point estimates lie outside "
          f"the monotone profile: {int(d.g1_above_one.sum())} above 1 and "
          f"{int(d.g1_below_zero.sum())} below 0.")
    print("   These are point-estimate model discrepancies, not significance tests.")
    print(f"   {int(d.requires_negative_a1.sum())}/{len(d)} require a_1 < 0; "
          f"{int(d.worse_than_mutant_homodimer.sum())}/{len(d)} require a_1 < m")
    print("\nwrote profile_identity.json, profile_scn5a.csv")


if __name__ == "__main__":
    main()
