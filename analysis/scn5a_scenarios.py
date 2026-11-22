"""Pure calculations for the SCN5A control, replication and scale scenarios.

Every test conditions on independently sampled cells and independent control
means. Unknown batch structure and shared-normalization covariance are not
estimated by these sensitivity calculations.
"""
from math import sqrt

import numpy as np
from scipy import stats
from scipy.optimize import brentq

TEST_NAMES = (
    "fixed floor, two-sided normal",
    "floor uncertainty propagated, two-sided normal",
    "floor uncertainty propagated, one-sided normal",
    "Welch-Satterthwaite, reported cell counts",
    "log-ratio delta method, two-sided normal",
)


def control_floor(reference, reference_se, surrogate, surrogate_se=0.,
                  reference_n=199, surrogate_n=None):
    """F=B²/W and p=1−B/W, with independent means B and W.

    B is either the normalized single-WT reference (fixed at 100) or the
    W822X mixed response, conditionally treated as an inert-product proxy.
    Welch degrees of freedom follow the two derivative variance components.
    """
    floor = surrogate ** 2 / reference
    parts = np.array([(2 * surrogate / reference * surrogate_se) ** 2,
                      (floor / reference * reference_se) ** 2])
    ns = np.array([surrogate_n or 2, reference_n], dtype=float)
    variance = float(parts.sum())
    dof = variance ** 2 / np.sum(parts ** 2 / (ns - 1)) if variance else float("inf")
    return {"p": 1 - surrogate / reference, "floor": floor,
            "floor_se": sqrt(variance), "floor_df": float(dof),
            "surrogate_variance_component": float(parts[0]),
            "reference_variance_component": float(parts[1])}


def test_floor(frame, floor, floor_se, floor_df=198., family=None):
    family = len(frame) if family is None else family
    variance = frame.het_se.to_numpy(float) ** 2
    gap = floor - frame.het.to_numpy(float)
    z_fixed = gap / np.sqrt(variance)
    z_prop = gap / np.sqrt(variance + floor_se ** 2)
    dof = (variance + floor_se ** 2) ** 2 / (
        variance ** 2 / (frame.het_n.to_numpy(float) - 1)
        + floor_se ** 4 / floor_df)
    z_log = np.log(floor / frame.het) / np.sqrt(
        (frame.het_se / frame.het) ** 2 + (floor_se / floor) ** 2)
    p = (2 * stats.norm.sf(z_fixed), 2 * stats.norm.sf(z_prop),
         stats.norm.sf(z_prop), 2 * stats.t.sf(z_prop, dof),
         2 * stats.norm.sf(z_log))
    adjusted = {name: np.minimum(1., family * raw) for name, raw in zip(TEST_NAMES, p)}
    sets = {name: sorted(frame.loc[q < .05, "variant"]) for name, q in adjusted.items()}
    intersection = sorted(set.intersection(*(set(v) for v in sets.values())))
    return {"counts": {name: len(v) for name, v in sets.items()}, "sets": sets,
            "intersection": intersection,
            "adjusted_p_values": {name: {v: float(q) for v, q in zip(frame.variant, values)}
                                  for name, values in adjusted.items()}}


def critical_effective_df(statistic, family):
    normal = stats.norm.isf(.05 / (2 * family))
    if not np.isfinite(statistic) or statistic <= normal:
        return None
    return float(brentq(lambda nu: stats.t.isf(.05 / (2 * family), nu) - statistic,
                        1., 1e8))


def robustness_rows(frame, scenarios, effective_dfs=(5, 7, 9, 10, 15, 19, 29, 49, 99, 199)):
    """Return scenario-specific tipping points and joint SE/df boundaries.

    The main inflation multiplier acts on all input SEs simultaneously.
    lambda_max=z/tcrit(df) is a strict boundary: rejection holds for a smaller
    multiplier. Values below one mean the original SEs already fail at that df.
    A second boundary inflates only the variant SE, leaving control SE fixed.
    Discrepancy SDs enter in quadrature, in original current percentage points.
    """
    family = len(frame)
    zc = stats.norm.isf(.05 / (2 * family))
    rows, surface = [], []
    for scenario, sc in scenarios.items():
        floor, floor_se = sc["floor"], sc["floor_se"]
        for r in frame.itertuples():
            gap = floor - r.het
            sigma_prop = sqrt(r.het_se ** 2 + floor_se ** 2)
            sigma_log = sqrt((r.het_se / r.het) ** 2 + (floor_se / floor) ** 2)
            z_prop, z_log = gap / sigma_prop, np.log(floor / r.het) / sigma_log
            prop_need = (gap / zc) ** 2 - sigma_prop ** 2
            log_need = (np.log(floor / r.het) / zc) ** 2 - sigma_log ** 2
            row = {"scenario": scenario, "variant": r.variant, "floor": floor,
                   "floor_se": floor_se, "mixed_mean": r.het, "mixed_se": r.het_se,
                   "mixed_cells": int(r.het_n), "z_propagated": float(z_prop),
                   "z_log_ratio": float(z_log),
                   "tau_star_propagated": sqrt(prop_need) if z_prop > zc else None,
                   "tau_star_log_ratio": floor * sqrt(log_need) if z_log > zc else None,
                   "effective_df_star_propagated": critical_effective_df(z_prop, family),
                   "effective_df_star_log_ratio": critical_effective_df(z_log, family),
                   "normal_max_all_se_inflation_log_ratio": float(z_log / zc),
                   "in_five_test_intersection": r.variant in sc["intersection"]}
            rows.append(row)
            for df in (*effective_dfs, float("inf")):
                critical = zc if np.isinf(df) else stats.t.isf(.05 / (2 * family), df)
                prop_variant = (gap / critical) ** 2 - floor_se ** 2
                log_variant = (np.log(floor / r.het) / critical) ** 2 - (floor_se / floor) ** 2
                surface.append({"scenario": scenario, "variant": r.variant,
                                "effective_df": None if np.isinf(df) else df,
                                "effective_df_label": "normal" if np.isinf(df) else str(df),
                                "critical_value": float(critical),
                                "max_all_se_inflation_propagated": max(0., float(z_prop / critical)),
                                "max_all_se_inflation_log_ratio": max(0., float(z_log / critical)),
                                "max_variant_se_inflation_propagated": sqrt(prop_variant) / r.het_se if gap > 0 and prop_variant >= 0 else 0.,
                                "max_variant_se_inflation_log_ratio": sqrt(log_variant) / (r.het_se / r.het) if floor > r.het and log_variant >= 0 else 0.})
    return rows, surface


def robustness_counts(frame, floor, floor_se, effective_dfs=(9,),
                      discrepancy_sds=(0., 10., 15., 20.)):
    """Counts and variant sets at a specified floor, across two sensitivities.

    Effective-df scenarios keep the reported log-ratio SE unchanged. Additional
    discrepancy enters the propagated variance as tau² and the log variance as
    (tau/F)², matching robustness_rows. Every scenario uses all variants for
    the same two-sided Bonferroni family.
    """
    family = len(frame)
    gap = floor - frame.het.to_numpy(float)
    log_gap = np.log(floor / frame.het.to_numpy(float))
    var_prop = frame.het_se.to_numpy(float) ** 2 + floor_se ** 2
    var_log = (frame.het_se / frame.het).to_numpy(float) ** 2 + (floor_se / floor) ** 2
    normal_critical = stats.norm.isf(.05 / (2 * family))
    df_rows, discrepancy_rows = [], []
    for df in effective_dfs:
        selected = sorted(frame.loc[log_gap / np.sqrt(var_log) >
                                    stats.t.isf(.05 / (2 * family), df), "variant"])
        df_rows.append({"effective_df": df, "log_ratio_count": len(selected),
                        "log_ratio_variants": selected})
    for tau in discrepancy_sds:
        propagated = sorted(frame.loc[gap / np.sqrt(var_prop + tau ** 2) >
                                      normal_critical, "variant"])
        log_ratio = sorted(frame.loc[log_gap / np.sqrt(var_log + (tau / floor) ** 2) >
                                     normal_critical, "variant"])
        discrepancy_rows.append({"additional_sd": tau,
                                 "propagated_normal_count": len(propagated),
                                 "propagated_normal_variants": propagated,
                                 "log_ratio_count": len(log_ratio),
                                 "log_ratio_variants": log_ratio})
    return {"floor": floor, "floor_se": floor_se, "family_size": family,
            "effective_df_scenarios": df_rows,
            "discrepancy_scenarios": discrepancy_rows}


def floor_interval_for_intersection(frame, target, relative_floor_se):
    """Exact floor interval giving the target five-test intersection.

    Vary F while keeping s_F/F fixed. Every constituent test is monotone in F
    on the below-floor branch. Find its entry threshold for each variant,
    then intersect the five test decisions. Boundaries use strict P<.05.
    """
    family = len(frame)
    entries = {}
    for _, row in frame.iterrows():
        one = frame.loc[[row.name]]
        thresholds = []
        for test in TEST_NAMES:
            def root(floor):
                q = test_floor(one, floor, relative_floor_se * floor, family=family)["adjusted_p_values"][test][row.variant]
                return q - .05
            thresholds.append(float(brentq(root, float(row.het), 1e5)))
        entries[row.variant] = max(thresholds)
    lo = max(entries[v] for v in target)
    upper_variant = min((v for v in entries if v not in target), key=entries.get)
    hi = entries[upper_variant]
    return {"lower_exclusive": lo, "upper_inclusive": hi,
            "relative_floor_se": relative_floor_se,
            "target": sorted(target), "next_variant": upper_variant,
            "entry_floors": entries,
            "limits": "Interval holds only for the specified relative floor SE, family51, reported cell counts and five tests; it is not a confidence interval for the floor."}
