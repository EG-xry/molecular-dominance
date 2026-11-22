#!/usr/bin/env python3
"""SCN5A specification-sensitivity driver.

One command emits every specification statistic quoted in the manuscript, the
intersection across specifications, the candidate table, and Figure 1.

    python3 run_scn5a.py

Input   scn5a_peak.csv      50 missense variants and W822X, and scn5a_controls.csv
                            (WT, WT+WT).  Both extracted by extract_scn5a.py from
                            O'Neill et al. 2022 Genet Med, supplementary
                            workbook NIHMS1819453-supplement-ONeill_GiM_
                            Supplemental_xlsx.xlsx (peak current density,
                            % of wild type).  Columns: variant, homo, homo_se,
                            homo_n (variant alone), het, het_se, het_n
                            (wild type + variant co-expression).
Output  scn5a_specifications.csv    per-variant statistics, every specification
        scn5a_specifications.json   counts, intersection, provenance
        Figure1.png / Figure1.pdf

Every number in the manuscript's section 4, its abstract, and Figure 1 comes
from this script's output files and from nowhere else.
"""
import json
import sys
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

from hgvs import scn5a_hgvs
from scn5a_scenarios import (control_floor, floor_interval_for_intersection,
                            robustness_counts, robustness_rows, test_floor)

CSV = "scn5a_peak.csv"
CONTROLS = "scn5a_controls.csv"
SEED_INFO = "no randomisation: all statistics here are closed-form"

# --------------------------------------------------------------------------
# Controls.  These are read from the file, not asserted, so that a change in
# the extraction propagates rather than silently disagreeing with the text.
# --------------------------------------------------------------------------
ctl = pd.read_csv(CONTROLS).set_index("variant")
WT = float(ctl.loc["WT", "het"])                      # 100.0 by normalisation
WW, WW_SE, WW_N = (float(ctl.loc["WT+WT", c]) for c in ("het", "het_se", "het_n"))

d = pd.read_csv(CSV).reset_index(drop=True)
FAMILY = len(d)                                        # Bonferroni family size

# --------------------------------------------------------------------------
# The floor.
#
# Under (i) equal abundance of the two alleles, (ii) random dimer assembly,
# (iii) total complex abundance unchanged from the WT+WT condition, and
# (iv) wild-type homodimer activity unchanged by the presence of the mutant,
# a dead-mutant dimer that kills every heterodimer leaves only the WW quarter
# of the pool:  floor = WW / 2**n_subunits, with n_subunits = 2.
#
# The floor therefore inherits the WT+WT standard error, scaled identically.
# It is a *conditional* prediction of that model, not a measured quantity.
# --------------------------------------------------------------------------
def floor_for(n_subunits: int) -> tuple[float, float, int]:
    k = 2 ** n_subunits
    return WW / k, WW_SE / k, int(WW_N)


FLOOR, FLOOR_SE, FLOOR_N = floor_for(2)                # 54.6 +- 1.925, n = 199

# --------------------------------------------------------------------------
# Specifications.  Each returns a boolean Series: does this variant's mixed
# state sit significantly BELOW the model's floor after Bonferroni correction
# across the 51-variant family?
#
# Written out explicitly because the previous draft described them in prose
# and two of them were then irreproducible (review finding 2).
# --------------------------------------------------------------------------
def bonf(p):
    return np.minimum(1.0, FAMILY * np.asarray(p, dtype=float))


def spec_fixed_floor(df, floor=FLOOR):
    """Floor treated as a known constant; only the variant's SE enters.
    z = (floor - het) / het_se,  two-sided normal."""
    z = (floor - df.het) / df.het_se
    return z, bonf(2 * stats.norm.sf(z))


def spec_propagated_two(df, floor=FLOOR, fse=FLOOR_SE):
    """Floor uncertainty propagated in quadrature, two-sided normal.
    z = (floor - het) / sqrt(het_se^2 + floor_se^2)."""
    z = (floor - df.het) / np.sqrt(df.het_se ** 2 + fse ** 2)
    return z, bonf(2 * stats.norm.sf(z))


def spec_propagated_one(df, floor=FLOOR, fse=FLOOR_SE):
    """As above, one-sided: the model predicts a floor, so only breaches below
    it are evidence against it."""
    z = (floor - df.het) / np.sqrt(df.het_se ** 2 + fse ** 2)
    return z, bonf(stats.norm.sf(z))


def spec_welch(df, floor=FLOOR, fse=FLOOR_SE, fn=FLOOR_N):
    """Welch's t with Satterthwaite degrees of freedom, using the reported
    cell counts for BOTH variances:

        df = (s1^2 + s2^2)^2 / ( s1^4/(n1-1) + s2^4/(n2-1) )

    The previous draft used het_n - 1 alone, which ignores the floor's own
    198 degrees of freedom and understates df; that is the source of the
    6-vs-7 discrepancy in review finding 2.  Valid only if cells are
    independent and the normalisation adds no variance."""
    s1, s2 = df.het_se ** 2, fse ** 2
    nu = (s1 + s2) ** 2 / (s1 ** 2 / (df.het_n - 1) + s2 ** 2 / (fn - 1))
    t = (floor - df.het) / np.sqrt(s1 + s2)
    return t, bonf(2 * stats.t.sf(t, nu)), nu


def spec_log_ratio(df, floor=FLOOR, fse=FLOOR_SE):
    """Delta method on log(floor / het), two-sided normal:

        z = log(floor/het) / sqrt( (het_se/het)^2 + (floor_se/floor)^2 )

    Each standard error is divided by its OWN mean.  The previous draft
    divided the floor's SE by het, which is not the delta-method variance;
    that is the source of the 3-vs-4 discrepancy in review finding 2.
    This specification is the relevant one if the assay's multiplicative
    error dominates, which percentage-of-wild-type readouts often show."""
    z = np.log(floor / df.het) / np.sqrt((df.het_se / df.het) ** 2 + (fse / floor) ** 2)
    return z, bonf(2 * stats.norm.sf(z))


# --------------------------------------------------------------------------
# Run them all.
# --------------------------------------------------------------------------
out = d.copy()
counts, sets = {}, {}

z, p = spec_fixed_floor(out);      out["z_fixed"], out["p_fixed"] = z, p
z, p = spec_propagated_two(out);   out["z_prop"], out["p_prop_two"] = z, p
_, p = spec_propagated_one(out);   out["p_prop_one"] = p
t, p, nu = spec_welch(out);        out["t_welch"], out["welch_df"], out["p_welch"] = t, nu, p
z, p = spec_log_ratio(out);        out["z_log"], out["p_log"] = z, p

PRIMARY = {
    "fixed floor, two-sided normal":                 "p_fixed",
    "floor uncertainty propagated, two-sided normal": "p_prop_two",
    "floor uncertainty propagated, one-sided normal": "p_prop_one",
    "Welch-Satterthwaite, reported cell counts":      "p_welch",
    "log-ratio delta method, two-sided normal":       "p_log",
}
for label, col in PRIMARY.items():
    hit = out.loc[out[col] < 0.05, "variant"]
    counts[label], sets[label] = len(hit), set(hit)

# Alternative assembly stoichiometries, fixed-floor specification.  These are
# reported as a sensitivity on the assembly model, not as competing tests.
for n_sub, name in ((3, "trimer"), (4, "tetramer")):
    fl, _, _ = floor_for(n_sub)
    z, p = spec_fixed_floor(out, floor=fl)
    out[f"p_fixed_{name}"] = p
    hit = out.loc[p < 0.05, "variant"]
    counts[f"{name} floor ({fl:.2f}%), fixed-floor"] = len(hit)
    sets[f"{name} floor ({fl:.2f}%), fixed-floor"] = set(hit)

INTERSECTION = sorted(set.intersection(*(sets[k] for k in PRIMARY)))
out["in_intersection"] = out.variant.isin(INTERSECTION)
# HGVS names on the MANE Select reference NP_000326.2 (source labels follow NP_932173.1).
out["hgvs_p_NP_000326_2"] = out.variant.map(scn5a_hgvs)

# Illustrative departures from a strict additive prediction. These four
# controls do not sample the error distribution of the conservative assembly
# floor. W822X retains 4.7% current alone; protein absence was not measured.
# The WW/2 + alone prediction additionally assumes linear dose scaling.
CONTROLS = ["W822X", "G1661R", "S1672Y", "R893C"]
resid = [float(out.loc[out.variant == v, "het"].iloc[0])
         - (WW / 2 + float(out.loc[out.variant == v, "homo"].iloc[0]))
         for v in CONTROLS if (out.variant == v).any()]
rms_model_error = float(np.sqrt(np.mean(np.square(resid)))) if resid else np.nan

# Sensitivity: add a hypothetical independent Gaussian discrepancy SD.
# These SDs are analyst-chosen scenarios, not estimated from the four controls.
model_error_rows = {}
for extra in (0, 5, 10, 15, 20):
    se = np.sqrt(FLOOR_SE ** 2 + extra ** 2)
    z = (FLOOR - out.het) / np.sqrt(out.het_se ** 2 + se ** 2)
    hit = out.loc[bonf(2 * stats.norm.sf(z)) < 0.05, "variant"]
    model_error_rows[extra] = sorted(hit)

# Descriptive mixed-minus-alone contrast. A mechanistic interpretation still
# assumes the mutant contribution is unchanged across conditions. The test
# also assumes independent errors and omits normalization/batch covariance.
mf = out[out.het < out.homo].copy()
mf["diff"] = mf.homo - mf.het
mf["diff_se"] = np.sqrt(mf.homo_se ** 2 + mf.het_se ** 2)
mf["z"] = mf["diff"] / mf["diff_se"]
mf["p"] = 2 * stats.norm.sf(mf.z.abs())

# Strict-additivity floor.  WT+WT = 218.4 is 2.39 SE ABOVE the 200 that strict
# additivity predicts, so anchoring the floor to the observed control rather
# than to 200 raises it from 50.0 to 54.6 -- the direction that produces MORE
# breaches.  Both are reported.
STRICT = 2 * WT / 4
z_strict = (STRICT - out.het) / out.het_se
out["p_strict"] = bonf(2 * stats.norm.sf(z_strict))
strict_hits = sorted(out.loc[out.p_strict < 0.05, "variant"])
supra = (WW - 2 * WT) / WW_SE

# Falconer mapping.  Falconer's x is the deviation of the mutant HOMOZYGOTE
# from the wild-type homozygote.  The assay measures the variant at ONE allele
# against a single-WT reference of 100, while the mixed state carries TWO
# alleles against 218.4, and peak current is extensive -- so the one-allele
# value cannot serve as 2a directly.  A hypothetical linear-dose model would rescale it to a two-allele
# equivalent; WT+WT does not validate that model for mutant constructs:
#     x = WW - 2*homo      (mutant homozygote, two-allele scale)
#     y = WW - het         (heterozygote,      two-allele scale)
fal = out[out.variant != "W822X"]
fx, fy = (WW - 2 * fal.homo).values, (WW - fal.het).values
b_fal = float(np.sum(fx * fy) / np.sum(fx * fx))
_rng = np.random.default_rng(20260928)
_pos = fal.variant.str.extract(r"(\d+)")[0].astype(int).values
_g = [np.where(_pos == u)[0] for u in np.unique(_pos)]
_bs = []
for _ in range(20000):
    _s = np.concatenate([_g[i] for i in _rng.integers(0, len(_g), len(_g))])
    _bs.append(np.sum(fx[_s] * fy[_s]) / np.sum(fx[_s] ** 2))
b_fal_ci = (float(np.quantile(_bs, 0.025)), float(np.quantile(_bs, 0.975)))

# Uncorrected breaches, for the descriptive sentence in the text.
n_below = int((out.het < FLOOR).sum())
n_nominal = int((2 * stats.norm.sf((FLOOR - out.het) / out.het_se) < 0.05).sum())

# p-to-rescue: the mutant fraction p that would place each intersection
# variant exactly at its own measured level, holding total abundance fixed.
# floor(p) = WW * (1-p)^2  =>  p = 1 - sqrt(het / WW)
out["p_to_rescue"] = 1 - np.sqrt(out.het / WW)
sub = out[out.in_intersection]
p_resc = sub.p_to_rescue
# also over the fixed-floor set, the most permissive primary specification
ff = out[out.variant.isin(sets["fixed floor, two-sided normal"])]
p_resc_ff = ff.p_to_rescue

# --------------------------------------------------------------------------
# Subunit-fraction sensitivity (revision of 8 November 2026).
#
# p = 0.5 is an assumption, not a measurement.  Every specification is
# recomputed with the dimer floor WW*(1-p)^2 and its SE WW_SE*(1-p)^2 over a
# grid of assumed mutant fractions.  Alternative assembly orders and the
# exactly additive 200% reference are recomputed on the same grid.
# --------------------------------------------------------------------------
P_GRID = (0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60)
ALPHA = 0.05
Z_TWO = float(stats.norm.isf(ALPHA / (2 * FAMILY)))    # Bonferroni two-sided
Z_ONE = float(stats.norm.isf(ALPHA / FAMILY))          # Bonferroni one-sided


def floor_at(p, n_subunits=2, reference=None):
    """Conditional floor and SE at assumed mutant fraction p."""
    ref, ref_se = (WW, WW_SE) if reference is None else (reference, 0.0)
    w = (1.0 - p) ** n_subunits
    return ref * w, ref_se * w


def primary_pvalues(df, floor, fse, fn=FLOOR_N):
    """Bonferroni-adjusted probabilities of the five primary specifications."""
    return {
        "fixed floor, two-sided normal": spec_fixed_floor(df, floor)[1],
        "floor uncertainty propagated, two-sided normal": spec_propagated_two(df, floor, fse)[1],
        "floor uncertainty propagated, one-sided normal": spec_propagated_one(df, floor, fse)[1],
        "Welch-Satterthwaite, reported cell counts": spec_welch(df, floor, fse, fn)[1],
        "log-ratio delta method, two-sided normal": spec_log_ratio(df, floor, fse)[1],
    }


# Control-based normalization scenarios. Mixed-arm current is normalized to
# the SB-only WT line (100%), not to the LP-only line used for the alone arm.
# Under linear dose scaling the extra LP WT contributes WW - WT, giving
# p_LP = (WW-WT)/WW and a surviving WW current WT**2 / WW. The inferred
# fraction depends on the same control, so propagate the derivative of that
# expression rather than treating the fraction as an independent estimate.
SOURCE_P = 1.0 - WT / WW
SOURCE_FLOOR = WT ** 2 / WW
SOURCE_FLOOR_SE = WT ** 2 * WW_SE / WW ** 2
POWER_EXPONENT = np.log(WW / WT) / np.log(2.0)
INERT = out.set_index("variant").loc["W822X"]
INERT_CALIBRATION = control_floor(WW, WW_SE, float(INERT.het), float(INERT.het_se),
                                  int(WW_N), int(INERT.het_n))
SOURCE_SCENARIOS = {
    "SB-normalized linear-dose interpretation": (SOURCE_P, SOURCE_FLOOR, SOURCE_FLOOR_SE),
    "equal-dose additive calibration sensitivity": (0.5, STRICT, 0.0),
    "inert p.Trp822Ter calibration, propagated control SE": (
        INERT_CALIBRATION["p"], INERT_CALIBRATION["floor"], INERT_CALIBRATION["floor_se"]),
    "inert p.Trp822Ter calibration, fixed floor": (
        INERT_CALIBRATION["p"], INERT_CALIBRATION["floor"], 0.),
}
source_scenarios = {}
source_sets = []
for name, (p_source, fl, fse) in SOURCE_SCENARIOS.items():
    floor_df = INERT_CALIBRATION["floor_df"] if name.startswith("inert") else WW_N - 1
    pv = primary_pvalues(out, fl, fse, floor_df + 1)
    hit_sets = {lab: sorted(out.loc[np.asarray(padj) < ALPHA, "variant"])
                for lab, padj in pv.items()}
    intersection = sorted(set.intersection(*(set(v) for v in hit_sets.values())))
    source_sets.append(set(intersection))
    source_scenarios[name] = {
        "assumed_mutant_fraction": p_source, "floor": fl, "floor_se": fse,
        "floor_effective_df": floor_df,
        "counts": {lab: len(v) for lab, v in hit_sets.items()},
        "sets": hit_sets, "intersection": intersection,
        "adjusted_p_values": {lab: {v: float(q) for v, q in zip(out.variant, padj)}
                              for lab, padj in pv.items()},
    }
SOURCE_INTERSECTION = sorted(set.intersection(*source_sets))
out["in_source_informed_intersection"] = out.variant.isin(SOURCE_INTERSECTION)

# Every reported floor has its own robustness table and joint uncertainty
# surface. Historical p=.5 numbers remain labeled as such rather than being
# reused for the smaller 45.79% floor.
reported_floor_scenarios = dict(source_scenarios)
reported_floor_scenarios["equal fraction with observed WT+WT reference"] = {
    "assumed_mutant_fraction": .5, "floor": FLOOR, "floor_se": FLOOR_SE,
    "floor_effective_df": WW_N - 1, **test_floor(out, FLOOR, FLOOR_SE)}
robustness, joint_surface = robustness_rows(out, reported_floor_scenarios)
higher_floor_robustness = robustness_counts(
    out, INERT_CALIBRATION["floor"], INERT_CALIBRATION["floor_se"])
for row in higher_floor_robustness["effective_df_scenarios"]:
    row["log_ratio_hgvs"] = [scn5a_hgvs(v) for v in row["log_ratio_variants"]]
for row in higher_floor_robustness["discrepancy_scenarios"]:
    for key in ("propagated_normal", "log_ratio"):
        row[key + "_hgvs"] = [scn5a_hgvs(v) for v in row[key + "_variants"]]
for row in (*robustness, *joint_surface):
    row["hgvs_p_NP_000326_2"] = scn5a_hgvs(row["variant"])
same_two_floor_intervals = {
    "observed_WT_WT_relative_SE": floor_interval_for_intersection(
        out, ["D356N", "G1712C"], WW_SE / WW),
    "fixed_floor": floor_interval_for_intersection(out, ["D356N", "G1712C"], 0.),
}
pd.DataFrame(robustness).to_csv("scn5a_robustness.csv", index=False)
pd.DataFrame(joint_surface).to_csv("scn5a_joint_uncertainty.csv", index=False)
json.dump({"family_size": FAMILY, "scenarios": reported_floor_scenarios,
           "inert_control": {"variant": "W822X", "mean": float(INERT.het),
                             "se": float(INERT.het_se), "cells": int(INERT.het_n),
                             **INERT_CALIBRATION},
           "same_two_floor_intervals": same_two_floor_intervals,
           "replication": {"unit": "cells treated as independent",
                           "variant_cell_range": [int(out.het_n.min()), int(out.het_n.max())],
                           "control_cells": int(WW_N),
                           "transfections_plates_days": "not reported; effective df and SE inflation are scenarios"},
           "power_law_equivalence": {"identity": "100*0.5**log2(W/100)=100**2/W",
                                     "independent_corroboration": False},
           "joint_surface_definition": "Largest common multiplier of all input SEs that retains rejection is z/tcrit(df); equality is the P=.05 boundary. Values <1 mean reported SEs already fail at that effective df. Variant-only inflation is also provided.",
           "robustness": robustness,
           "higher_floor_robustness": higher_floor_robustness,
           "joint_uncertainty": joint_surface},
          open("scn5a_floor_scenarios.json", "w"), indent=2, allow_nan=False)

fraction_rows = []
fraction_sets = {}
for p_ in P_GRID:
    fl, fse_ = floor_at(p_)
    pv = primary_pvalues(out, fl, fse_)
    row = {"p": p_, "dimer_floor": fl, "dimer_floor_se": fse_}
    psets = {}
    for label, padj in pv.items():
        psets[label] = set(out.loc[np.asarray(padj) < ALPHA, "variant"])
        row[label] = len(psets[label])
    inter = sorted(set.intersection(*psets.values()))
    row["intersection"] = len(inter)
    for n_sub, name in ((3, "trimer"), (4, "tetramer")):
        fl_n, _ = floor_at(p_, n_sub)
        row[f"{name}, fixed floor"] = int((spec_fixed_floor(out, fl_n)[1] < ALPHA).sum())
    fl_add, _ = floor_at(p_, 2, reference=2 * WT)
    row["additive 200% reference, dimer fixed floor"] = int(
        (spec_fixed_floor(out, fl_add)[1] < ALPHA).sum())
    fraction_rows.append(row)
    fraction_sets[f"{p_:.2f}"] = {"intersection": inter,
                                  "log_ratio": sorted(psets["log-ratio delta method, two-sided normal"]),
                                  "intersection_equals_log_ratio": inter == sorted(
                                      psets["log-ratio delta method, two-sided normal"])}
fraction_table = pd.DataFrame(fraction_rows)

# Nesting of the five specifications at p = 0.5.
_order = ["log-ratio delta method, two-sided normal",
          "Welch-Satterthwaite, reported cell counts",
          "floor uncertainty propagated, two-sided normal",
          "floor uncertainty propagated, one-sided normal",
          "fixed floor, two-sided normal"]
nesting = {f"{a} subset of {b}": sets[a] <= sets[b] for a, b in zip(_order, _order[1:])}
nesting["one-sided propagated equals fixed floor"] = (
    sets["floor uncertainty propagated, one-sided normal"] == sets["fixed floor, two-sided normal"])
nesting["intersection equals log-ratio set"] = (
    set(INTERSECTION) == sets["log-ratio delta method, two-sided normal"])


def _significance_fraction(row, label):
    """Largest assumed fraction at which the variant remains significant."""
    one = out.loc[[row.name]]

    def h(p_):
        fl, fse_ = floor_at(p_)
        padj = float(np.asarray(primary_pvalues(one, fl, fse_)[label])[0])
        return np.log(max(padj, 1e-300)) - np.log(ALPHA)

    from scipy.optimize import brentq
    lo, hi = 0.01, 0.99
    if h(lo) >= 0 or h(hi) <= 0:
        return float("nan")
    return float(brentq(h, lo, hi, xtol=1e-10))


def _tau_star(row):
    """Discrepancy SD (percentage points, quadrature) at which the propagated
    two-sided test stops rejecting at p = 0.5."""
    gap = FLOOR - row.het
    if gap <= 0:
        return float("nan")
    need = (gap / Z_TWO) ** 2 - row.het_se ** 2 - FLOOR_SE ** 2
    return float(np.sqrt(need)) if need > 0 else float("nan")


def _tau_star_log(row):
    """Discrepancy SD (percentage points at the floor) at which the log-ratio
    test stops rejecting at p = 0.5; tau enters the log scale as tau/FLOOR."""
    if row.het >= FLOOR:
        return float("nan")
    need = ((np.log(FLOOR / row.het) / Z_TWO) ** 2
            - (row.het_se / row.het) ** 2 - (FLOOR_SE / FLOOR) ** 2)
    return float(FLOOR * np.sqrt(need)) if need > 0 else float("nan")


def _nu_star(stat):
    """Smallest effective degrees of freedom at which a statistic still exceeds
    the Bonferroni two-sided t critical value."""
    from scipy.optimize import brentq
    if not np.isfinite(stat) or stat <= Z_TWO:
        return float("nan")
    f = lambda nu: float(stats.t.isf(ALPHA / (2 * FAMILY), nu)) - stat
    return float(brentq(f, 1.0 + 1e-9, 1e8, xtol=1e-8))


tipping = {}
for idx, r in out[out.variant.isin(sets["fixed floor, two-sided normal"])].iterrows():
    tipping[r.variant] = {
        "het": float(r.het), "het_se": float(r.het_se),
        "p_point_estimate_rescue": float(r.p_to_rescue),
        "p_significance_rescue": {lab: _significance_fraction(r, lab) for lab in PRIMARY},
        "tau_star_propagated_two_sided": _tau_star(r),
        "tau_star_log_ratio": _tau_star_log(r),
        "nu_star_log_ratio": _nu_star(float(r.z_log)),
        "nu_star_propagated_two_sided": _nu_star(float(r.z_prop)),
        "in_intersection": bool(r.in_intersection),
    }
t_crit_examples = {str(nu): float(stats.t.isf(ALPHA / (2 * FAMILY), nu)) for nu in (9, 19, 29, 49)}
t_crit_examples["normal"] = Z_TWO

# --------------------------------------------------------------------------
# Report.
# --------------------------------------------------------------------------
def emit(fh):
    w = fh.write
    w("SCN5A specification sensitivity\n")
    w(f"generated {date.today().isoformat()} by run_scn5a.py from {CSV}\n")
    w(f"python {sys.version.split()[0]}  numpy {np.__version__}  "
      f"pandas {pd.__version__}  scipy {stats.__name__.split('.')[0]} "
      f"{__import__('scipy').__version__}\n")
    w(f"{SEED_INFO}\n\n")
    w(f"controls read from file: WT = {WT:.1f}%,  "
      f"WT+WT = {WW:.1f} +- {WW_SE:.1f}% (n = {WW_N})\n")
    w(f"dimer floor = WT+WT / 4 = {FLOOR:.4g} +- {FLOOR_SE:.4g}%  "
      f"(conditional on equal abundance, random assembly,\n"
      f"  unchanged total abundance, unchanged WW activity)\n")
    w(f"family size for Bonferroni = {FAMILY}\n\n")
    w(f"{n_below}/{FAMILY} variants sit numerically below the floor; "
      f"{n_nominal} of those are nominally significant before correction.\n\n")
    w("specification counts after Bonferroni correction (alpha = 0.05):\n")
    for k, v in counts.items():
        w(f"  {v:>3d}   {k}\n")
    w(f"\nintersection of the {len(PRIMARY)} primary specifications "
      f"({len(INTERSECTION)}): {', '.join(INTERSECTION)}\n")
    w(f"\nWT+WT = {WW:.1f} is {supra:.2f} SE above the {2*WT:.0f} that strict "
      f"additivity predicts.\n"
      f"  floor anchored to the observed control: {FLOOR:.2f}%  -> "
      f"{counts['fixed floor, two-sided normal']} survivors\n"
      f"  floor under strict additivity:          {STRICT:.2f}%  -> "
      f"{len(strict_hits)} survivors  ({', '.join(strict_hits)})\n")
    w(f"\nFalconer mapping (x = WT+WT - 2*alone, y = WT+WT - mixed):\n"
      f"  b = {b_fal:.4f}  95% CI [{b_fal_ci[0]:.4f}, {b_fal_ci[1]:.4f}]"
      f"   d/a = {1-2*b_fal:+.4f}\n"
      "  conditional on commensurate LP/SB references and linear mutant dose scaling; not a fitness coefficient\n")
    w("\nillustrative additive-model residuals for "
      f"({', '.join(CONTROLS)}):\n  residuals "
      + ", ".join(f"{r:+.1f}" for r in resid)
      + f" (RMS {rms_model_error:.1f}) against a floor SE of {FLOOR_SE:.2f}\n")
    w("significant variants under hypothetical discrepancy SDs added in quadrature:\n")
    for k, v in model_error_rows.items():
        w(f"  +{k:2d} points: {len(v):2d}  {', '.join(v) if len(v) <= 6 else ''}\n")
    w("\ndescriptive mixed-below-alone contrast (unchanged mutant contribution\n"
      "and independent errors would be needed for mechanistic interpretation):\n")
    if len(mf):
        for _, r in mf.iterrows():
            w(f"  {r.variant}: alone {r.homo:.1f}+/-{r.homo_se:.1f}, "
              f"mixed {r.het:.1f}+/-{r.het_se:.1f}, difference "
              f"{r['diff']:+.1f}+/-{r.diff_se:.1f}, z = {r.z:.2f}, "
              f"p = {r.p:.3f} (uncorrected)\n")
        w("  the numerical direction is descriptive; the contrast is not significant\n")
    else:
        w("  none\n")
    w("\np that would place each variant exactly at its own measured level,\n"
      "holding total complex abundance fixed (p = 1 - sqrt(het/WT+WT)).\n"
      "This is the mutant fraction an abundance imbalance would have to reach\n"
      "to rescue the random-assembly model for that variant:\n")
    w(f"  intersection ({len(sub)}): " + ", ".join(
        f"{v} {pv:.3f}" for v, pv in zip(sub.variant, p_resc)) + "\n")
    w(f"    range {p_resc.min():.3f} to {p_resc.max():.3f}, "
      f"median {np.median(p_resc):.3f}\n")
    w(f"  fixed-floor set ({len(ff)}): range {p_resc_ff.min():.3f} to "
      f"{p_resc_ff.max():.3f}, median {np.median(p_resc_ff):.3f}\n")
    w("\nsubunit-fraction sensitivity (dimer floor WT+WT*(1-p)^2; counts after "
      "Bonferroni):\n")
    w(fraction_table.to_string(index=False, float_format=lambda v: f"{v:.3f}") + "\n")
    for k, v in fraction_sets.items():
        w(f"  p = {k}: intersection {', '.join(v['intersection']) or 'none'}"
          f" (equals log-ratio set: {v['intersection_equals_log_ratio']})\n")
    w("\nnesting at p = 0.5:\n")
    for k, v in nesting.items():
        w(f"  {k}: {v}\n")
    w("\ncontrol-based normalization scenarios (cells treated as independent):\n")
    w(f"  WT+WT/additive expectation = {WW / (2 * WT):.4f}; "
      "this discrepancy does not identify an unknown nonlinear response curve.\n")
    for name, result in source_scenarios.items():
        w(f"  {name}: p = {result['assumed_mutant_fraction']:.6f}, "
          f"floor = {result['floor']:.4f} +/- {result['floor_se']:.4f}; "
          f"intersection {result['intersection']}\n")
        for label, count in result["counts"].items():
            w(f"    {count}: {label}\n")
    w("\nper-variant tipping points (fixed-floor set at p = 0.5):\n"
      "  p_sig = largest assumed fraction still significant (log-ratio test);\n"
      "  tau* = discrepancy SD removing the propagated two-sided rejection;\n"
      "  nu* = smallest effective df keeping the log-ratio rejection\n")
    for v, tp in tipping.items():
        w(f"  {v:7s} point-rescue p {tp['p_point_estimate_rescue']:.3f}  "
          f"p_sig {tp['p_significance_rescue']['log-ratio delta method, two-sided normal']:.3f}  "
          f"tau* {tp['tau_star_propagated_two_sided']:.2f} (log {tp['tau_star_log_ratio']:.2f})  "
          f"nu* {tp['nu_star_log_ratio']:.1f}\n")
    w("  Bonferroni two-sided critical values: " + ", ".join(
        f"df {k}: {v:.2f}" for k, v in t_crit_examples.items()) + "\n")

emit(sys.stdout)
with open("scn5a_specifications.log", "w") as fh:
    emit(fh)

out.to_csv("scn5a_specifications.csv", index=False)
json.dump({
    "generated": date.today().isoformat(),
    "input": CSV,
    "controls": {"WT": WT, "WT_WT": WW, "WT_WT_se": WW_SE, "WT_WT_n": WW_N},
    "floor": {"value": FLOOR, "se": FLOOR_SE, "n": FLOOR_N, "n_subunits": 2},
    "family_size": FAMILY,
    "n_below_floor": n_below,
    "n_nominal_uncorrected": n_nominal,
    "counts": counts,
    "sets": {k: sorted(v) for k, v in sets.items()},
    "intersection": INTERSECTION,
    "intersection_hgvs_NP_000326_2": [scn5a_hgvs(v) for v in INTERSECTION],
    "direct_arm_comparison_scope": "Alone and mixed percentages use LP-only and SB-only references; subtraction requires an assumed reference conversion and is not an uncalibrated current contrast.",
    "source_normalization": {
        "mixed_reference": "SB-only WT cell line; the alone arm uses LP-only WT",
        "source": "O'Neill et al. 2022, Data Analysis; doi:10.1016/j.gim.2022.02.010",
        "WT_WT_over_additive_expectation": WW / (2 * WT),
        "linear_dose_fraction": SOURCE_P,
        "equal_dose_power_law_exponent": float(POWER_EXPONENT),
        "power_law_equivalence": "100*0.5**log2(W/100)=100**2/W; identical floor by construction, not an independent scenario or corroboration",
        "inert_control_calibration": INERT_CALIBRATION,
        "reported_floor_scenarios_file": "scn5a_floor_scenarios.json",
        "replication_unit": "cells treated as independent; 27–164 per mixed variant; transfections, plates and days not reported",
        "scenarios": source_scenarios,
        "intersection_across_scenarios": SOURCE_INTERSECTION,
        "intersection_hgvs_NP_000326_2": [scn5a_hgvs(v) for v in SOURCE_INTERSECTION],
        "limits": ["The inferred fraction assumes linear current versus assembly-competent abundance and unchanged WT function.",
                   "The equal-dose 50% floor is an additive calibration sensitivity; an unknown nonlinear response function does not identify a floor.",
                   "The power-law floor equals the linear control floor algebraically and supplies no independent corroboration.",
                   "Treating p.Trp822Ter as inert uses its134.2% mixed response as a WT contribution proxy; the source reports4.7% current alone and suggests a poorly conducting multimer. Its biological validity is unresolved.",
                   "The normalized SB reference is treated as fixed; shared-reference/batch covariance is unavailable."],
    },
    "strict_additivity": {"floor": STRICT, "survivors": strict_hits,
                          "supra_additive_SE": supra},
    # Falconer places the mutant homozygote at -a and wild type at +a. With
    # x = phi_WT - phi_hom = 2a and y = phi_WT - phi_het = a - d, the ratio is
    # y/x = 1/2 - d/(2a), so d/a = 1 - 2b. An earlier version wrote 2b - 1,
    # which reports the sign of the deviation backwards; the magnitude and the
    # position relative to the sign-invariant additive point of 1/2 are
    # unaffected.
    "falconer": {"b": b_fal, "ci": list(b_fal_ci), "d_over_a": 1 - 2 * b_fal,
                 "convention": "x = phi_WT - phi_hom = 2a; y = phi_WT - phi_het = a - d"},
    "illustrative_additive_residuals": resid,
    "illustrative_additive_residual_rms": rms_model_error,
    "discrepancy_interpretation": "analyst-chosen Gaussian SDs; not an estimated floor-error distribution",
    "survivors_by_hypothetical_discrepancy_sd": {str(k): v for k, v in model_error_rows.items()},
    "descriptive_mixed_below_alone": (
        mf[["variant", "homo", "het", "z", "p"]].to_dict("records")),
    "p_to_rescue_intersection": {v: float(p) for v, p in zip(sub.variant, p_resc)},
    "p_to_rescue_fixed_floor": {v: float(p) for v, p in zip(ff.variant, p_resc_ff)},
    "fraction_sensitivity": {
        "interpretation": ("assumed mutant fraction in the assembly-competent pool; "
                           "analyst-chosen grid, not an estimate of p"),
        "table": fraction_table.to_dict("records"),
        "sets": fraction_sets,
    },
    "nesting_at_p_0_5": nesting,
    "tipping_points": tipping,
    "tipping_definitions": {
        "p_point_estimate_rescue": "fraction at which the floor equals the measured mixed current",
        "p_significance_rescue": "largest assumed fraction at which the Bonferroni-adjusted test still rejects",
        "tau_star_propagated_two_sided": "Gaussian discrepancy SD (percentage points) at which the propagated two-sided test stops rejecting, p = 0.5",
        "tau_star_log_ratio": "discrepancy SD (percentage points at the floor, entering the log scale as tau/floor) at which the log-ratio test stops rejecting, p = 0.5",
        "nu_star_log_ratio": "smallest effective degrees of freedom at which the log-ratio statistic still exceeds the Bonferroni two-sided t critical value, p = 0.5",
    },
    "bonferroni_t_critical_values": t_crit_examples,
    "versions": {"python": sys.version.split()[0], "numpy": np.__version__,
                 "pandas": pd.__version__, "scipy": __import__("scipy").__version__},
}, open("scn5a_specifications.json", "w"), indent=2)

# --------------------------------------------------------------------------
# Figure 1 (submitted as Fig. 3), drawn at its printed width.
# --------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from palette import (BLUE, BLUE_DARK, BLUE_LIGHT, BLUE_MID, CONTEXT, FIG_WIDTH, GOLD, GOLD_DARK,
                     GOLD_LIGHT, INK, MUTED, TEAL, marker_handle, panel_label, rich_text,
                     rounded_bar, soft_grid, tint, use_figure_style)

use_figure_style()
NUMBER_WORDS = {4: "four", 5: "five", 6: "six"}
FIG_H = 6.5
fig = plt.figure(figsize=(FIG_WIDTH, FIG_H))
# Fixed axes positions (figure fractions): panel B needs a label column on its
# left for the full specification names.
ax1 = fig.add_axes([0.085, 0.625, 0.90, 0.335])
ax2 = fig.add_axes([0.27, 0.11, 0.24, 0.325])
ax3 = fig.add_axes([0.655, 0.11, 0.32, 0.325])
P_HALF = [("p", True), (" = 0.5", False)]            # italic p in panel titles

# ---- A: ranked mixed currents against the floor ---------------------------
s = out.sort_values("het").reset_index(drop=True)
xs = np.arange(len(s))
core = s.variant.isin(source_scenarios["SB-normalized linear-dose interpretation"]["intersection"]).values
inert_extra = s.variant.isin(source_scenarios["inert p.Trp822Ter calibration, propagated control SE"]["intersection"]).values & ~core
other = ~core & ~inert_extra
ax1.axhspan(SOURCE_FLOOR - SOURCE_FLOOR_SE, SOURCE_FLOOR + SOURCE_FLOOR_SE, color=BLUE_LIGHT, alpha=0.55, lw=0, zorder=0)
ax1.axhline(SOURCE_FLOOR, color=BLUE, lw=1.3, zorder=1)
ax1.axhspan(INERT_CALIBRATION["floor"] - INERT_CALIBRATION["floor_se"],
            INERT_CALIBRATION["floor"] + INERT_CALIBRATION["floor_se"],
            color=tint(TEAL, 0.12), lw=0, zorder=0)
ax1.axhline(INERT_CALIBRATION["floor"], color=TEAL, lw=1.1, ls="--", zorder=1)
ax1.axhline(WT, color=INK, lw=0.7, ls=":", zorder=1)
ax1.axhline(0.75 * WT, color=TEAL, lw=0.9, ls="-.", zorder=1)
ax1.errorbar(xs[other], s.het[other], yerr=s.het_se[other], fmt="o", ms=3.8,
             color=CONTEXT, mec="white", mew=0.4, lw=0.7, zorder=3)
ax1.errorbar(xs[inert_extra], s.het[inert_extra], yerr=s.het_se[inert_extra], fmt="o", ms=4.6,
             color=TEAL, mec="white", mew=0.4, lw=0.8, zorder=4)
# The flagged variants have standard errors smaller than their markers, so
# their error bars are drawn over the markers in a darker gold.
ax1.plot(xs[core], s.het[core], "o", ms=6.2, color=GOLD, mec="white", mew=0.6, zorder=5)
ax1.errorbar(xs[core], s.het[core], yerr=s.het_se[core], fmt="none", ecolor=GOLD_DARK,
             elinewidth=1.0, zorder=6)
ax1.set_xlim(-1, len(s))
ax1.set_ylim(0, 250)
ax1.set_ylabel("peak current density,\n% of single wild type")
panel_label(ax1, "A", "wild type + variant against the conditional floor")
handles = [marker_handle("o", GOLD, 6.2, mec="white", mew=0.6),
           marker_handle("o", TEAL, 4.6, mec="white", mew=0.4),
           marker_handle("o", CONTEXT, 3.8, mec="white", mew=0.4)]
labels = ["significantly below both floors",
          "significantly below the higher floor only",
          "other variants"]
ax1.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.02, 1.0),
           handlelength=1.2, labelspacing=0.4, borderaxespad=0.2)

# Each reference line is named at its right end, under the line. The check
# below confirms that no variant or error bar lies under a label.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
to_data = ax1.transData.inverted()
for level, text, offset, va in ((WT, "single wild type, 100%", 4, "bottom"),
                    (0.75 * WT, "O'Neill et al. threshold, 75%", -4, "top"),
                    (SOURCE_FLOOR, f"WT+WT-derived calibration floor, {SOURCE_FLOOR:.1f}%", -4, "top")):
    anchor_x, dx, ha = (.02, 2, "left") if level == WT else (1, -2, "right")
    note = ax1.annotate(text, (anchor_x, level), xycoords=("axes fraction", "data"), xytext=(dx, offset),
                        textcoords="offset points", ha=ha, va=va, fontsize=8, color=INK)
    box = note.get_window_extent(renderer)
    (x_lo, y_lo), (x_hi, y_hi) = to_data.transform([(box.x0, box.y0), (box.x1, box.y1)])
    under = ((xs > x_lo - 0.6) & (xs < x_hi + 0.6)
             & (s.het + s.het_se > y_lo) & (s.het - s.het_se < y_hi))
    if under.any():
        raise ValueError(f"the label '{text}' would cover a variant; move it")

ax1.annotate(f"inert p.Trp822Ter calibration floor, {INERT_CALIBRATION['floor']:.1f}%",
             (0., INERT_CALIBRATION["floor"]), xycoords=("axes fraction", "data"),
             xytext=(.02, 140), textcoords=("axes fraction", "data"),
             ha="left", va="bottom", fontsize=8, color=TEAL,
             arrowprops={"arrowstyle": "-", "color": TEAL, "lw": .65, "relpos": (0, 0)})

# Variant names sit below the axis on stepped rows, each joined to its own
# position by a short vertical leader, so no line crosses the plotting area.
# Labels are placed from right to left and left aligned at their leader; a
# label drops one row whenever it would run into the label to its right.
blend = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
axis_height_pt = ax1.get_window_extent().height * 72 / fig.dpi
row_pt, first_pt, gap = 10.5, 8.5, 0.8
row_ends = []                           # leftmost occupied x in each row
placed = []
for i in sorted(xs[core], reverse=True):
    name = scn5a_hgvs(s.variant[i])
    probe = ax1.text(0, 0, name, fontsize=8)
    width_px = probe.get_window_extent(renderer).width
    probe.remove()
    x0 = to_data.transform((0, 0))[0]
    width = to_data.transform((width_px, 0))[0] - x0 + 0.55
    row = next((r for r, left in enumerate(row_ends) if i + width + gap <= left), len(row_ends))
    if row == len(row_ends):
        row_ends.append(i)
    else:
        row_ends[row] = i
    placed.append((i, name, row))
for i, name, row in placed:
    depth = (first_pt + row * row_pt) / axis_height_pt
    ax1.plot([i, i], [0, -depth], transform=blend, color=GOLD, lw=0.8, clip_on=False,
             solid_capstyle="butt")
    ax1.text(i + 0.45, -depth, name, transform=blend, fontsize=8, color=INK,
             ha="left", va="center")
ax1.set_xticks([])
label_y = -(first_pt + len(row_ends) * row_pt + 6) / axis_height_pt
rich_text(ax1, 0.5, label_y, [("51 ", False), ("SCN5A", True),
                              (" variants (50 missense and p.Trp822Ter), ranked by mixed state current", False)],
          ha="center", va="top", fontsize=8)

# ---- B: the control assumptions determine the intersection count --------
scenario_order = ["SB-normalized linear-dose interpretation",
                  "equal-dose additive calibration sensitivity",
                  "equal fraction with observed WT+WT reference",
                  "inert p.Trp822Ter calibration, propagated control SE",
                  "inert p.Trp822Ter calibration, fixed floor"]
names_b = ["WT+WT: 45.8 ± 1.6", "additive: 50, fixed",
           "equal fractions: 54.6 ± 1.9", "inert p.Trp822Ter: 82.5 ± 7.0",
           "inert p.Trp822Ter: 82.5, fixed"]
vals = [len(reported_floor_scenarios[name]["intersection"]) for name in scenario_order]
ax2.set_ylim(len(vals)-.3, -.7)
ax2.set_xlim(0, 12)
ax2.set_xticks(range(0,13,2))
soft_grid(ax2, "x")
for yy, value, colour in zip(range(len(vals)), vals, [BLUE, BLUE_MID, CONTEXT, TEAL, TEAL]):
    rounded_bar(ax2, 0, value, yy, .60, radius=2.5, color=colour, zorder=2)
    ax2.text(value+.2, yy, str(value), va="center", ha="left", fontsize=8, color=INK)
ax2.set_yticks(range(len(vals)), names_b)
ax2.tick_params(axis="y", length=0)
ax2.set_xlabel("five-test intersection (Bonferroni across 51)")
panel_label(ax2, "B", "control and floor assumptions",
            x=0.012, y=ax2.get_position().y1, transform=fig.transFigure)

# ---- C: counts across assumed fractions -----------------------------------
fine = np.round(np.arange(0.38, 0.6201, 0.0025), 6)
curves = {lab: [] for lab in PRIMARY}
inter_curve = []
for p_ in fine:
    fl, fse_ = floor_at(p_)
    pv = primary_pvalues(out, fl, fse_)
    cur = {}
    for lab in PRIMARY:
        cur[lab] = set(out.loc[np.asarray(pv[lab]) < ALPHA, "variant"])
        curves[lab].append(len(cur[lab]))
    inter_curve.append(len(set.intersection(*cur.values())))
# The log ratio set equals the five specification intersection at every grid
# fraction, so one gold line shows both, over a light fill; the other four
# specifications take the blue ramp from the most permissive (light) to the
# strictest (dark).
names = dict(zip(PRIMARY, ["fixed, two sided", "propagated, two sided", "propagated, one sided",
                           "Welch", "log ratio"]))
log_key = list(PRIMARY)[4]
if curves[log_key] != inter_curve:
    raise ValueError("log ratio counts no longer equal the intersection; restore a separate line")
top = max(max(v) for v in curves.values()) + 1.5
soft_grid(ax3, "y")
ax3.fill_between(fine, inter_curve, step="mid", color=tint(GOLD_LIGHT, 0.32), lw=0, zorder=1)
ramp = {list(PRIMARY)[2]: (BLUE_LIGHT, "-", 1.5), list(PRIMARY)[0]: (BLUE_MID, "--", 1.0),
        list(PRIMARY)[1]: (BLUE, "-", 1.0), list(PRIMARY)[3]: (BLUE_DARK, "-", 1.0)}
for lab in [list(PRIMARY)[i] for i in (2, 0, 1, 3)]:
    c, ls, lw = ramp[lab]
    ax3.step(fine, curves[lab], where="mid", color=c, ls=ls, lw=lw, label=names[lab], zorder=3)
ax3.step(fine, inter_curve, where="mid", color=GOLD, lw=2.0, zorder=4,
         label="log ratio = intersection")
# The p = 0.5 marker stops below the legend so that the two never overlap.
ax3.plot([0.5, 0.5], [0, 13.2], color=MUTED, lw=0.7, ls=":", zorder=2)
ax3.plot([SOURCE_P, SOURCE_P], [0, 4], color=TEAL, lw=0.9, ls=":", zorder=3)
ax3.annotate("WT+WT $p$ = 0.542", (SOURCE_P, 4), xytext=(3, 2),
             textcoords="offset points", fontsize=7, color=TEAL, rotation=90,
             va="bottom", ha="left")
# Intersection counts at four fractions: beside the line where the blue curves
# leave room, and at the right end of the line for the last step.
for p_mark, off, align in ((0.40, (0, 5), "center"), (0.50, (-6, 5), "center"),
                           (0.52, (3, -9.5), "center")):
    k = int(np.argmin(np.abs(fine - p_mark)))
    ax3.annotate(str(inter_curve[k]), (fine[k], inter_curve[k]), xytext=off,
                 textcoords="offset points", ha=align, fontsize=8, color=INK)
k_end = int(np.argmin(np.abs(fine - 0.60)))
if inter_curve[k_end] != inter_curve[-1]:
    raise ValueError("the intersection changes after p = 0.60; label that fraction on the line")
ax3.annotate(str(inter_curve[-1]), (fine[-1], inter_curve[-1]), xytext=(4, 0),
             textcoords="offset points", ha="left", va="center", fontsize=8, color=INK,
             annotation_clip=False)
ax3.plot(INERT_CALIBRATION["p"], 10, "s", color=TEAL, mec="white", mew=.5, ms=5, zorder=6)
ax3.annotate("inert control: 10\nwith control SE", (INERT_CALIBRATION["p"], 10),
             xytext=(.405, 3.1), textcoords="data", fontsize=7, color=TEAL,
             arrowprops={"arrowstyle": "-", "color": TEAL, "lw": .7}, ha="left", va="center")
ax3.set_xlabel("assumed mutant fraction $p$\nin the assembly competent pool")
ax3.set_ylabel("variants significantly below\nthe dimer floor (Bonferroni)")
panel_label(ax3, "C", "sensitivity to the assumed fraction")
ax3.legend(loc="upper right", handlelength=1.8, labelspacing=0.3, borderaxespad=0.2)
ax3.set_xlim(fine[0], fine[-1])
ax3.set_ylim(0, top)

fig.savefig("Figure1.png", dpi=300)
fig.savefig("Figure1.pdf")

# A reviewable view of the joint uncertainty surface for the two variants
# defining the lower-floor intersection. Complete51-variant surfaces are CSV.
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
joint_frame = pd.DataFrame(joint_surface)
joint_dfs = ["5", "7", "9", "19", "49", "199", "normal"]
surface_rows, surface_labels = [], []
for scenario, short in zip(scenario_order,
                           ["45.79 ± 1.61", "50, fixed", "54.6 ± 1.925",
                            "82.46 ± 7.02", "82.46, fixed"]):
    for variant in ("D356N", "G1712C"):
        values = joint_frame[(joint_frame.scenario == scenario) &
                             (joint_frame.variant == variant)].set_index("effective_df_label")
        surface_rows.append(values.loc[joint_dfs, "max_all_se_inflation_log_ratio"].to_numpy())
        surface_labels.append(f"{short}: {scn5a_hgvs(variant)}")
surface_fig, surface_ax = plt.subplots(figsize=(FIG_WIDTH, 4.5))
colors = LinearSegmentedColormap.from_list("inflation", [GOLD_LIGHT, "white", BLUE_DARK])
matrix = np.array(surface_rows)
surface_ax.imshow(matrix, cmap=colors, norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=2.5), aspect="auto")
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        value = matrix[i, j]
        surface_ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if value > 1.9 else INK)
for edge in np.arange(1.5, len(surface_rows)-1, 2):
    surface_ax.axhline(edge, color="white", lw=2)
surface_ax.set_yticks(range(len(surface_labels)), surface_labels)
surface_ax.set_xticks(range(len(joint_dfs)), joint_dfs)
surface_ax.set_xlabel("assumed effective degrees of freedom")
surface_ax.set_ylabel("floor (% ± SE) and variant")
surface_ax.set_title("Largest joint SE multiplier retaining the log-ratio rejection", loc="left", fontsize=10)
surface_fig.text(.02, .015,
                 "Bonferroni across 51 variants. Rejection requires a multiplier below the printed boundary.\n"
                 "Values below 1 fail at the reported SEs. Effective degrees of freedom are scenarios, not measured batch counts.",
                 fontsize=7, color=MUTED)
surface_fig.tight_layout(rect=(0,.09,1,1))
surface_fig.savefig("SCN5A_robustness.png", dpi=300)
surface_fig.savefig("SCN5A_robustness.pdf")
plt.close(surface_fig)
print("\nwrote scn5a_specifications.{csv,json,log}, Figure1.{png,pdf}")
