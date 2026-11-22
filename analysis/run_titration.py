#!/usr/bin/env python3
"""Inversion of published mixed-ratio titrations (revision of 21 November 2026).

    python3 run_titration.py

Applies the identification result of the main text to two published
titrations digitized from their figures (data/titration/):

  * Hosy and Vivaudou (2014): KATP channels with four SUR2A subunits and six
    interior wild-type:mutant cRNA proportions (n = 4, three unknown profile
    entries, three surplus fractions);
  * Clatot et al. (2017): NaV1.5 L325R and R104W titrated against wild type
    (treated as a dimer, n = 2, one unknown entry, several surplus fractions).

For each series, measured responses are converted to transfer ratios
d(p) = (R0 - R(p)) / (R0 - R1) with the pure wild-type (R0) and pure mutant
(R1) endpoints, where p is the input mutant fraction. The unrestricted
interior profile is estimated by generalized least squares with a Monte Carlo
covariance matrix that includes the shared endpoints, and a chi-square test
asks whether one profile fits every interior fraction. Candidate profiles
considered by the original authors are compared on the same scale.

The input fractions are cRNA or cDNA proportions, not measured
assembly-competent fractions, and digitized values carry reading error; the
analysis therefore demonstrates the calculation under the paper's stated
assumptions rather than recovering established molecular activities.
Writes titration.json and Titration.{png,pdf} (S2 Fig).
"""
import json
from datetime import date
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "titration" / "titration_digitized.csv"
NMC = 20000
SEED = 20261111
SUBUNITS = {"Hosy": 4, "Clatot": 2}             # SUR2A in the KATP octamer; NaV1.5 treated as a dimer
DIGITIZATION_FLOOR = {"Hosy": 0.01, "Clatot": 0.0003}   # inter-rendering scatter (raster); vector resolution
CANDIDATE_FAMILY_SIZE = 26
SERIES_FAMILY_SIZE = 7
FRACTION_SHIFTS = (-0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05)
MUTANT_ENDPOINT_SHIFTS = (0.0, 0.01, 0.02, 0.05, 0.10)  # fraction of the WT endpoint
PANELS = {"Hosy": "Fig 3I"}                     # Fig 4B is derived from Fig 3I and is not pooled


def se_of(row, floor):
    """Standard error from plotted error bars plus digitization uncertainty.

    Error-bar ends are given as err_low and err_high; where only one end is
    drawn (Hosy) or one end is clipped (Clatot S460A), the other is mirrored."""
    lo, hi, v = row.err_low, row.err_high, row.value
    if np.isfinite(lo) and np.isfinite(hi):
        half = (hi - lo) / 2
    elif np.isfinite(hi):
        half = hi - v
    elif np.isfinite(lo):
        half = v - lo
    else:
        half = 0.0
    if str(row.err_type).upper().startswith("SD") and np.isfinite(row.n) and row.n > 1:
        half = half / np.sqrt(row.n)
    res = max(row.resolution_data_units if np.isfinite(row.resolution_data_units) else 0.0, floor)
    return float(np.sqrt(half ** 2 + res ** 2))


def design(p, n):
    return np.array([[comb(n, k) * q ** k * (1 - q) ** (n - k) for k in range(1, n)] for q in p])


def candidates(n):
    if n == 2:
        return {"heterodimer inactive (g1 = 1)": [1.0],
                "independent subunits (g1 = 0.5)": [0.5],
                "heterodimer fully active (g1 = 0)": [0.0]}
    out = {}
    for j in range(1, n + 1):          # activity requires at least j wild-type subunits
        out[f"requires >= {j} wild-type subunits"] = [0.0 if n - k >= j else 1.0 for k in range(1, n)]
    out["proportional independent-subunits (g_k = k/n)"] = [k / n for k in range(1, n)]
    if n == 4:
        # Two wild-type subunits suffice only in some arrangements. Hosy and
        # Vivaudou weighted adjacent and opposite pairs equally ("2*"); random
        # assembly gives 4 adjacent and 2 opposite arrangements of two subunits.
        out["requires >= 2 wild-type, equal arrangement weights (Hosy 2*)"] = [0.0, 0.5, 1.0]
        out["requires >= 2 wild-type, adjacent pairs only (random placement)"] = [0.0, 1 / 3, 1.0]
        out["requires >= 2 wild-type, opposite pairs only (random placement)"] = [0.0, 2 / 3, 1.0]
    return out


def fit_series(df, n, rng, floor, *, fraction_shift=0.0,
               mutant_endpoint_shift=0.0):
    df = df.sort_values("mutant_fraction_input")
    r0 = df[np.isclose(df.mutant_fraction_input, 0.0)]
    r1 = df[np.isclose(df.mutant_fraction_input, 1.0)]
    inner = df[(df.mutant_fraction_input > 0) & (df.mutant_fraction_input < 1)]
    if len(r0) != 1 or len(r1) != 1 or len(inner) < n - 1:
        raise ValueError("series needs one wild-type endpoint, one mutant endpoint, and >= n-1 interior fractions")
    p = inner.mutant_fraction_input.to_numpy(float) + fraction_shift
    if not np.all((p > 0) & (p < 1)):
        raise ValueError("shifted interior fractions must stay strictly between zero and one")
    vals = np.concatenate([[r0.value.iloc[0], r1.value.iloc[0]], inner.value.to_numpy(float)])
    vals[1] += mutant_endpoint_shift * vals[0]
    ses = np.array([se_of(r, floor) for r in pd.concat([r0, r1, inner]).itertuples()])
    draws = vals + ses * rng.standard_normal((NMC, len(vals)))
    D = (draws[:, :1] - draws[:, 2:]) / (draws[:, :1] - draws[:, 1:2])
    d = (vals[0] - vals[2:]) / (vals[0] - vals[1])
    cov = np.cov(D, rowvar=False)
    W = np.linalg.inv(cov)
    X = design(p, n)
    target = d - p ** n
    A = X.T @ W @ X
    g = np.linalg.solve(A, X.T @ W @ target)
    resid = target - X @ g
    chi2 = float(resid @ W @ resid)
    dof = len(p) - (n - 1)
    # Monte Carlo distribution of the GLS estimate (fixed weights).
    G = np.linalg.solve(A, X.T @ W @ (D - p ** n).T).T
    # Each candidate is tested against the unrestricted fit: the increase in
    # chi2 when the profile is fixed at the candidate, with n-1 degrees of
    # freedom. Candidate interpretation is conditional on one profile fitting
    # all fractions; pointwise MC intervals need not match adjusted tests. The candidate's own
    # goodness of fit over all J fractions is kept for reference.
    cand = {}
    for name, prof in candidates(n).items():
        pred = p ** n + X @ np.array(prof)
        r = d - pred
        c2 = float(r @ W @ r)
        delta = c2 - chi2
        cand[name] = {"profile": prof, "delta_chi2": delta, "dof": int(n - 1),
                      "p_value": float(stats.chi2.sf(max(delta, 0), n - 1)),
                      "p_bonferroni_26": min(1.0, CANDIDATE_FAMILY_SIZE * float(stats.chi2.sf(max(delta, 0), n - 1))),
                      "fit_chi2": c2, "fit_dof": int(len(p)),
                      "fit_p_value": float(stats.chi2.sf(c2, len(p)))}
    out = {"n_subunits": n,
           "assumptions": {"fraction_shift": fraction_shift,
                           "mutant_endpoint_shift_relative_to_WT": mutant_endpoint_shift,
                           "mutant_endpoint_value": float(vals[1])},
           "mutant_fractions": p.tolist(), "d_observed": d.tolist(),
           "d_se": np.sqrt(np.diag(cov)).tolist(),
           "d_covariance": cov.tolist(),
           "profile_estimate": g.tolist(),
           "profile_ci": [np.quantile(G[:, k], [0.025, 0.975]).tolist() for k in range(n - 1)],
           "one_profile_fit": {"chi2": chi2, "dof": int(dof),
                               "p_value": float(stats.chi2.sf(chi2, dof)) if dof > 0 else None,
                               "p_bonferroni_7": min(1.0, SERIES_FAMILY_SIZE * float(stats.chi2.sf(chi2, dof))) if dof > 0 else None,
                               "lack_of_fit_after_bonferroni": bool(SERIES_FAMILY_SIZE * stats.chi2.sf(chi2, dof) < .05) if dof > 0 else None},
           "candidates": cand}
    if n == 2:
        out["g1_by_fraction"] = ((d - p ** 2) / (2 * p * (1 - p))).tolist()
    return out


def short_label(key, construct):
    if key == "Hosy":
        return "Hosy & Vivaudou: SUR2A wild type + DNN"
    c = construct.replace("Nav1.5 ", "").split(" (")[0]
    extra = " + difopein" if "difopein" in construct else ""
    return f"Clatot et al.: {c}{extra}"


def coded_fractions(df, key, coding):
    """Use stated input fractions, or a literal-ratio sensitivity.

    Supplementary Table 1 and Supplementary Fig. 3 state fractions .10, .25,
    .50, .75 and .90, matching the plotted coordinates. The text's 1:10, 1:4,
    4:1 and 10:1 labels are shorthand for these mixtures. Converting those
    labels literally to 1/11, .20, .80 and 10/11 is a sensitivity only.
    Legacy machine keys are retained to preserve existing archived references.
    """
    df = df.copy()
    if key == "Clatot" and coding == "stated_cDNA_ratios":
        df["mutant_fraction_input"] = df.mutant_fraction_input.replace({0.1: 1 / 11, 0.25: 0.20, 0.75: 0.80, 0.9: 10 / 11})
        df["wt_fraction_input"] = 1 - df.mutant_fraction_input
    return df


def compact_sensitivity(result):
    return {k: result[k] for k in ("assumptions", "mutant_fractions", "d_covariance", "profile_estimate",
                                   "profile_ci", "one_profile_fit", "candidates")}



PRIMARY_CODING = "plotted_positions"
SENSITIVITY_CODING = "stated_cDNA_ratios"
CODINGS = (PRIMARY_CODING, SENSITIVITY_CODING)
CODING_LABELS = {PRIMARY_CODING: "Source-stated input fractions",
                 SENSITIVITY_CODING: "Literal-ratio sensitivity"}
CODING_ROLES = {PRIMARY_CODING: "primary documented input design",
                SENSITIVITY_CODING: "sensitivity to literal conversion of ratio shorthand"}


def holm_adjust(pvalues):
    """Holm step-down familywise adjusted probabilities, in input order."""
    pvalues = np.asarray(pvalues, float)
    order = np.argsort(pvalues, kind="stable")
    sorted_adjusted = np.minimum(1.0, np.maximum.accumulate((len(order) - np.arange(len(order))) * pvalues[order]))
    out = np.empty(len(order))
    out[order] = sorted_adjusted
    return out


def adjust_family(fits):
    """Correct 26 candidates and seven fit tests within a fixed scenario."""
    candidates_flat = [c for fit in fits.values() for c in fit["candidates"].values()]
    if len(candidates_flat) != CANDIDATE_FAMILY_SIZE or len(fits) != SERIES_FAMILY_SIZE:
        raise ValueError("Every correction family must retain 26 candidates and seven series")
    for c, adjusted in zip(candidates_flat, holm_adjust([c["p_value"] for c in candidates_flat])):
        c["p_holm_26"] = float(adjusted)
    for fit, adjusted in zip(fits.values(), holm_adjust([f["one_profile_fit"]["p_value"] for f in fits.values()])):
        f = fit["one_profile_fit"]
        f["p_holm_7"] = float(adjusted)
        f["lack_of_fit_after_holm"] = bool(adjusted < .05)


def scenario_key(coding, shift, endpoint):
    return f"{coding}; fraction_shift={shift:+.2f}; mutant_endpoint_shift={endpoint:.2f}"


def decision_ranges(scenarios, coding, candidate, correction, endpoint_grid):
    """Contiguous tested decision stability; no search for favorable islands.

    For candidates the state is (one-profile rejection, candidate rejection).
    For fits it is only one-profile rejection. Ranges compare each grid point
    with baseline, and therefore describe decisions, not a validated biology.
    """
    fit_key = f"p_{correction}_7"
    candidate_key = f"p_{correction}_26"

    def state(shift, endpoint):
        fit = scenarios[scenario_key(coding, shift, endpoint)]
        values = [fit["one_profile_fit"][fit_key] < .05]
        if candidate is not None:
            values.append(fit["candidates"][candidate][candidate_key] < .05)
        return tuple(values)

    baseline = state(0.0, 0.0)
    neg, pos = 0.0, 0.0
    for shift in sorted((q for q in FRACTION_SHIFTS if q < 0), reverse=True):
        if state(shift, 0.0) != baseline:
            break
        neg = shift
    for shift in sorted(q for q in FRACTION_SHIFTS if q > 0):
        if state(shift, 0.0) != baseline:
            break
        pos = shift
    max_endpoint = 0.0
    for endpoint in endpoint_grid:
        if state(0.0, endpoint) != baseline:
            break
        max_endpoint = endpoint
    frontier = []
    for radius in sorted({abs(q) for q in FRACTION_SHIFTS}):
        stable_endpoint = None
        for endpoint in endpoint_grid:
            stable = all(state(q, e) == baseline for q in FRACTION_SHIFTS if abs(q) <= radius
                         for e in endpoint_grid if e <= endpoint)
            if not stable:
                break
            stable_endpoint = endpoint
        frontier.append({"symmetric_fraction_radius": radius,
                         "max_endpoint_increment_preserving_every_included_grid_point": stable_endpoint})
    return {
        "correction": correction, "baseline_profile_rejected": baseline[0],
        "baseline_candidate_rejected": baseline[1] if candidate is not None else None,
        "baseline_candidate_interpretable": not baseline[0] if candidate is not None else None,
        "fraction_shift_low_at_zero_endpoint": neg,
        "fraction_shift_high_at_zero_endpoint": pos,
        "largest_endpoint_increment_at_zero_fraction_shift": max_endpoint,
        "all_grid_decisions_match_baseline": all(state(q, e) == baseline for q in FRACTION_SHIFTS for e in endpoint_grid),
        "joint_rectangle_frontier": frontier,
    }


def replication_rows(df, key):
    rows = []
    for r in df.sort_values("mutant_fraction_input").itertuples():
        row = {"plotted_fraction": float(r.mutant_fraction_input), "n_recordings": int(r.n),
               "unit": "inside-out macropatches from Xenopus oocytes" if key == "Hosy" else "whole-cell HEK293 recordings (source labels n as experiments)",
               "error_bar": str(r.err_type), "reported_sem_available": bool(np.isfinite(r.err_low) or np.isfinite(r.err_high)),
               "independent_batches": "not reported by condition", "digitization_resolution": float(r.resolution_data_units)}
        rows.append(row)
    return rows


def main():
    raw = pd.read_csv(DATA)
    for c in ("err_low", "err_high", "n", "resolution_data_units", "value", "wt_fraction_input", "mutant_fraction_input"):
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw["mutant_fraction_input"] = raw["mutant_fraction_input"].fillna(1 - raw["wt_fraction_input"])
    results, frames = {}, {}
    for i, ((study, panel, construct), df) in enumerate(raw.groupby(["study", "figure_panel", "construct"], sort=False)):
        key = "Hosy" if study.startswith("Hosy") else "Clatot"
        if key in PANELS and panel != PANELS[key]:
            continue
        label = short_label(key, construct)
        results[label] = {"study": study, "figure_panel": panel, "construct": construct,
                          "monte_carlo_seed": SEED + i, "measure": str(df.measure.iloc[0]),
                          "n_subunits": SUBUNITS[key], "replication": replication_rows(df, key),
                          "codings": {}, "sensitivity": {}}
        frames[label] = (key, df)
    # The same seven series and all 26 candidates enter each coding/scenario
    # family. KATP coding is unchanged and its endpoint is never shifted.
    for coding in CODINGS:
        for shift in FRACTION_SHIFTS:
            for endpoint in MUTANT_ENDPOINT_SHIFTS:
                fits = {}
                for label, (key, df) in frames.items():
                    endpoint_used = endpoint if key == "Clatot" else 0.0
                    fits[label] = fit_series(coded_fractions(df, key, coding), SUBUNITS[key],
                        np.random.default_rng(results[label]["monte_carlo_seed"]), DIGITIZATION_FLOOR[key],
                        fraction_shift=shift, mutant_endpoint_shift=endpoint_used)
                adjust_family(fits)
                for label, fit in fits.items():
                    key, _ = frames[label]
                    # For KATP, Holm still depends on the Clatot endpoint scenario:
                    # store the full family context even though its own endpoint is fixed.
                    fit["family_scenario"] = {"coding": coding, "coding_label": CODING_LABELS[coding],
                                               "coding_role": CODING_ROLES[coding], "fraction_shift": shift,
                                               "Clatot_endpoint_increment": endpoint}
                    results[label]["sensitivity"][scenario_key(coding, shift, endpoint)] = fit
                    if shift == 0.0 and endpoint == 0.0:
                        results[label]["codings"][coding] = fit
    profile_rows, candidate_rows, range_rows, fit_range_rows, replicate_rows = [], [], [], [], []
    for label, r in results.items():
        for row in r["replication"]:
            replicate_rows.append({"series": label, **row})
        r["decision_stability"] = {}
        endpoint_grid = (0.0,) if frames[label][0] == "Hosy" else MUTANT_ENDPOINT_SHIFTS
        for coding, fit in r["codings"].items():
            f = fit["one_profile_fit"]
            profile_rows.append({"series": label, "coding": coding, "coding_label": CODING_LABELS[coding],
                "coding_role": CODING_ROLES[coding], "n_subunits": r["n_subunits"],
                "fractions": ", ".join(f"{v:.6g}" for v in fit["mutant_fractions"]),
                "profile": "; ".join(f"{v:.3f} ({lo:.3f}, {hi:.3f})" for v, (lo, hi) in zip(fit["profile_estimate"], fit["profile_ci"])),
                **f})
            for name, c in fit["candidates"].items():
                candidate_rows.append({"series": label, "coding": coding, "coding_label": CODING_LABELS[coding],
                    "coding_role": CODING_ROLES[coding], "candidate": name,
                    "delta_chi2": c["delta_chi2"], "df": c["dof"], "nominal_P": c["p_value"],
                    "Bonferroni_P": c["p_bonferroni_26"], "Holm_P": c["p_holm_26"],
                    "profile_rejected_Bonferroni": f["lack_of_fit_after_bonferroni"],
                    "profile_rejected_Holm": f["lack_of_fit_after_holm"]})
            for correction in ("bonferroni", "holm"):
                fit_range = decision_ranges(r["sensitivity"], coding, None, correction, endpoint_grid)
                fit_range_rows.append({"series": label, "coding": coding, "coding_label": CODING_LABELS[coding],
                    "coding_role": CODING_ROLES[coding], "endpoint_sensitivity_applicable": frames[label][0] == "Clatot", **fit_range})
                for name in fit["candidates"]:
                    ranges = decision_ranges(r["sensitivity"], coding, name, correction, endpoint_grid)
                    r["decision_stability"][f"{coding}; {correction}; {name}"] = ranges
                    range_rows.append({"series": label, "coding": coding, "coding_label": CODING_LABELS[coding],
                        "coding_role": CODING_ROLES[coding], "candidate": name,
                        "endpoint_sensitivity_applicable": frames[label][0] == "Clatot", **ranges})
    out = {"generated": date.today().isoformat(), "n_mc": NMC, "seed": SEED,
        "input": str(DATA.relative_to(HERE)), "coding_order": list(CODINGS),
        "primary_coding": PRIMARY_CODING, "coding_labels": CODING_LABELS, "coding_roles": CODING_ROLES,
        "fraction_definition": "Clatot Supplementary Table 1 and Supplementary Fig. 3 state input mutant fractions 0.10, 0.25, 0.50, 0.75 and 0.90, matching the plotted coordinates. The ratio labels in the text are shorthand for these percentages. Literal conversion to 1/11, 0.20, 0.50, 0.80 and 10/11 is a sensitivity, retained under the legacy stated_cDNA_ratios key. Input fractions are assumed to represent the assembly pool. KATP uses its stated cRNA fractions in both comparison families.",
        "unit_of_analysis": "Recordings treated as independent; digitized SEMs are treated as known error scales.",
        "source_locators": {"Hosy": "Fig. 3I; current in 300 uM MgADP divided by nucleotide-free current, inside-out oocyte macropatches at -50 mV",
                            "Clatot": "Fig. 1c and Supplementary Fig. 3C; peak current density at -20 mV normalized to each series' WT or S460A endpoint. Supplementary Table 1 supplies recording counts only, not means or SEMs."},
        "digitization_error_floor": DIGITIZATION_FLOOR,
        "fit_degrees_of_freedom": "J - (n_subunits - 1), where J is the number of interior fractions; candidate contrasts use n_subunits - 1.",
        "candidate_multiplicity": {"methods": ["Bonferroni", "Holm"], "family_size": 26,
            "scope": "all candidates across seven series within each coding and fixed scenario; none removed because the one-profile fit fails"},
        "profile_fit_multiplicity": {"methods": ["Bonferroni", "Holm"], "family_size": 7},
        "sensitivity_grid": {"interior_fraction_shifts": list(FRACTION_SHIFTS),
            "Clatot_mutant_endpoint_additions_relative_to_WT": list(MUTANT_ENDPOINT_SHIFTS),
            "boundary_interpretation": "Largest contiguous tested shifts/increments preserving the baseline joint fit/candidate decision, not exact continuous thresholds. Joint rectangle frontiers require stability at every included grid point."},
        "limits": ["The source states its input proportions; identifying them with assembly-competent fractions remains an assumption.",
            "Recordings treated as independent. Digitized SEMs and endpoint errors are treated as known scales; the unknown number of independent cultures, batches, or days cannot be recovered from cell/patch counts.",
            "Both Bonferroni and Holm results are reported. Correction-dependent decisions do not establish compatibility or incompatibility of a biological mechanism.",
            "A profile rejection rejects the conjunction of assumed subunit number, fractions, fixed abundance, linear readout and one profile, not a profile independent of those assumptions.",
            "Scenario grids are descriptive checks, not measured uncertainty distributions; no favorable isolated scenario is selected.",
            "The proportional g_k=k/n line is also an independent-subunits response and cannot distinguish coassembly by itself.",
            "Candidate non-rejection does not establish equivalence; inference remains conditional on an adequate one-profile fit."],
        "series": results}
    (HERE / "titration.json").write_text(json.dumps(out, indent=2) + "\n")
    for name, rows in (("titration_profiles", profile_rows), ("titration_candidate_tests", candidate_rows),
                       ("titration_sensitivity_ranges", range_rows), ("titration_fit_sensitivity_ranges", fit_range_rows),
                       ("titration_replication", replicate_rows)):
        pd.DataFrame(rows).to_csv(HERE / f"{name}.csv", index=False)
    figure(results)
    print(pd.DataFrame(profile_rows).to_string(index=False))
    print("\nAll 52 baseline candidate tests (26 per coding):")
    print(pd.DataFrame(candidate_rows).to_string(index=False))
    print("\nWrote titration.json, five result CSVs and Titration.png/pdf/tiff.")


def figure(results):
    """Source-stated input fractions, with literal-ratio fits as sensitivity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from palette import BLUE_DARK, FIG_WIDTH, GOLD_LIGHT, MUTED, TEAL, halo, panel_label, use_figure_style
    use_figure_style()
    fig, axes = plt.subplots(4, 2, figsize=(FIG_WIDTH, 9.1), squeeze=False)
    grid = np.linspace(0, 1, 201)
    for i, (label, r) in enumerate(results.items()):
        ax = axes.flat[i]
        n = r["n_subunits"]
        # Candidate curves are independent of the fraction coding.
        for c in r["codings"][PRIMARY_CODING]["candidates"].values():
            ax.plot(grid, grid ** n + design(grid, n) @ np.array(c["profile"]), color=GOLD_LIGHT, lw=.8, alpha=.7)
        fit = r["codings"][PRIMARY_CODING]
        ax.plot(grid, grid ** n + design(grid, n) @ np.array(fit["profile_estimate"]), color=BLUE_DARK, lw=1.5)
        ax.errorbar(fit["mutant_fractions"], fit["d_observed"], yerr=1.96*np.array(fit["d_se"]),
            fmt="o", color=BLUE_DARK, ms=3.8, capsize=1.5, lw=.7, zorder=3)
        f = fit["one_profile_fit"]
        notes = [f"Stated inputs: χ²={f['chi2']:.2f} ({f['dof']} df)",
                 f"P(B/H)={f['p_bonferroni_7']:.3g}/{f['p_holm_7']:.3g}"]
        if i != 0:
            alt = r["codings"][SENSITIVITY_CODING]
            ax.plot(grid, grid ** n + design(grid, n) @ np.array(alt["profile_estimate"]), color=TEAL, lw=1.0, ls="--")
        title = label.split(": ", 1)[1]
        panel_label(ax, "ABCDEFG"[i], title, size=8, letter_size=9)
        ax.text(.02, .97, f"$n_{{subunits}}={n}$", transform=ax.transAxes,
                fontsize=7, va="top", color=MUTED, path_effects=halo(2))
        ax.text(.98, .03, "\n".join(notes), transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6.7, color=MUTED, path_effects=halo(2), linespacing=1.3)
        ax.set(xlim=(-.02, 1.02), ylim=(-.25,1.1), xlabel="input mutant fraction $p$", ylabel="transfer ratio $d(p)$")
        ax.tick_params(labelsize=7)
    leg=axes.flat[-1]
    leg.axis("off")
    leg.legend(handles=[Line2D([0],[0],color=BLUE_DARK,marker="o",label="Source-stated inputs and fitted profile"),
        Line2D([0],[0],color=TEAL,ls="--",label="Literal-ratio sensitivity fit"),
        Line2D([0],[0],color=GOLD_LIGHT,label="Candidate profiles")],loc="upper left",fontsize=8)
    fig.tight_layout(pad=.6,h_pad=1.1,w_pad=1.1)
    fig.savefig(HERE/"Titration.png",dpi=300)
    fig.savefig(HERE/"Titration.pdf")
    fig.savefig(HERE/"Titration.tiff",dpi=300,pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    main()
