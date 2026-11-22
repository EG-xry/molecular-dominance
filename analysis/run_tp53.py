

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Section 3: are transfer slopes portable across studies?

    python3 run_tp53.py

Two wild-type-background TP53 screens are regressed on a common alone-state
axis: Giacomelli's p53-null nutlin-3 panel (urn:mavedb:00000068-b-1). Their
matched-cohort comparison uses the missense variants present in these two
mixed-state screens and the alone-state screen.

HCT116 mut/Delta is a monoallelic mutant state. Its comparison is reported
separately on its own coverage as a cross-study alone-state comparison; it
does not define or belong to the mixed-screen common cohort.

Nonsense and stop-loss are excluded from every regression panel. Both axes are the
screens' published scores, which are separately normalised, so these slopes
are descriptive score relationships, not calibrated dominance coefficients.
Synonymous-control slopes are accompanied by centered statistics because
score offsets can produce a through-origin slope without an association.

Separately, synonymous/nonsense diagnostics assess the applicability of
Funk et al.'s two-anchor scaling to these exact deposited screens. Nonsense
controls are not assumed to define maximal dominant-negative activity in a
WT-background assay. No two-anchor mixed-screen calibration is imposed.

Writes tp53_portability.csv/json, tp53_controls.csv, and tp53_anchors.csv/json.
"""
import json
import re
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

SRC = _data("tp53", "../tp53-transfer")
rng = np.random.default_rng(20260928)
NBOOT = 20000
NANCHORBOOT = 5000
ANCHOR_SEED = 20261007

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}


def parse(h):
    """HGVS protein -> (key, position, synonymous). Ter is not an amino acid."""
    if not isinstance(h, str):
        return None
    m = re.fullmatch(r"(?:[A-Za-z0-9_.]+:)?p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)", h)
    if not m:
        return None
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if wt == "Ter" or mut == "Ter" or wt not in AA3:
        return None
    if mut == "=" or mut == wt:
        return f"{AA3[wt]}{pos}=", pos, True
    if mut not in AA3:
        return None
    return f"{AA3[wt]}{pos}{AA3[mut]}", pos, False


def load(fname, procol="hgvs_pro"):
    d = pd.read_csv(f"{SRC}/{fname}")
    p = d[procol].map(parse)
    k = p.notna()
    out = pd.DataFrame({"key": [v[0] for v in p[k]],
                        "pos": [v[1] for v in p[k]],
                        "syn": [v[2] for v in p[k]],
                        "score": d.loc[k, "score"].values}).dropna(subset=["score"])
    # HCT116 tiles at nucleotide level, so several variants share one protein
    # change; without this they would be counted repeatedly and their residues
    # over-weighted.
    return (out.groupby("key", as_index=False)
            .agg(pos=("pos", "first"), syn=("syn", "first"),
                 score=("score", "mean"), n_nt=("score", "size")))


def parse_anchor(h):
    """Classify single protein substitutions, retaining stops for diagnostics only.

    Frameshifts, multi-substitutions, uncertain calls, and stop-loss variants
    are excluded. Amino-acid synonyms use the same residue key as load().
    """
    parsed = parse(h)
    if parsed is not None:
        key, pos, syn = parsed
        return key, pos, "synonymous" if syn else "missense"
    if not isinstance(h, str):
        return None
    match = re.fullmatch(r"(?:[A-Za-z0-9_.]+:)?p\.([A-Z][a-z]{2})(\d+)Ter", h)
    if match is not None and match.group(1) in AA3:
        wt, pos = AA3[match.group(1)], int(match.group(2))
        return f"{wt}{pos}*", pos, "nonsense"
    return None


def load_anchor_variants(fname, procol="hgvs_pro"):
    """Read finite deposited scores and average nucleotide alleles by protein key.

    This matches the protein-level weighting of the regressions; it is not
    a reconstruction of Funk's original per-exon, per-replicate normalization.
    """
    source = pd.read_csv(f"{SRC}/{fname}")
    parsed = source[procol].map(parse_anchor)
    scores = pd.to_numeric(source["score"], errors="coerce")
    keep = parsed.notna() & np.isfinite(scores)
    frame = pd.DataFrame({"key": [v[0] for v in parsed[keep]],
                          "pos": [v[1] for v in parsed[keep]],
                          "variant_class": [v[2] for v in parsed[keep]],
                          "score": scores[keep].to_numpy()})
    return (frame.groupby("key", as_index=False)
            .agg(pos=("pos", "first"), variant_class=("variant_class", "first"),
                 score=("score", "mean"), n_input=("score", "size")))


def anchor_gap_ci(frame, generator, n_boot=NANCHORBOOT):
    """Resample residues jointly across the two control classes.

    Controls absent at a sampled residue remain absent. The interval is
    conditional on deposited scores, not on raw-count or replicate uncertainty.
    """
    controls = frame[frame.variant_class.isin(["synonymous", "nonsense"])]
    wide = controls.pivot(index="pos", columns="variant_class", values="score")
    values = wide.reindex(columns=["synonymous", "nonsense"]).to_numpy()
    samples = values[generator.integers(0, len(values), size=(n_boot, len(values)))]
    valid = np.isfinite(samples).any(axis=1).all(axis=1)
    medians = np.nanmedian(samples[valid], axis=1)
    gaps = medians[:, 1] - medians[:, 0]
    return (*np.quantile(gaps, [0.025, 0.975]).astype(float), int(valid.sum()))


def assess_anchors(n_boot=NANCHORBOOT, seed=ANCHOR_SEED):
    """Describe candidate anchors without treating class separation as calibration.

    Both mixed screens measure suppression of WT activity. Truncating TP53
    variants can lack the ability to inhibit WT, so their LOF state does not
    establish the same functional endpoint as a strong dominant-negative
    missense variant. A negative gap is diagnostic of these deposited controls,
    not a general mathematical prohibition on affine rescaling.
    """
    alone_name = "Giacomelli, p53-null + nutlin-3"
    specs = {alone_name: ("mavedb_00000068-b-1.csv", "hgvs_pro"),
             **SCREENS, **ALONE_STATE_SCREENS}
    generator = np.random.default_rng(seed)
    rows = []
    for name, (filename, column) in specs.items():
        frame = load_anchor_variants(filename, column)
        for scope, panel in [("all deposited positions", frame),
                             ("Funk coding region 126-307",
                              frame[frame.pos.between(126, 307)])]:
            syn = panel.loc[panel.variant_class == "synonymous", "score"]
            stop = panel.loc[panel.variant_class == "nonsense", "score"]
            if len(syn) < 2 or len(stop) < 2:
                raise ValueError(f"insufficient anchor controls: {name}, {scope}")
            lo, hi, n_valid = anchor_gap_ci(panel, generator, n_boot)
            rows.append({"screen": name,
                         "accession": "urn:mavedb:" + filename.removeprefix("mavedb_").removesuffix(".csv"),
                         "state": ("mixed" if name in SCREENS else "alone"),
                         "scope": scope,
                         "n_synonymous": len(syn), "n_nonsense": len(stop),
                         "synonymous_median": float(syn.median()),
                         "synonymous_q25": float(syn.quantile(0.25)),
                         "synonymous_q75": float(syn.quantile(0.75)),
                         "nonsense_median": float(stop.median()),
                         "nonsense_q25": float(stop.quantile(0.25)),
                         "nonsense_q75": float(stop.quantile(0.75)),
                         "nonsense_minus_synonymous": float(stop.median() - syn.median()),
                         "gap_ci_lo": lo, "gap_ci_hi": hi,
                         "n_valid_boot": n_valid,
                         "prob_nonsense_above_synonymous": float(
                             stats.mannwhitneyu(stop, syn).statistic / (len(stop) * len(syn))),
                         "anchor_interpretation": (
                             "not validated as a common dominant-negative endpoint"
                             if name in SCREENS else
                             "LOF reference-class contrast; not molecular calibration")})
    return rows


def boot_ci(x, y, pos, n=NBOOT):
    groups = [np.array(g) for g in pd.Series(range(len(x))).groupby(pos.values).apply(list)]
    out = []
    for _ in range(n):
        s = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        out.append(np.sum(x[s] * y[s]) / np.sum(x[s] ** 2))
    return float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


# MIXED-STATE screens only. An earlier version of this analysis included
# urn:mavedb:00001213 (HCT116) as a third wild-type-background screen. It is
# not one: the authors "haploidized the cells for TP53 by deleting intronic
# splicing branch points in one of the two TP53 alleles" and then induced
# "monoallelic p53 variant expression (HCT116 mut/Delta)". The variant is
# expressed over a deleted allele, so no wild-type protein is present and the
# measurement is an ALONE state, not a mixed one. It is reported below as a
# cross-study alone-state comparison, which is a different question.
SCREENS = {
    "Giacomelli, WT background + nutlin-3": ("mavedb_00000068-a-1.csv", "hgvs_pro"),
    "Boettcher, K562 WT-activity reporter": ("mavedb_00001235-a-1.csv", "hgvs_pro"),
}
ALONE_STATE_SCREENS = {
    "HCT116, mut/Delta monoallelic": ("mavedb_00001213-a-1.csv", "HGVS(protein)"),
}

def cohorts(alone, mixed, alone_state):
    """Keep the two mixed screens' common missense cohort distinct from HCT116."""
    common = set(alone.loc[~alone.syn, "key"])
    for frame in mixed.values():
        common &= set(frame.loc[~frame.syn, "key"])
    result = []
    for scope in ("own coverage", "common variant set"):
        screens = {**mixed, **alone_state} if scope == "own coverage" else mixed
        for name, frame in screens.items():
            paired = alone.merge(frame, on="key", suffixes=("_x", "_y"),
                                 validate="one_to_one").dropna(subset=["score_x", "score_y"])
            if scope == "common variant set":
                paired = paired[paired.key.isin(common)]
            result.append((scope, name, paired))
    return common, result


def main():
    ALONE = load("mavedb_00000068-b-1.csv")
    loaded = {name: load(fn, col) for name, (fn, col) in SCREENS.items()}
    alone_state = {name: load(fn, col) for name, (fn, col) in ALONE_STATE_SCREENS.items()}
    common, panels = cohorts(ALONE, loaded, alone_state)
    rows = []
    for scope, name, j in panels:
        mis = j[~j.syn_x].reset_index(drop=True)
        x, y = mis.score_x.values, mis.score_y.values
        b = float(np.sum(x * y) / np.sum(x * x))
        lo, hi = boot_ci(x, y, mis.pos_x)
        rows.append({"scope": scope, "screen": name,
                     "state": ("mixed" if name in loaded else "alone (mut/Delta)"),
                     "n": len(mis),
                     "n_residues": int(mis.pos_x.nunique()),
                     "b": round(b, 4), "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                     "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 4)})

    # Synonymous variants provide a diagnostic of the origin and covariance on
    # each pair of score axes. They are not established zero-effect observations.
    ctl_rows = []
    for name, m in {**loaded, **alone_state}.items():
        j = ALONE.merge(m, on="key", suffixes=("_x", "_y"), validate="one_to_one").dropna(
            subset=["score_x", "score_y"])
        c = j[j.syn_x & j.syn_y]
        if len(c) < 20:
            ctl_rows.append({"screen": name, "n_control": len(c),
                             "b_control": None, "note": "too few synonymous controls"})
            continue
        x, y = c.score_x.values, c.score_y.values
        b = float(np.sum(x * y) / np.sum(x * x))
        lo, hi = boot_ci(x, y, c.pos_x)
        # A through-origin slope is E[xy]/E[x^2], and E[xy] = Cov(x,y) + E[x]E[y].
        # Where the axes have nonzero means the slope can be large with no
        # association at all, so the CENTRED statistics are reported beside it --
        # the same distinction this paper makes about panel slopes, applied to its
        # own control diagnostic.
        r_c, p_c = stats.pearsonr(x, y)
        slope_c = float(np.cov(x, y)[0, 1] / np.var(x, ddof=1))
        mean_product = float(x.mean() * y.mean())
        ctl_rows.append({"screen": name,
                         "mean_x": round(float(x.mean()), 4),
                         "mean_y": round(float(y.mean()), 4),
                         "mean_product_share_of_Exy":
                             round(mean_product / float(np.mean(x * y)), 3),
                         "centred_pearson_r": round(float(r_c), 4),
                         "centred_p": round(float(p_c), 4),
                         "centred_slope": round(slope_c, 4),
                         "pairing": ("within study" if "Giacomelli" in name
                                     else "across studies"),
                         "n_control": len(c), "b_control": round(b, 4),
                         "ci_lo": round(lo, 4), "ci_hi": round(hi, 4)})
    ctl = pd.DataFrame(ctl_rows)
    ctl.to_csv("tp53_controls.csv", index=False)

    df = pd.DataFrame(rows)
    df.to_csv("tp53_portability.csv", index=False)
    json.dump({"generated": date.today().isoformat(),
               "common_alone_axis": "urn:mavedb:00000068-b-1",
               "n_boot": NBOOT, "seed": 20260928,
               "common_mixed_missense_n": len(common),
               "common_mixed_definition": "intersection of two mixed screens and Giacomelli alone axis; HCT116 excluded",
               "bootstrap_scope": "residues conditional on published normalized scores",
               "input_files": {"alone_axis": "mavedb_00000068-b-1.csv",
                               "mixed_screens": {k: v[0] for k, v in SCREENS.items()},
                               "separate_alone_screen": {k: v[0] for k, v in ALONE_STATE_SCREENS.items()}},
               "rows": rows,
               "controls": ctl_rows},
              open("tp53_portability.json", "w"), indent=2)

    print("Screens regressed on the common alone-state axis "
          "(urn:mavedb:00000068-b-1).\n"
          "HCT116 is shown for reference but is a mut/Delta MONOALLELIC state, not a\n"
          "mixed state, so it does not enter the portability comparison.\n")
    for scope, sub in df.groupby("scope", sort=False):
        print(f"--- {scope} ---")
        print(sub.drop(columns="scope").to_string(index=False))
        mixed_only = sub[sub.state == "mixed"]
        lo, hi = mixed_only.b.min(), mixed_only.b.max()
        iv = sorted(zip(mixed_only.ci_lo, mixed_only.ci_hi))
        ov = [f"{a}-{b}" for (a, b), (c, d) in zip(iv, iv[1:]) if b >= c]
        print(f"    mixed-state spread {lo:.3f} to {hi:.3f} = {hi/lo:.1f}-fold; "
              f"{'some intervals overlap' if ov else 'all intervals disjoint'}\n")
    print(f"common variant set: {len(common)} missense variants in both mixed screens "
          "and the Giacomelli alone axis; HCT116 is excluded from this cohort.")
    print("\nSynonymous controls (protein-level averages matched on each pair of score axes).")
    print("The through-origin slope is reported beside the CENTRED statistics,")
    print("because E[xy] = Cov(x,y) + E[x]E[y] and these axes are not centred:")
    print(ctl[["screen", "n_control", "b_control", "mean_x", "mean_y",
               "mean_product_share_of_Exy", "centred_pearson_r", "centred_p"]]
          .to_string(index=False))
    print("""
    Cross-study normalization offsets and correlated systematic effects can
    produce nonzero control slopes even without a shared physical reference.
    Neither a larger cross-study slope nor independent laboratory provenance
    establishes biological reproducibility or excludes normalization artifacts.

    The HCT116 synonymous-control slope is dominated by the product of means;
    its centered association is approximately zero. The control slopes therefore
    do not identify a shared-noise fraction or support biological deconvolution.
    The mixed-screen slope spread is descriptive and remains confounded by
    published score scales, origin choices, treatment, and assay context.
    """)

    anchors = assess_anchors()
    pd.DataFrame(anchors).to_csv("tp53_anchors.csv", index=False)
    anchor_report = {
        "generated": date.today().isoformat(),
        "n_boot": NANCHORBOOT, "seed": ANCHOR_SEED,
        "bootstrap_scope": "residues jointly across synonymous and nonsense classes, conditional on deposited scores",
        "weighting": "finite deposited scores averaged per protein substitution; one synonymous or nonsense key per residue",
        "scope_rationale": "full available coverage and coding residues 126-307 covered by Funk's scan; the latter limits late-truncation heterogeneity",
        "probability_definition": "P(nonsense score > synonymous score) + 0.5 P(tie), Mann-Whitney U divided by the number of pairs",
        "two_anchor_formula": "z = (score - median_synonymous) / (median_nonsense - median_synonymous); RFS = 2*z - 1",
        "origin_caution": "RFS=0 is the midpoint between anchors, not neutral. A neutral-centred through-origin comparison would use z=(RFS+1)/2 on both axes, conditional on appropriate anchors.",
        "mixed_two_anchor_calibration_performed": False,
        "reason": "Nonsense LOF is not a validated maximal dominant-negative endpoint in WT-background screens. The observed negative stop-minus-synonymous contrasts reinforce the distinction from mutant-alone LOF screens; sign alone is not a general test of whether affine rescaling is mathematically possible.",
        "limitations": "Two-anchor scaling removes arbitrary affine score offsets and scales conditional on valid reference classes; it does not establish equivalent biological endpoints or nonlinear response functions. This diagnostic does not re-estimate Funk's published per-exon normalization, raw-count uncertainty, or a molecular dominance coefficient.",
        "sources": {
            "Funk_2025": "https://doi.org/10.1038/s41588-024-02039-4",
            "Giacomelli_2018": "https://doi.org/10.1038/s41588-018-0204-y",
            "Boettcher_2019": "https://doi.org/10.1126/science.aax3649"},
        "rows": anchors}
    with open("tp53_anchors.json", "w") as handle:
        json.dump(anchor_report, handle, indent=2, allow_nan=False)
    print("\nCandidate synonymous/nonsense anchors (not a forced calibration):")
    print(pd.DataFrame(anchors)[["screen", "scope", "n_synonymous", "n_nonsense",
                                "synonymous_median", "nonsense_median",
                                "nonsense_minus_synonymous", "gap_ci_lo", "gap_ci_hi"]]
          .to_string(index=False))
    print("Nonsense LOF does not establish a common dominant-negative endpoint in WT-background assays.")
    print("No two-anchor calibration of the mixed-screen slopes was performed.")
    print("wrote tp53_portability.csv/json, tp53_controls.csv, tp53_anchors.csv/json")


if __name__ == "__main__":
    main()
