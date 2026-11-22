"""
dominance.py - the effect-transfer slope from paired multiplexed assays.

NAMING.  The quantity this module computes, b = sum(xy)/sum(x^2), is an
EFFECT-WEIGHTED PANEL SUMMARY of per-variant transfer ratios.  It is NOT a
per-variant dominance coefficient, not the population-genetic h, and not the
molecular coefficient g -- it is an affine transform of g under equal
expression and random dimer assembly.  Earlier versions of this module called
it a "dominance coefficient", which is the conflation the accompanying
manuscript exists to argue against, so the name is not used here.

The coefficient is the through-origin slope relating a variant's effect measured ALONE
(null background / homozygous) to its effect measured MIXED (wild-type present):

    b = sum_i(x_i * y_i) / sum_i(x_i ** 2)

Interpretation requires two things that are easy to get wrong, and this module enforces both:

1. PHYSICAL SCALE AND ORIGIN. b is a through-origin slope, so it depends on where each axis
   puts zero. Multiplying either axis by a POSITIVE constant rescales b by a constant factor
   and leaves its sign, and the ordering of per-variant ratios, intact. CENTRING an axis does
   not: separate Z-scoring moves both origins and can change the sign and the ordering
   outright -- x=(1,2), y=(2,1) gives b=+0.8 as given and b=-1.0 after each axis is Z-scored.
   Pass `units="physical"` only when both axes are in the same physical unit (mV, normalised
   current, fraction of wild-type) AND share an origin at no effect. Otherwise b is reported
   but flagged, and stoichiometric interpretation is refused.

2. SUBUNIT RATIO p. The stoichiometric reading of b depends on the mutant fraction of the
   subunit pool. Without a measured p, ladder comparison is refused rather than assumed.

Also enforced, because each caused a real error in development:

  * variant-class filtering - nonsense (Ter) and stop-loss variants are excluded from a
    "missense" panel because they are different variant classes. Their functional effects
    and their ability to inhibit wild type must be assessed in the relevant assay.
  * class-wise diagnostics - a panel with no true effect need NOT return b ~ 0, and treating
    that as a correctness check is wrong. Offsets, correlated errors, and biological class
    effects can each influence the slope. The accompanying manuscript reports benign and
    synonymous controls with large slopes but does not establish that they have zero true
    effects or identify the cause. The historical class_gates name now supplies descriptive
    diagnostics only; a class label does not define an expected slope.
  * residue clustering - substitutions at one position are not independent observations.
  * a baseline - the mixed score adjusted for the alone score captures the same signal; the
    coefficient must be reported against it, not instead of it.
  * correlation screening - a through-origin slope through an uncorrelated cloud is a ratio of
    noise. The |r| threshold below is a conservative reporting convention, not a validity
    criterion: clearing it does not make b interpretable, and no conclusion in the
    accompanying manuscript turns on its exact value.

No dependencies beyond numpy/pandas/scipy.
"""
from __future__ import annotations
import re
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

class Uninterpretable(ValueError):
    """Raised when the estimate would not be a measurement. Fails closed by design."""


AA3 = {"Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
       "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"}
_HGVS = re.compile(r"^p\.\(?([A-Za-z]{3})(\d+)([A-Za-z]{3}|=)\)?$")
# X and * are legacy one-letter spellings of a stop codon; O'Neill et al.
# write W822X, which an alphabet without X silently classified as "other".
_ONE = re.compile(r"^p?\.?([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWYX*=])$")
_1TO3 = {"A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu","G":"Gly",
         "H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser",
         "T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Ter"}


def classify_variant(hgvs_pro: str):
    """Return (wt, pos, mt, klass). klass in {missense, synonymous, nonsense, stoploss, other}."""
    s = str(hgvs_pro).strip()
    # "p.=" is the MaveDB spelling for a synonymous change whose residue is not
    # named.  It carries no position, but it is a synonymous variant, and
    # classifying it as "other" silently removed every control from the
    # zero-effect gate.
    if s in ("p.=", "p.(=)", "="):
        return (None, None, None, "synonymous")
    m = _HGVS.match(s)
    if m:
        wt, pos, mt = m.group(1), int(m.group(2)), m.group(3)
    else:
        m1 = _ONE.match(s)          # single-letter form, e.g. "R144Q" or "L81P"
        if not m1:
            return (None, None, None, "other")
        wt = _1TO3[m1.group(1)]
        pos = int(m1.group(2))
        g3 = m1.group(3)
        mt = "=" if g3 == "=" else _1TO3.get(g3, g3)
    if mt == "=" or wt == mt:
        return (wt, pos, wt, "synonymous")
    # X and * are legacy spellings of Ter; O'Neill et al. write W822X.
    if mt in ("Ter", "X", "*"):
        return (wt, pos, "Ter", "nonsense")
    if wt in ("Ter", "X", "*"):
        return ("Ter", pos, mt, "stoploss")
    if wt in AA3 and mt in AA3:
        return (wt, pos, mt, "missense")
    return (wt, pos, mt, "other")


def slope_through_origin(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = float((x * x).sum())
    if d == 0:
        return np.nan
    return float((x * y).sum() / d)


@dataclass
class Result:
    b: float
    ci: tuple
    n: int
    n_residues: int
    pearson_r: float
    units: str
    baseline_partial_r: float = np.nan
    coefficient_partial_r: float = np.nan
    interpretable: bool = True
    p_measured: bool = False
    warnings_: list = field(default_factory=list)

    def __str__(self):
        s = (f"effect-transfer slope b = {self.b:.4f}  95% CI [{self.ci[0]:.4f}, {self.ci[1]:.4f}]\n"
             f"  n = {self.n} variants at {self.n_residues} residues;  r(x,y) = {self.pearson_r:+.3f}\n"
             f"  axis units: {self.units}")
        if not np.isnan(self.coefficient_partial_r):
            s += (f"\n  vs baseline (mixed adjusted for alone): coefficient {self.coefficient_partial_r:+.3f},"
                  f" baseline {self.baseline_partial_r:+.3f}")
        for w in self.warnings_:
            s += f"\n  ! {w}"
        return s


def estimate(df, x="x", y="y", pos="pos", units="unknown", min_pearson=0.30,
             n_boot=20000, seed=20260926, target=None, allow_uncorrelated=False):
    """Estimate the effect-transfer slope with a residue-clustered bootstrap.

    Returns an effect-weighted panel summary, not a per-variant coefficient.

    df must already be restricted to one variant class (use classify_variant) and both axes
    centred so that zero means 'no effect' (e.g. on the synonymous mean).
    `target`, if given, is a per-variant binary label used for the baseline comparison.
    """
    rng = np.random.default_rng(seed)
    d = df[[c for c in (x, y, pos) if c is not None]].dropna()
    if len(d) < 20:
        raise ValueError(f"need >=20 paired observations, got {len(d)}")
    xv, yv, pv = d[x].to_numpy(float), d[y].to_numpy(float), d[pos].to_numpy()
    b = slope_through_origin(xv, yv)
    r = float(stats.pearsonr(xv, yv).statistic)

    uniq = np.unique(pv)
    idx_by_res = {u: np.flatnonzero(pv == u) for u in uniq}
    draws = np.empty(n_boot)
    for k in range(n_boot):
        pick = rng.choice(uniq, size=uniq.size, replace=True)
        idx = np.concatenate([idx_by_res[u] for u in pick])
        draws[k] = slope_through_origin(xv[idx], yv[idx])
    ci = (float(np.nanpercentile(draws, 2.5)), float(np.nanpercentile(draws, 97.5)))

    warns = []
    refused = abs(r) < min_pearson
    if refused and not allow_uncorrelated:
        raise Uninterpretable(
            f"|r(x,y)| = {abs(r):.3f} < {min_pearson}: a through-origin slope is not a "
            "measurement where the axes are uncorrelated. Pass allow_uncorrelated=True to "
            "obtain the number anyway; it will be flagged as not interpretable.")
    if refused:
        warns.append(f"NOT INTERPRETABLE: |r(x,y)| = {abs(r):.3f} < {min_pearson}.")
    if units != "physical":
        warns.append(f"axes are '{units}', not 'physical': the ABSOLUTE value of b is "
                     "scale-indeterminate. Sign and relative ordering only.")

    res = Result(b, ci, len(d), int(uniq.size), r, units, warnings_=warns,
                 interpretable=(not refused) and units == "physical")

    if target is not None:
        t = df.loc[d.index, target].to_numpy(float)
        resid = yv - b * xv
        res.coefficient_partial_r = _partial_spearman(resid, t, xv)
        res.baseline_partial_r = _partial_spearman(yv, t, xv)
        if res.baseline_partial_r >= res.coefficient_partial_r:
            warns.append("the baseline (mixed adjusted for alone) matches or beats the "
                         "coefficient here; report both.")
    return res


def _partial_spearman(a, b_, ctrl):
    A = np.column_stack([np.ones(len(a)), np.asarray(ctrl, float)])
    ra = a - A @ np.linalg.lstsq(A, np.asarray(a, float), rcond=None)[0]
    rb = b_ - A @ np.linalg.lstsq(A, np.asarray(b_, float), rcond=None)[0]
    return float(stats.spearmanr(ra, rb).statistic)


def class_gates(df, klass="klass", x="x", y="y", pos="pos"):
    """Per-class through-origin slopes, as a DIAGNOSTIC rather than a pass/fail test.

    A class with no true effect need not return b ~ 0. Score origins, correlated errors,
    and class-specific effects can influence the slope. Neither the class label nor the
    slope alone identifies those contributions. Returns a descriptive table alongside
    the panel result; it does not assign a universal expected slope to a variant class.
    """
    rows = []
    for k, g in df.groupby(klass):
        if len(g) < 20:
            rows.append(dict(klass=k, n=len(g), b=np.nan, note="too few"))
            continue
        b = slope_through_origin(g[x], g[y])
        rows.append(dict(klass=k, n=len(g), b=b,
                         note="descriptive class slope; origin and assay context required"))
    return pd.DataFrame(rows)


def stoichiometric_ladder(n_subunits: int, p: float):
    """Predicted coefficient if >=t mutant subunits inactivate an n-mer, UNDER RANDOM ASSEMBLY.

    The floor (1-p)**n is the all-wild-type complex fraction when subunits assemble at random.
    It is NOT a lower bound over all pairwise mechanisms: preferential pairing preserves the
    subunit fraction while lowering the all-WT fraction, and can reproduce observations below
    the floor with no loss of protein. Breaching the floor therefore rejects a CONJUNCTION
    (random assembly AND equal abundance AND this stoichiometry AND unchanged all-WT activity),
    and does not identify which conjunct failed.
    """
    from math import comb
    P = [comb(n_subunits, k) * p**k * (1 - p)**(n_subunits - k) for k in range(n_subunits + 1)]
    return {t: float(sum(P[k] for k in range(t, n_subunits + 1)))
            for t in range(1, n_subunits + 1)}


def compare_to_ladder(result: Result, n_subunits: int, p=None, p_measured=False,
                      assume_p_because=None):
    """Compare a coefficient to the stoichiometric ladder.

    Refuses unless the axes are physical AND p was either measured or declared
    an explicit assumption.

    p_measured=True asserts that the mutant subunit fraction was measured, and
    the caller must say what measured it.  assume_p_because="<reason>" instead
    runs the comparison as an avowed conditional: the reason is carried into
    the returned table and every row is flagged `conditional`, so a sensitivity
    analysis can be published without the result being mistaken for a
    measurement.  Supplying p with neither flag is still refused, because a
    design ratio silently promoted to a measurement is the error this guard
    exists to prevent.
    """
    if result.units != "physical":
        raise ValueError("ladder comparison refused: axes are not in physical units, so the "
                         "absolute coefficient is not comparable to a predicted fraction.")
    if p is None:
        raise Uninterpretable(
            "ladder comparison refused: the mutant subunit fraction p was not measured. "
            "Without it the model has a free parameter and any coefficient is reproduced by "
            "some (t, p) pair.")
    if not p_measured and not assume_p_because:
        raise Uninterpretable(
            "ladder comparison refused: p was supplied but not flagged as MEASURED. A design "
            "ratio (e.g. 1:1 plasmid) is an assumption about DNA, not a measurement of the "
            "subunit pool. Pass p_measured=True only with allele-specific evidence, and state "
            "what measured it; or pass assume_p_because='<reason>' to publish the comparison "
            "as an avowed conditional.")
    if p_measured and assume_p_because:
        raise ValueError("p cannot be both measured and assumed; pass one.")
    lad = stoichiometric_ladder(n_subunits, p)
    lo, hi = result.ci
    out = []
    for t, pred in lad.items():
        out.append(dict(t=t, predicted=pred, excluded_by_CI=not (lo <= pred <= hi),
                        conditional=not p_measured,
                        p_basis=assume_p_because or "measured"))
    floor = (1 - p) ** n_subunits
    return pd.DataFrame(out), floor
