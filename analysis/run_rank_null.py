#!/usr/bin/env python3
"""Describe KCNQ2 endpoint coupling; retire the independent-endpoint null.

Current at -40 mV depends on conductance and activation, so independent true
transfer ratios across these endpoints are not an admissible null. This
script retains the observed descriptive correlations and an explicitly
illustrative Boltzmann calculation. It emits no cross-endpoint P values.
Neither the Boltzmann parameters nor their equality across arms are fitted
or validated by this illustration.
"""
import json
from datetime import date

import numpy as np
from scipy import stats

from run_kv7_rank import CONTROL, PANEL, paired_measurements

PAIRS = [("dV1/2", "+40 mV"), ("dV1/2", "-40 mV"), ("+40 mV", "-40 mV")]


def open_probability(voltage, half_activation=-25.0, slope=10.0):
    """Illustrative Boltzmann activation, not a fit to the deposited data."""
    if slope <= 0:
        raise ValueError("slope must be positive")
    return 1.0 / (1.0 + np.exp((half_activation - np.asarray(voltage)) / slope))


def normalized_current(voltage, shift, conductance_ratio=1.0,
                       half_activation=-25.0, slope=10.0):
    return (conductance_ratio * open_probability(voltage, half_activation + shift, slope)
            / open_probability(voltage, half_activation, slope))


def main():
    panel = {k: paired_measurements(*v) for k, v in PANEL.items()}
    control = {k: paired_measurements(*v) for k, v in CONTROL.items()}
    noise = {k: float(np.sqrt(np.mean(c.x ** 2))) for k, c in control.items()}
    results = {}
    for p, q in PAIRS:
        fp = panel[p].loc[lambda d: d.x.abs() >= noise[p]].assign(r=lambda d: d.y / d.x)
        fq = panel[q].loc[lambda d: d.x.abs() >= noise[q]].assign(r=lambda d: d.y / d.x)
        j = fp.merge(fq, on="id", suffixes=("_p", "_q"), validate="one_to_one")
        results[f"{p} vs {q}"] = {
            "observed_rho": float(stats.spearmanr(j.r_p, j.r_q).statistic),
            "n_matched_filtered": len(j),
            "interpretation": "descriptive only; endpoint dependence is not evidence of transferable dominance",
        }
    out = {
        "generated": date.today().isoformat(),
        "status": "independent-endpoint null withdrawn",
        "reason": "Current at -40 mV is coupled to conductance and activation; drawing all three true endpoint ratios independently breaks the measurement relation.",
        "screening_scale": noise,
        "pairs": results,
        "boltzmann_illustration": {
            "half_activation_mV": -25.0, "slope_mV": 10.0, "shift_mV": 10.0,
            "conductance_ratio": 1.0,
            "normalized_current_minus40": float(normalized_current(-40, 10)),
            "normalized_current_plus40": float(normalized_current(40, 10)),
            "relation": "I_v/I_WT,v = (G/G_WT) P_open(v; Vhalf+shift,k)/P_open(v; Vhalf,k)",
            "scope": "Illustrates structural coupling only; not fitted predictions of the observed correlations.",
        },
        "limits": [
            "No inferential cross-endpoint null is fitted and no P values are reported.",
            "A valid residual comparison would require activation-curve parameters, their joint uncertainty, and a validated measurement model.",
            "Descriptive rank intervals condition on the measured endpoints and do not establish a biological parameter common to them.",
        ],
    }
    with open("rank_null.json", "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    for key, row in results.items():
        print(f"{key}: rho {row['observed_rho']:.3f}, n={row['n_matched_filtered']} (descriptive)")
    b = out["boltzmann_illustration"]
    print(f"Illustrative +10 mV activation shift: relative currents {b['normalized_current_minus40']:.3f} at -40 mV and {b['normalized_current_plus40']:.3f} at +40 mV.")
    print("Independent-endpoint null withdrawn; no inferential rank P values. Wrote rank_null.json.")


if __name__ == "__main__":
    main()
