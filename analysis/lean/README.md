# Machine checked proofs for the dominance manuscript

This Lean 4 project states the mathematical results of "What paired variant assays can and cannot identify about molecular dominance" as theorems and proves them. Lean's kernel checks every step of every proof.

## What is proved

- 100 public theorems in seven files under `DominanceProofs/`
- Every theorem depends only on Lean's three standard axioms: `propext`, `Classical.choice`, and `Quot.sound`
- No proof contains `sorry`, `admit`, `native_decide`, or any added axiom
- S1 Text, Section S13 and its Table O map each manuscript claim to its theorem

| File | Theorems | Content |
|:--|--:|:--|
| `Basic.lean` | 15 | Binomial weights, equations (2), (3), and (6), endpoint values, additive subunits, dimer formulas |
| `Identifiability.lean` | 7 | Vandermonde rescaling, nonsingular design, necessary and sufficient identification |
| `Bounds.lean` | 14 | Sharp bounds of Section S9, feasibility, monotonicity in the fraction, hull over a fraction scenario, the synthetic tetramer |
| `Floor.lean` | 24 | Equation (4) for any order, dimer compatibility boundaries, ratio contrasts, log ratio tipping points |
| `Arrangements.lean` | 13 | Count reduction of configuration activities, adjacent and opposite tetramer pairs, nonrandom dimer assembly |
| `Slope.lean` | 14 | Equation (1) as a weighted mean, equation (5) under its assumptions, two anchor invariance, shared offset rank counterexample |
| `Examples.lean` | 13 | Constructed examples and floor arithmetic quoted in the text |

## What is not proved

- Sampling distributions, P values, and confidence intervals
- Bootstrap, Monte Carlo, and simulation results
- Empirical values, except the arithmetic checks in `Examples.lean`
- Whether the model's assumptions hold in any cell or assay

## How to check

- Install the Lean toolchain manager `elan` (https://github.com/leanprover/elan)
- From this directory, download the prebuilt Mathlib library: `lake exe cache get`
- Build every proof: `lake build`
- Print the axioms of every public theorem: `lake env lean Audit.lean`
- The Python test `test_lean_proofs.py` in the parent directory runs the same checks and skips the build when Lean or Mathlib is unavailable

## Versions

- Lean `leanprover/lean4:v4.36.0-rc2` (see `lean-toolchain`)
- Mathlib tag `v4.36.0-rc2`, commit `9d581dfff2e41de989072ea63abf3daf69311520` (see `lake-manifest.json`)
- The Mathlib download needs network access and about 7 GB of disk; the analysis reproduction does not need it

## Authorship of the proofs

- The author's manuscript fixed what each theorem must say
- Anthropic Claude Opus 5.5 wrote the formal statements and proofs on 8 November 2026
- Each statement was written down before its proof and compared afterwards with that first version, so no proof could quietly weaken what the theorem claims
