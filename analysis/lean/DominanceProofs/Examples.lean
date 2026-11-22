/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Worked examples and arithmetic quoted in the manuscript

Exact checks of the constructed examples in S1 Text, Section S1 and the main text, and of the
floor arithmetic in the main-text Methods and Table 2.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

/-! ### A nonmonotone profile with activities in `[0, 1]` (S1 Text, Section S1) -/

/-- Activities `a = (1, 0.2, 0.9, 0.5)` for `n = 3`. -/
noncomputable def exampleActivity : ℕ → ℝ := fun k =>
  if k = 0 then 1 else if k = 1 then 0.2 else if k = 2 then 0.9 else 0.5

/-- The normalized profile is `g = (0, 1.6, 0.2, 1)`. -/
theorem nonmonotone_profile :
    normalize 3 exampleActivity 0 = 0 ∧ normalize 3 exampleActivity 1 = 1.6 ∧
      normalize 3 exampleActivity 2 = 0.2 ∧ normalize 3 exampleActivity 3 = 1 := by
  norm_num [normalize, exampleActivity]

/-- Exact transfer ratios at fractions 0.1 and 0.9 are 0.3952 and 0.8208. -/
theorem nonmonotone_ratios :
    transfer 3 (normalize 3 exampleActivity) 0.1 = 0.3952 ∧
      transfer 3 (normalize 3 exampleActivity) 0.9 = 0.8208 := by
  norm_num [transfer, w, normalize, exampleActivity, Finset.sum_range_succ, Nat.choose]

/-- Profile with interior entries `g₁, g₂` and endpoints `0` and `1`, for `n = 3`. -/
noncomputable def profile3 (g₁ g₂ : ℝ) : ℕ → ℝ := fun k =>
  if k = 0 then 0 else if k = 1 then g₁ else if k = 2 then g₂ else 1

/-- Least-squares objective for the two exact ratios. -/
noncomputable def monotoneObjective (g₁ g₂ : ℝ) : ℝ :=
  (transfer 3 (profile3 g₁ g₂) 0.1 - 0.3952) ^ 2 + (transfer 3 (profile3 g₁ g₂) 0.9 - 0.8208) ^ 2

/-- Unconstrained inversion recovers the true interior entries exactly. -/
theorem unconstrained_recovers : monotoneObjective 1.6 0.2 = 0 := by
  norm_num [monotoneObjective, transfer, w, profile3, Finset.sum_range_succ, Nat.choose]

/-- The monotone fit over `0 ≤ g₁ ≤ g₂ ≤ 1` is `(0.9, 0.9)`. -/
theorem monotone_fit_minimum (g₁ g₂ : ℝ) (_h0 : 0 ≤ g₁) (h12 : g₁ ≤ g₂) (_h1 : g₂ ≤ 1) :
    monotoneObjective 0.9 0.9 ≤ monotoneObjective g₁ g₂ := by
  norm_num [monotoneObjective, transfer, w, profile3, Finset.sum_range_succ, Nat.choose]
  nlinarith [sq_nonneg (0.243 * (g₁ - 0.9) + 0.027 * (g₂ - 0.9)),
    sq_nonneg (0.027 * (g₁ - 0.9) + 0.243 * (g₂ - 0.9)), h12]

/-- Its profile RMSE is 0.700 even without noise. -/
theorem monotone_fit_rmse : Real.sqrt (((1.6 - 0.9) ^ 2 + (0.2 - 0.9) ^ 2) / 2 : ℝ) = 0.7 := by
  rw [show (((1.6 - 0.9) ^ 2 + (0.2 - 0.9) ^ 2) / 2 : ℝ) = 0.7 ^ 2 by norm_num]
  exact Real.sqrt_sq (by norm_num)

/-! ### Other constructed examples -/

/-- Main text: `a = (1, 0, 0.5)` gives `g = (0, 2, 1)`, a valid nonnegative profile that is not
monotone. -/
theorem monotonicity_example :
    (1 - (1 : ℝ)) / (1 - 0.5) = 0 ∧ (1 - (0 : ℝ)) / (1 - 0.5) = 2 ∧
      (1 - (0.5 : ℝ)) / (1 - 0.5) = 1 := by
  norm_num

/-- Pal et al.'s homodimer normalization `m = 1 - r`, `a₁ = 1 - g r` gives `g₁ = g`. -/
theorem pal_mapping (r gP : ℝ) (hr : r ≠ 0) : (1 - (1 - gP * r)) / (1 - (1 - r)) = gP := by
  rw [sub_sub_cancel, sub_sub_cancel, mul_div_assoc, div_self hr, mul_one]

/-- Quantitative-genetic notation: homozygotes `M + a` and `M - a` and heterozygote `M + δ`
give `x = 2a`, `y = a - δ`, and `d = 1/2 - δ/(2a)`. -/
theorem falconer_mapping (M a δ : ℝ) (ha : a ≠ 0) :
    ((M + a) - (M + δ)) / ((M + a) - (M - a)) = 1 / 2 - δ / (2 * a) := by
  have h2a : (2 * a) ≠ 0 := mul_ne_zero two_ne_zero ha
  rw [show (M + a) - (M - a) = 2 * a by ring, show (M + a) - (M + δ) = a - δ by ring, sub_div,
    show a / (2 * a) = 1 / 2 by rw [div_eq_div_iff h2a two_ne_zero]; ring]

/-! ### Floor arithmetic (main-text Methods and Table 2) -/

/-- At `p = 0.5` the reference `218.4 ± 7.7` gives `F = 54.6 ± 1.925`. -/
theorem floor_at_half :
    (218.4 : ℝ) * (1 - 0.5) ^ 2 = 54.6 ∧ (7.7 : ℝ) * (1 - 0.5) ^ 2 = 1.925 := by
  norm_num

/-- The dimer-floor row of S1 Table B (SCN5A), to the displayed precision. -/
theorem scn5a_dimer_floor_row :
    |(218.4 : ℝ) * (1 - 0.40) ^ 2 - 78.6| < 0.05 ∧ |(218.4 : ℝ) * (1 - 0.45) ^ 2 - 66.1| < 0.05 ∧
      |(218.4 : ℝ) * (1 - 0.48) ^ 2 - 59.1| < 0.05 ∧ |(218.4 : ℝ) * (1 - 0.52) ^ 2 - 50.3| < 0.05 ∧
      |(218.4 : ℝ) * (1 - 0.55) ^ 2 - 44.2| < 0.05 ∧ |(218.4 : ℝ) * (1 - 0.60) ^ 2 - 34.9| < 0.05 := by
  norm_num [abs_lt]

/-- Trimer, tetramer, and additive-reference floors at `p = 0.5`. -/
theorem other_order_floors :
    (218.4 : ℝ) * (1 - 0.5) ^ 3 = 27.3 ∧ (218.4 : ℝ) * (1 - 0.5) ^ 4 = 13.65 ∧
      (200 : ℝ) * (1 - 0.5) ^ 2 = 50 := by
  norm_num

/-- The two-wild-type reference exceeds 200% by 18.4 points, 9.2% of that expectation. -/
theorem reference_excess : (218.4 : ℝ) - 200 = 18.4 ∧ (18.4 : ℝ) / 200 = 0.092 := by
  norm_num

/-- p.Asp785Asn: mixed `24.7 ± 5.6` versus alone `36.9 ± 6.7` gives `12.2 ± 8.7`. -/
theorem d785n_contrast :
    (36.9 : ℝ) - 24.7 = 12.2 ∧ |Real.sqrt ((5.6 : ℝ) ^ 2 + 6.7 ^ 2) - 8.7| < 0.05 := by
  refine ⟨by norm_num, ?_⟩
  have hlo : (8.65 : ℝ) < Real.sqrt ((5.6 : ℝ) ^ 2 + 6.7 ^ 2) := by
    rw [Real.lt_sqrt (by norm_num)]
    norm_num
  have hhi : Real.sqrt ((5.6 : ℝ) ^ 2 + 6.7 ^ 2) < 8.75 := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  rw [abs_lt]
  constructor <;> linarith

end Dominance
