/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# The conditional assembly floor, dimer compatibility, and ratio contrasts

Main text, "SCN5A conditional assembly bound" (equation (4)) and "Unknown subunit fractions and
mutant-dose scaling have distinct consequences"; S1 Text, Sections S3 and S9.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

/-! ### Equation (4): the conditional floor -/

/-- Main-text equation (4) for any order `n`: with nonnegative composition-specific activities
and the wild-type complex at its reference activity, the mixed current is at least
`W (1-p)^n`, the all-wild-type contribution. -/
theorem assembly_floor (n : ℕ) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) (W : ℝ) (hW : 0 ≤ W)
    (a : ℕ → ℝ) (ha0 : a 0 = 1) (ha : ∀ k ∈ range (n + 1), 0 ≤ a k) :
    W * (1 - p) ^ n ≤ W * meanActivity n a p := by
  refine mul_le_mul_of_nonneg_left ?_ hW
  unfold meanActivity
  rw [Finset.sum_range_succ']
  have h0 : w n 0 p * a 0 = (1 - p) ^ n := by simp [w, ha0]
  have hrest : 0 ≤ ∑ k ∈ range n, w n (k + 1) p * a (k + 1) := by
    refine Finset.sum_nonneg fun k hk => mul_nonneg (w_nonneg hp.1 hp.2) (ha (k + 1) ?_)
    rw [Finset.mem_range] at hk ⊢
    omega
  linarith

/-- Complete poisoning (`a k = 0` for every `k ≥ 1`) attains the floor. -/
theorem complete_poisoning_attains_floor (n : ℕ) (p : ℝ) (a : ℕ → ℝ) (ha0 : a 0 = 1)
    (ha : ∀ k, 1 ≤ k → a k = 0) :
    meanActivity n a p = (1 - p) ^ n := by
  unfold meanActivity
  rw [Finset.sum_eq_single 0]
  · simp [w, ha0]
  · intro k _ hk
    rw [ha k (Nat.one_le_iff_ne_zero.mpr hk), mul_zero]
  · intro h
    exact absurd (Finset.mem_range.mpr (Nat.succ_pos n)) h

/-- Dimer: with nonnegative mutant-homodimer activity, a mixed current below the floor implies a
negative inferred heterodimer activity. -/
theorem below_floor_negative_heterodimer (p : ℝ) (hp : 0 < p ∧ p < 1) (W : ℝ) (hW : 0 < W)
    (a : ℕ → ℝ) (ha0 : a 0 = 1) (ha2 : 0 ≤ a 2)
    (h : W * meanActivity 2 a p < W * (1 - p) ^ 2) : a 1 < 0 := by
  obtain ⟨hp0, hp1⟩ := hp
  have h' : meanActivity 2 a p < (1 - p) ^ 2 := lt_of_mul_lt_mul_left h hW.le
  have hexp : meanActivity 2 a p = (1 - p) ^ 2 * a 0 + 2 * p * (1 - p) * a 1 + p ^ 2 * a 2 := by
    simp [meanActivity, w, Finset.sum_range_succ, Nat.choose]
  rw [hexp, ha0] at h'
  have hq : 0 < 1 - p := by linarith
  have hpp : 0 < 2 * p * (1 - p) := by positivity
  have hp2 : 0 ≤ p ^ 2 * a 2 := by positivity
  refine not_le.mp fun hcon => ?_
  have : 0 ≤ 2 * p * (1 - p) * a 1 := mul_nonneg hpp.le hcon
  linarith

/-- Negative heterodimer activity does not necessarily breach the floor: the mutant-homodimer
contribution can offset it. -/
theorem negative_heterodimer_need_not_breach :
    ∃ a : ℕ → ℝ, a 0 = 1 ∧ a 1 < 0 ∧ 0 ≤ a 2 ∧ (1 - 1/2 : ℝ) ^ 2 ≤ meanActivity 2 a (1/2) := by
  refine ⟨fun k => if k = 1 then -1/10 else 1, by simp, by norm_num, by norm_num, ?_⟩
  norm_num [meanActivity, w, Finset.sum_range_succ]

/-! ### `d > 1` and convex combinations -/

/-- If every composition entry is at most one, no assembly weights give `d > 1`: a ratio above one
requires some entry above one. -/
theorem gt_one_needs_entry_gt_one (n : ℕ) (π g : ℕ → ℝ) (hπ : ∀ k ∈ range (n + 1), 0 ≤ π k)
    (hsum : ∑ k ∈ range (n + 1), π k = 1) (hd : 1 < ∑ k ∈ range (n + 1), π k * g k) :
    ∃ k ∈ range (n + 1), 1 < g k := by
  by_contra hcon
  have hle : ∑ k ∈ range (n + 1), π k * g k ≤ ∑ k ∈ range (n + 1), π k := by
    refine Finset.sum_le_sum fun k hk => ?_
    have hgk : g k ≤ 1 := not_lt.mp fun h => hcon ⟨k, hk, h⟩
    calc π k * g k ≤ π k * 1 := mul_le_mul_of_nonneg_left hgk (hπ k hk)
      _ = π k := mul_one _
  linarith

/-- With `g 0 = 0` and `g n = 1`, the entry above one is interior; this holds for any assembly
weights, so redistributing assembly cannot reconcile `d > 1`. -/
theorem gt_one_needs_interior_entry (n : ℕ) (π g : ℕ → ℝ) (hπ : ∀ k ∈ range (n + 1), 0 ≤ π k)
    (hsum : ∑ k ∈ range (n + 1), π k = 1) (h0 : g 0 = 0) (hn : g n = 1)
    (hd : 1 < ∑ k ∈ range (n + 1), π k * g k) :
    ∃ k ∈ Ico 1 n, 1 < g k := by
  obtain ⟨k, hk, hgk⟩ := gt_one_needs_entry_gt_one n π g hπ hsum hd
  refine ⟨k, ?_, hgk⟩
  have hk' : k < n + 1 := Finset.mem_range.mp hk
  have hk0 : k ≠ 0 := by
    intro hk0
    rw [hk0, h0] at hgk
    linarith
  have hkn : k ≠ n := by
    intro hkn
    rw [hkn, hn] at hgk
    linarith
  rw [Finset.mem_Ico]
  omega

/-! ### Dimer compatibility boundaries (S1 Text, Section S9) -/

/-- Inferred dimer heterodimer entry from ratio `d` at fraction `p`. -/
noncomputable def g1 (p d : ℝ) : ℝ := (d - p ^ 2) / (2 * p * (1 - p))

/-- For `0 ≤ d ≤ 1`, `g₁(p)` decreases strictly on `(0, 1)`. -/
theorem g1_strictAntiOn (d : ℝ) (hd : 0 ≤ d ∧ d ≤ 1) :
    StrictAntiOn (fun p => g1 p d) (Set.Ioo 0 1) := by
  intro p hp q hq hpq
  obtain ⟨hp0, hp1⟩ := hp
  obtain ⟨hq0, hq1⟩ := hq
  obtain ⟨hd0, hd1⟩ := hd
  have hP : 0 < 2 * p * (1 - p) := mul_pos (mul_pos two_pos hp0) (by linarith)
  have hQ : 0 < 2 * q * (1 - q) := mul_pos (mul_pos two_pos hq0) (by linarith)
  have key : 0 < d * (1 - p - q) + p * q := by
    rcases le_or_gt 0 (1 - p - q) with h | h
    · have h1 : 0 ≤ d * (1 - p - q) := mul_nonneg hd0 h
      have h2 : 0 < p * q := mul_pos hp0 hq0
      linarith
    · have h1 : 0 ≤ (1 - d) * (-(1 - p - q)) := mul_nonneg (by linarith) (by linarith)
      have h2 : 0 < (1 - p) * (1 - q) := mul_pos (by linarith) (by linarith)
      linarith
  have hprod : 0 < (q - p) * (d * (1 - p - q) + p * q) := mul_pos (by linarith) key
  change (d - q ^ 2) / (2 * q * (1 - q)) < (d - p ^ 2) / (2 * p * (1 - p))
  rw [div_lt_div_iff₀ hQ hP]
  linarith

/-- The monotone hypothesis `0 ≤ g₁ ≤ 1` holds exactly for `1 - √(1-d) ≤ p ≤ √d`. -/
theorem g1_compatible_iff (d : ℝ) (hd : 0 ≤ d ∧ d ≤ 1) (p : ℝ) (hp : 0 < p ∧ p < 1) :
    (0 ≤ g1 p d ∧ g1 p d ≤ 1) ↔ (1 - Real.sqrt (1 - d) ≤ p ∧ p ≤ Real.sqrt d) := by
  obtain ⟨hd0, hd1⟩ := hd
  obtain ⟨hp0, hp1⟩ := hp
  have hq : 0 < 1 - p := by linarith
  have hden : 0 < 2 * p * (1 - p) := mul_pos (mul_pos two_pos hp0) hq
  have hA : 0 ≤ g1 p d ↔ p ≤ Real.sqrt d := by
    rw [g1, le_div_iff₀ hden, zero_mul, sub_nonneg, Real.le_sqrt hp0.le hd0]
  have hB : g1 p d ≤ 1 ↔ 1 - Real.sqrt (1 - d) ≤ p := by
    have h1 : 1 - Real.sqrt (1 - d) ≤ p ↔ (1 - p) ^ 2 ≤ 1 - d := by
      rw [sub_le_comm, Real.le_sqrt hq.le (by linarith)]
    rw [h1, g1, div_le_iff₀ hden, one_mul]
    constructor <;> intro h <;> linarith
  rw [hA, hB]
  exact and_comm

/-- The lower boundary `p⋆ = 1 - √(1-d)` satisfies `g₁ = 1`. -/
theorem g1_at_lower_boundary (d : ℝ) (hd : 0 < d ∧ d < 1) :
    g1 (1 - Real.sqrt (1 - d)) d = 1 := by
  obtain ⟨hd0, hd1⟩ := hd
  have hs2 : Real.sqrt (1 - d) ^ 2 = 1 - d := Real.sq_sqrt (by linarith)
  have hs0 : 0 < Real.sqrt (1 - d) := Real.sqrt_pos.mpr (by linarith)
  have hs1 : Real.sqrt (1 - d) < 1 := by
    rw [Real.sqrt_lt' one_pos]
    linarith
  have h1 : 0 < 1 - Real.sqrt (1 - d) := by linarith
  have h2 : 0 < 1 - (1 - Real.sqrt (1 - d)) := by linarith
  have hden : 2 * (1 - Real.sqrt (1 - d)) * (1 - (1 - Real.sqrt (1 - d))) ≠ 0 :=
    mul_ne_zero (mul_ne_zero two_ne_zero h1.ne') h2.ne'
  rw [g1, div_eq_iff hden, one_mul]
  linear_combination hs2

/-- Below the lower boundary, `g₁ > 1`. -/
theorem g1_gt_one_below (d : ℝ) (hd : 0 ≤ d ∧ d ≤ 1) (p : ℝ)
    (hp : 0 < p ∧ p < 1 - Real.sqrt (1 - d)) : 1 < g1 p d := by
  obtain ⟨hd0, hd1⟩ := hd
  obtain ⟨hp0, hp1⟩ := hp
  have hs0 : 0 ≤ Real.sqrt (1 - d) := Real.sqrt_nonneg _
  have hq0 : 0 < 1 - p := by linarith
  have hsq : 1 - d < (1 - p) ^ 2 := (Real.sqrt_lt' hq0).mp (by linarith)
  have hden : 0 < 2 * p * (1 - p) := mul_pos (mul_pos two_pos hp0) hq0
  rw [g1, one_lt_div hden]
  linarith

/-- Above the upper boundary, `g₁ < 0`. -/
theorem g1_neg_above (d : ℝ) (hd : 0 ≤ d ∧ d ≤ 1) (p : ℝ) (hp : Real.sqrt d < p ∧ p < 1) :
    g1 p d < 0 := by
  obtain ⟨hp0, hp1⟩ := hp
  have hpos : 0 < p := lt_of_le_of_lt (Real.sqrt_nonneg d) hp0
  have hsq : d < p ^ 2 := (Real.sqrt_lt hd.1 hpos.le).mp hp0
  have hden : 0 < 2 * p * (1 - p) := mul_pos (mul_pos two_pos hpos) (by linarith)
  rw [g1]
  exact div_neg_of_neg_of_pos (by linarith) hden

/-- At every interior fraction, `g₁` increases strictly with `d` (used for the uncertainty
envelopes). -/
theorem g1_strictMono_in_d (p : ℝ) (hp : 0 < p ∧ p < 1) : StrictMono (fun d => g1 p d) := by
  obtain ⟨hp0, hp1⟩ := hp
  have hden : 0 < 2 * p * (1 - p) := mul_pos (mul_pos two_pos hp0) (by linarith)
  intro d₁ d₂ h
  change (d₁ - p ^ 2) / (2 * p * (1 - p)) < (d₂ - p ^ 2) / (2 * p * (1 - p))
  exact div_lt_div_of_pos_right (by linarith) hden

/-- Over the scenario `p ∈ [0.4, 0.6]`, some monotone dimer profile is possible iff
`d ∈ [0.16, 0.84]`. -/
theorem dimer_scenario_window (d : ℝ) :
    (∃ p ∈ Set.Icc (2/5 : ℝ) (3/5), p ^ 2 ≤ d ∧ d ≤ 1 - (1 - p) ^ 2) ↔
      (4/25 ≤ d ∧ d ≤ 21/25) := by
  constructor
  · rintro ⟨p, ⟨hp1, hp2⟩, hlo, hhi⟩
    constructor
    · nlinarith [mul_self_nonneg (p - 2/5)]
    · nlinarith [mul_self_nonneg (3/5 - p)]
  · rintro ⟨hlo, hhi⟩
    rcases le_or_gt d (16/25) with h | h
    · refine ⟨2/5, ⟨le_rfl, by norm_num⟩, ?_, ?_⟩
      · norm_num
        linarith
      · norm_num
        linarith
    · refine ⟨3/5, ⟨by norm_num, le_rfl⟩, ?_, ?_⟩
      · norm_num
        linarith
      · norm_num
        linarith

/-! ### Ordering and range of normalized profiles (S1 Text, Section S9; main-text Methods) -/

/-- For `m < 1`, a larger normalized entry means a smaller activity. -/
theorem profile_order_m_lt_one (a b m : ℝ) (hm : m < 1) :
    (1 - a) / (1 - m) ≤ (1 - b) / (1 - m) ↔ b ≤ a := by
  rw [div_le_div_iff_of_pos_right (by linarith)]
  constructor <;> intro h <;> linarith

/-- For `m > 1`, the order reverses. -/
theorem profile_order_m_gt_one (a b m : ℝ) (hm : 1 < m) :
    (1 - a) / (1 - m) ≤ (1 - b) / (1 - m) ↔ a ≤ b := by
  rw [div_le_div_right_of_neg (by linarith)]
  constructor <;> intro h <;> linarith

/-- `g₁ > 1` means `a₁ < m` when `m < 1`. -/
theorem g_gt_one_iff (a₁ m : ℝ) (hm : m < 1) : 1 < (1 - a₁) / (1 - m) ↔ a₁ < m := by
  rw [one_lt_div (by linarith)]
  constructor <;> intro h <;> linarith

/-- At `m = 0.5`, the bounded model `0 ≤ g ≤ 2` is the same as `0 ≤ a ≤ 1`. -/
theorem bounded_profile_half (a : ℝ) :
    (0 ≤ (1 - a) / (1 - 1/2) ∧ (1 - a) / (1 - 1/2) ≤ 2) ↔ (0 ≤ a ∧ a ≤ 1) := by
  have h : (1 - a) / (1 - 1/2 : ℝ) = 2 * (1 - a) := by ring
  rw [h]
  constructor
  · rintro ⟨h1, h2⟩
    constructor <;> linarith
  · rintro ⟨h1, h2⟩
    constructor <;> linarith

/-- Nonnegative activity alone does not imply `g ≥ 0` at `m = 0.5`. -/
theorem nonneg_activity_not_enough : ∃ a : ℝ, 0 ≤ a ∧ (1 - a) / (1 - 1/2) < 0 := by
  exact ⟨2, by norm_num, by norm_num⟩

/-! ### Ratio contrasts (main-text Methods; S1 Text, Section S9) -/

/-- For `x = W - sA > 0`, `d = (W - H)/x > q` iff the linear contrast `(1-q)W + qsA - H` is
positive. -/
theorem ratio_gt_iff_contrast (W H A s q : ℝ) (hx : 0 < W - s * A) :
    q < (W - H) / (W - s * A) ↔ 0 < (1 - q) * W + q * s * A - H := by
  rw [lt_div_iff₀ hx]
  constructor <;> intro h <;> linarith

/-- At `q = 1`, `W` cancels: `d > 1` iff `sA - H > 0`. -/
theorem ratio_gt_one_iff (W H A s : ℝ) (hx : 0 < W - s * A) :
    1 < (W - H) / (W - s * A) ↔ 0 < s * A - H := by
  rw [one_lt_div hx]
  constructor <;> intro h <;> linarith

/-! ### Log-ratio tipping points (S1 Text, Section S3) -/

/-- The floor `W(1-p)²` and its standard error `s_W(1-p)²` have a relative standard error that
does not depend on the fraction. -/
theorem floor_relative_se (W sW p : ℝ) (hW : W ≠ 0) (hp : p ≠ 1) :
    (sW * (1 - p) ^ 2) / (W * (1 - p) ^ 2) = sW / W := by
  have hq : (1 - p) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.mpr (Ne.symm hp))
  rw [div_eq_div_iff (mul_ne_zero hW hq) hW]
  ring

/-- Under the log-ratio specification the test rejects exactly for fractions below
`p⋆ = 1 - exp((z_c σ - log(W/I))/2)`. -/
theorem logratio_tipping (W I σ zc p : ℝ) (hW : 0 < W) (hI : 0 < I) (hσ : 0 < σ) (hp : p < 1) :
    zc < Real.log (W * (1 - p) ^ 2 / I) / σ ↔
      p < 1 - Real.exp ((zc * σ - Real.log (W / I)) / 2) := by
  have hq : 0 < 1 - p := by linarith
  have hlog : Real.log (W * (1 - p) ^ 2 / I) = Real.log (W / I) + 2 * Real.log (1 - p) := by
    have h1 : W * (1 - p) ^ 2 / I = (W / I) * ((1 - p) * (1 - p)) := by ring
    rw [h1, Real.log_mul (div_pos hW hI).ne' (mul_pos hq hq).ne', Real.log_mul hq.ne' hq.ne']
    ring
  rw [lt_div_iff₀ hσ, hlog]
  constructor
  · intro h
    have h2 : (zc * σ - Real.log (W / I)) / 2 < Real.log (1 - p) := by linarith
    have h3 := (Real.lt_log_iff_exp_lt hq).mp h2
    linarith
  · intro h
    have h2 : Real.exp ((zc * σ - Real.log (W / I)) / 2) < 1 - p := by linarith
    have h3 := (Real.lt_log_iff_exp_lt hq).mpr h2
    linarith

/-- The point-estimate fraction `1 - √(I/W)` places the floor at the measured current. -/
theorem point_estimate_fraction (W I : ℝ) (hW : 0 < W) (hI : 0 ≤ I) :
    W * (1 - (1 - Real.sqrt (I / W))) ^ 2 = I := by
  have h1 : (1 - (1 - Real.sqrt (I / W))) = Real.sqrt (I / W) := by ring
  rw [h1, Real.sq_sqrt (div_nonneg hI hW.le)]
  field_simp

/-- Propagating the floor's standard error can only shrink the z statistic of a breach. -/
theorem propagated_z_le_fixed (F I sI sF : ℝ) (hbreach : 0 ≤ F - I) (hsI : 0 < sI) :
    (F - I) / Real.sqrt (sI ^ 2 + sF ^ 2) ≤ (F - I) / sI := by
  refine div_le_div_of_nonneg_left hbreach hsI ?_
  calc sI = Real.sqrt (sI ^ 2) := (Real.sqrt_sq hsI.le).symm
    _ ≤ Real.sqrt (sI ^ 2 + sF ^ 2) := Real.sqrt_le_sqrt (by linarith [sq_nonneg sF])

end Dominance
