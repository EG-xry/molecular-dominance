/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Composition-dependent activity model: definitions and basic identities

This file fixes the objects used throughout the formal verification of
"What paired variant assays can and cannot identify about molecular dominance".

* `w n k p` is the binomial assembly weight: the probability that a complex of `n`
  subunits, drawn independently from a pool with mutant fraction `p`, contains exactly
  `k` mutant subunits.
* `transfer n g p` is the transfer ratio `d(p) = ∑ₖ w n k p * g k` of main-text equation (3).
* `normalize n a` turns composition-specific activities `a k` (with `a 0 = 1`, `a n = m ≠ 1`)
  into the normalized activity-effect profile `g k = (1 - a k) / (1 - m)` of equation (2).

Profiles are functions `ℕ → ℝ`; only the values at `0, …, n` matter.
-/
import Mathlib

open Finset

namespace Dominance

/-- Binomial assembly weight `C(n,k) p^k (1-p)^(n-k)`. -/
noncomputable def w (n k : ℕ) (p : ℝ) : ℝ :=
  (n.choose k : ℝ) * p ^ k * (1 - p) ^ (n - k)

/-- Transfer ratio implied at mutant fraction `p` by a count-indexed profile `g`
(main-text equation (3), right-hand side). -/
noncomputable def transfer (n : ℕ) (g : ℕ → ℝ) (p : ℝ) : ℝ :=
  ∑ k ∈ range (n + 1), w n k p * g k

/-- Normalized activity-effect profile `g k = (1 - a k) / (1 - a n)` (main-text equation (2)). -/
noncomputable def normalize (n : ℕ) (a : ℕ → ℝ) : ℕ → ℝ :=
  fun k => (1 - a k) / (1 - a n)

/-- Mean activity of the assembled population, `∑ₖ w n k p * a k`. -/
noncomputable def meanActivity (n : ℕ) (a : ℕ → ℝ) (p : ℝ) : ℝ :=
  ∑ k ∈ range (n + 1), w n k p * a k

theorem w_nonneg {n k : ℕ} {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : 0 ≤ w n k p := by
  unfold w
  have h1 : 0 ≤ 1 - p := by linarith
  positivity

theorem w_pos {n k : ℕ} {p : ℝ} (hk : k ≤ n) (hp0 : 0 < p) (hp1 : p < 1) : 0 < w n k p := by
  unfold w
  have h1 : 0 < 1 - p := by linarith
  have hc : (0 : ℝ) < n.choose k := by exact_mod_cast Nat.choose_pos hk
  positivity

/-- The binomial weights sum to one (binomial theorem). -/
theorem sum_w (n : ℕ) (p : ℝ) : ∑ k ∈ range (n + 1), w n k p = 1 := by
  have h := (add_pow p (1 - p) n).symm
  rw [show p + (1 - p) = 1 by ring, one_pow] at h
  rw [← h]
  refine Finset.sum_congr rfl fun k _ => ?_
  unfold w
  ring

theorem w_at_zero (n k : ℕ) : w n k 0 = if k = 0 then 1 else 0 := by
  unfold w
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp
  · simp [hk.ne', zero_pow hk.ne']

theorem w_at_one (n k : ℕ) (hk : k ≤ n) : w n k 1 = if k = n then 1 else 0 := by
  unfold w
  by_cases h : k = n
  · subst h; simp
  · have : n - k ≠ 0 := by omega
    simp [h, zero_pow this]

/-- At `p = 0` the transfer ratio equals the wild-type entry. -/
theorem transfer_at_zero (n : ℕ) (g : ℕ → ℝ) : transfer n g 0 = g 0 := by
  unfold transfer
  rw [Finset.sum_eq_single 0]
  · simp [w_at_zero]
  · intro k _ hk
    simp [w_at_zero, hk]
  · intro h
    exact absurd (Finset.mem_range.mpr (Nat.succ_pos n)) h

/-- At `p = 1` the transfer ratio equals the all-mutant entry. -/
theorem transfer_at_one (n : ℕ) (g : ℕ → ℝ) : transfer n g 1 = g n := by
  unfold transfer
  rw [Finset.sum_eq_single n]
  · simp [w_at_one n n le_rfl]
  · intro k hk hkn
    have hk' : k ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
    simp [w_at_one n k hk', hkn]
  · intro h
    exact absurd (Finset.mem_range.mpr (Nat.lt_succ_self n)) h

theorem normalize_zero (n : ℕ) (a : ℕ → ℝ) (h0 : a 0 = 1) : normalize n a 0 = 0 := by
  simp [normalize, h0]

theorem normalize_self (n : ℕ) (a : ℕ → ℝ) (hm : a n ≠ 1) : normalize n a n = 1 := by
  have : 1 - a n ≠ 0 := sub_ne_zero.mpr (Ne.symm hm)
  simp [normalize, div_self this]

/-- Main-text equations (2)–(3): the loss of mean activity, divided by the loss of
all-mutant activity, is the binomial expectation of the normalized profile. -/
theorem transfer_normalize (n : ℕ) (a : ℕ → ℝ) (p : ℝ) (hm : a n ≠ 1) :
    (1 - meanActivity n a p) / (1 - a n) = transfer n (normalize n a) p := by
  have hden : 1 - a n ≠ 0 := sub_ne_zero.mpr (Ne.symm hm)
  unfold meanActivity transfer normalize
  rw [eq_comm, ← sub_eq_zero]
  have key : ∑ k ∈ range (n + 1), w n k p * ((1 - a k) / (1 - a n))
      = (∑ k ∈ range (n + 1), w n k p - ∑ k ∈ range (n + 1), w n k p * a k) / (1 - a n) := by
    rw [← Finset.sum_sub_distrib, Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  rw [key, sum_w]
  ring

/-- Main-text equation (6): with `g 0 = 0` and `g n = 1`, subtracting the all-mutant
contribution `p^n` leaves the interior terms. -/
theorem transfer_sub_endpoint (n : ℕ) (hn : 1 ≤ n) (g : ℕ → ℝ) (h0 : g 0 = 0) (h1 : g n = 1)
    (p : ℝ) : transfer n g p - p ^ n = ∑ k ∈ Ico 1 n, w n k p * g k := by
  unfold transfer
  have hsplit : range (n + 1) = insert 0 (insert n (Ico 1 n)) := by
    ext k
    simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_Ico]
    omega
  rw [hsplit, Finset.sum_insert, Finset.sum_insert]
  · simp [h0, h1, w]
  · simp
  · simp only [Finset.mem_insert, Finset.mem_Ico, not_or]
    omega

/-- Independent (additive) subunits: the linear profile `g k = k / n` gives `d(p) = p`. -/
theorem transfer_linear (n : ℕ) (hn : 1 ≤ n) (p : ℝ) :
    transfer n (fun k => (k : ℝ) / n) p = p := by
  have hsum := congrArg (Polynomial.eval p) (bernsteinPolynomial.sum_smul (R := ℝ) n)
  simp only [Polynomial.eval_finsetSum, nsmul_eq_mul,
    Polynomial.eval_X, bernsteinPolynomial, Polynomial.eval_mul, Polynomial.eval_pow,
    Polynomial.eval_sub, Polynomial.eval_one, Polynomial.eval_natCast] at hsum
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast (show n ≠ 0 by omega)
  unfold transfer w
  calc ∑ k ∈ range (n + 1), (n.choose k : ℝ) * p ^ k * (1 - p) ^ (n - k) * ((k : ℝ) / n)
      = (∑ k ∈ range (n + 1), (k : ℝ) * ((n.choose k : ℝ) * p ^ k * (1 - p) ^ (n - k))) / n := by
        rw [Finset.sum_div]
        refine Finset.sum_congr rfl fun k _ => ?_
        ring
    _ = (n * p) / n := by rw [hsum]
    _ = p := by field_simp

/-- Dimer form of equation (3). -/
theorem dimer_transfer (g : ℕ → ℝ) (h0 : g 0 = 0) (h2 : g 2 = 1) (p : ℝ) :
    transfer 2 g p = 2 * p * (1 - p) * g 1 + p ^ 2 := by
  simp [transfer, w, Finset.sum_range_succ, h0, h2, Nat.choose]

/-- Dimer at equal fractions: `d(1/2) = 1/4 + g₁/2`. -/
theorem dimer_half (g : ℕ → ℝ) (h0 : g 0 = 0) (h2 : g 2 = 1) :
    transfer 2 g (1 / 2) = 1 / 4 + g 1 / 2 := by
  rw [dimer_transfer g h0 h2]
  ring

/-- Dimer inversion at a known interior fraction: `g₁ = (d - p²) / (2p(1-p))`. -/
theorem dimer_inversion (g : ℕ → ℝ) (h0 : g 0 = 0) (h2 : g 2 = 1) {p : ℝ}
    (hp0 : 0 < p) (hp1 : p < 1) :
    g 1 = (transfer 2 g p - p ^ 2) / (2 * p * (1 - p)) := by
  have hden : 2 * p * (1 - p) ≠ 0 := by
    have : 0 < 1 - p := by linarith
    positivity
  rw [dimer_transfer g h0 h2, eq_div_iff hden]
  ring

end Dominance
