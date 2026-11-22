/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Panel slopes, moment adjustment, anchors, and rank correlations

Main text, equations (1) and (5) and "A panel slope is an effect-weighted summary"; S1 Text,
Sections S4 and S5.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

variable {ι κ : Type*}

/-- Through-origin least-squares slope, main-text equation (1). -/
noncomputable def slope (s : Finset ι) (x y : ι → ℝ) : ℝ :=
  (∑ i ∈ s, x i * y i) / ∑ i ∈ s, x i ^ 2

/-- Equation (1): the slope is the `x²`-weighted sum of the transfer ratios `d i = y i / x i`;
variants with `x i = 0` contribute nothing to numerator or denominator. -/
theorem slope_eq_weighted (s : Finset ι) (x y : ι → ℝ) :
    slope s x y =
      ∑ i ∈ s.filter (fun i => x i ≠ 0), (x i ^ 2 / ∑ j ∈ s, x j ^ 2) * (y i / x i) := by
  unfold slope
  have hnum : ∑ i ∈ s.filter (fun i => x i ≠ 0), x i * y i = ∑ i ∈ s, x i * y i := by
    apply Finset.sum_filter_of_ne
    intro i _ h hx
    exact h (by rw [hx, zero_mul])
  rw [← hnum, Finset.sum_div]
  refine Finset.sum_congr rfl fun i hi => ?_
  have hx : x i ≠ 0 := (Finset.mem_filter.mp hi).2
  rw [div_mul_eq_mul_div]
  congr 1
  rw [div_eq_mul_inv, show x i ^ 2 * (y i * (x i)⁻¹) = x i * y i * (x i * (x i)⁻¹) by ring,
    mul_inv_cancel₀ hx, mul_one]

/-- The weights `x i² / ∑ x j²` sum to one when some `x i ≠ 0`. -/
theorem slope_weights_sum (s : Finset ι) (x : ι → ℝ) (h : ∃ i ∈ s, x i ≠ 0) :
    ∑ i ∈ s.filter (fun i => x i ≠ 0), x i ^ 2 / ∑ j ∈ s, x j ^ 2 = 1 := by
  rw [← Finset.sum_div]
  have hfil : ∑ i ∈ s.filter (fun i => x i ≠ 0), x i ^ 2 = ∑ i ∈ s, x i ^ 2 := by
    apply Finset.sum_filter_of_ne
    intro i _ h' hx
    exact h' (by rw [hx]; norm_num)
  rw [hfil]
  obtain ⟨i, hi, hxi⟩ := h
  apply div_self
  exact (Finset.sum_pos' (fun j _ => sq_nonneg (x j))
    ⟨i, hi, lt_of_le_of_ne (sq_nonneg _) (pow_ne_zero 2 hxi).symm⟩).ne'

/-- A panel slope lies between the smallest and largest transfer ratios it averages. -/
theorem slope_mem_Icc (s : Finset ι) (x y : ι → ℝ) (h : ∃ i ∈ s, x i ≠ 0) (lo hi : ℝ)
    (hr : ∀ i ∈ s, x i ≠ 0 → lo ≤ y i / x i ∧ y i / x i ≤ hi) :
    lo ≤ slope s x y ∧ slope s x y ≤ hi := by
  have hw := slope_weights_sum s x h
  have hnn : ∀ i ∈ s.filter (fun i => x i ≠ 0), 0 ≤ x i ^ 2 / ∑ j ∈ s, x j ^ 2 :=
    fun i _ => div_nonneg (sq_nonneg _) (Finset.sum_nonneg fun j _ => sq_nonneg _)
  rw [slope_eq_weighted]
  constructor
  · calc lo = ∑ i ∈ s.filter (fun i => x i ≠ 0), (x i ^ 2 / ∑ j ∈ s, x j ^ 2) * lo := by
          rw [← Finset.sum_mul, hw, one_mul]
      _ ≤ ∑ i ∈ s.filter (fun i => x i ≠ 0), (x i ^ 2 / ∑ j ∈ s, x j ^ 2) * (y i / x i) := by
          apply Finset.sum_le_sum
          intro i hmem
          have hm := Finset.mem_filter.mp hmem
          exact mul_le_mul_of_nonneg_left (hr i hm.1 hm.2).1 (hnn i hmem)
  · calc ∑ i ∈ s.filter (fun i => x i ≠ 0), (x i ^ 2 / ∑ j ∈ s, x j ^ 2) * (y i / x i)
          ≤ ∑ i ∈ s.filter (fun i => x i ≠ 0), (x i ^ 2 / ∑ j ∈ s, x j ^ 2) * hi := by
          apply Finset.sum_le_sum
          intro i hmem
          have hm := Finset.mem_filter.mp hmem
          exact mul_le_mul_of_nonneg_left (hr i hm.1 hm.2).2 (hnn i hmem)
      _ = hi := by rw [← Finset.sum_mul, hw, one_mul]

/-- Rescaling the alone-state axis by `c` divides the slope by `c`. -/
theorem slope_scale_left (s : Finset ι) (x y : ι → ℝ) (c : ℝ) (hc : c ≠ 0) :
    slope s (fun i => c * x i) y = slope s x y / c := by
  simp only [slope]
  have h1 : ∑ i ∈ s, c * x i * y i = c * ∑ i ∈ s, x i * y i := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  have h2 : ∑ i ∈ s, (c * x i) ^ 2 = c * (c * ∑ i ∈ s, x i ^ 2) := by
    rw [Finset.mul_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [h1, h2, div_div, mul_comm (∑ i ∈ s, x i ^ 2) c, mul_div_mul_left _ _ hc]

/-- Rescaling the mixed-state axis by `c` multiplies the slope by `c`. -/
theorem slope_scale_right (s : Finset ι) (x y : ι → ℝ) (c : ℝ) :
    slope s x (fun i => c * y i) = c * slope s x y := by
  simp only [slope]
  have h1 : ∑ i ∈ s, x i * (c * y i) = c * ∑ i ∈ s, x i * y i := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [h1, mul_div_assoc]

/-- Sample mean. -/
noncomputable def mean (s : Finset ι) (x : ι → ℝ) : ℝ := (∑ i ∈ s, x i) / s.card

/-- A raw second moment is the centred moment plus the product of the means, so a nonzero
through-origin slope need not indicate centred association. -/
theorem raw_moment_decomposition (s : Finset ι) (hs : s.Nonempty) (x y : ι → ℝ) :
    mean s (fun i => x i * y i) =
      mean s (fun i => (x i - mean s x) * (y i - mean s y)) + mean s x * mean s y := by
  have hN : (s.card : ℝ) ≠ 0 := by exact_mod_cast hs.card_pos.ne'
  have hexp : ∀ a b : ℝ, ∑ i ∈ s, (x i - a) * (y i - b)
      = ∑ i ∈ s, x i * y i - b * ∑ i ∈ s, x i - a * ∑ i ∈ s, y i + s.card * (a * b) := by
    intro a b
    have h1 : ∀ i ∈ s, (x i - a) * (y i - b) = x i * y i - b * x i - a * y i + a * b :=
      fun i _ => by ring
    rw [Finset.sum_congr rfl h1, Finset.sum_add_distrib, Finset.sum_sub_distrib,
      Finset.sum_sub_distrib, Finset.sum_const, nsmul_eq_mul, Finset.mul_sum, Finset.mul_sum]
  simp only [mean]
  rw [hexp]
  field_simp
  ring

/-- A concrete case: a constant response has zero centred covariance but a nonzero
through-origin slope. -/
theorem slope_without_association :
    slope (univ : Finset (Fin 2)) (fun i => (i : ℝ) + 1) (fun _ => 1) = 3 / 5 ∧
      mean (univ : Finset (Fin 2))
        (fun i => ((i : ℝ) + 1 - mean univ (fun j : Fin 2 => (j : ℝ) + 1)) *
          (1 - mean univ (fun _ : Fin 2 => (1 : ℝ)))) = 0 := by
  constructor
  · norm_num [slope, Fin.sum_univ_two]
  · norm_num [mean, Fin.sum_univ_two]

/-- Equation (5) under its stated assumptions recovers the biological through-origin slope.
Panel `P` has observed effects `X + U` and `Y + V`; controls `C` have observed effects `U'` and
`V'` (negligible biological effect). Signal–nuisance cross terms vanish and nuisance second
moments transport from controls to panel. -/
theorem moment_adjusted_recovers (P : Finset ι) (C : Finset κ)
    (X Y U V : ι → ℝ) (U' V' : κ → ℝ)
    (hXV : ∑ i ∈ P, X i * V i = 0) (hUY : ∑ i ∈ P, U i * Y i = 0)
    (hXU : ∑ i ∈ P, X i * U i = 0)
    (hUV : (∑ i ∈ P, U i * V i) / P.card = (∑ j ∈ C, U' j * V' j) / C.card)
    (hUU : (∑ i ∈ P, U i ^ 2) / P.card = (∑ j ∈ C, U' j ^ 2) / C.card) :
    ((∑ i ∈ P, (X i + U i) * (Y i + V i)) / P.card - (∑ j ∈ C, U' j * V' j) / C.card) /
        ((∑ i ∈ P, (X i + U i) ^ 2) / P.card - (∑ j ∈ C, U' j ^ 2) / C.card) =
      slope P X Y := by
  have hnum : ∑ i ∈ P, (X i + U i) * (Y i + V i)
      = ∑ i ∈ P, X i * Y i + ∑ i ∈ P, U i * V i := by
    have h1 : ∀ i ∈ P, (X i + U i) * (Y i + V i)
        = X i * Y i + X i * V i + U i * Y i + U i * V i := fun i _ => by ring
    rw [Finset.sum_congr rfl h1, Finset.sum_add_distrib, Finset.sum_add_distrib,
      Finset.sum_add_distrib, hXV, hUY, add_zero, add_zero]
  have hden : ∑ i ∈ P, (X i + U i) ^ 2 = ∑ i ∈ P, X i ^ 2 + ∑ i ∈ P, U i ^ 2 := by
    have h1 : ∀ i ∈ P, (X i + U i) ^ 2 = X i ^ 2 + 2 * (X i * U i) + U i ^ 2 :=
      fun i _ => by ring
    rw [Finset.sum_congr rfl h1, Finset.sum_add_distrib, Finset.sum_add_distrib,
      ← Finset.mul_sum, hXU, mul_zero, add_zero]
  rw [hnum, hden, add_div, add_div, ← hUV, ← hUU, add_sub_cancel_right, add_sub_cancel_right]
  unfold slope
  by_cases hP : (P.card : ℝ) = 0
  · rw [Nat.cast_eq_zero, Finset.card_eq_zero] at hP
    subst hP
    simp
  · rw [div_div_div_cancel_right₀ hP]

/-! ### Two-anchor normalization (S1 Text, Section S5) -/

/-- Two-anchor transform `z = (s - lo)/(hi - lo)`. -/
noncomputable def anchorZ (lo hi s : ℝ) : ℝ := (s - lo) / (hi - lo)

/-- Two-anchor normalization removes any affine score scaling, provided the anchors transform
with the scores (as medians do under increasing affine maps). -/
theorem anchorZ_affine (α β lo hi s : ℝ) (hα : α ≠ 0) :
    anchorZ (α * lo + β) (α * hi + β) (α * s + β) = anchorZ lo hi s := by
  unfold anchorZ
  rw [show α * s + β - (α * lo + β) = α * (s - lo) by ring,
    show α * hi + β - (α * lo + β) = α * (hi - lo) by ring, mul_div_mul_left _ _ hα]

/-- Hence anchored slopes are unchanged by separate affine rescaling of either axis. -/
theorem anchored_slope_invariant (S : Finset ι) (x y : ι → ℝ) (lx hx ly hy α β γ δ : ℝ)
    (hα : α ≠ 0) (hγ : γ ≠ 0) :
    slope S (fun i => anchorZ (α * lx + β) (α * hx + β) (α * x i + β))
        (fun i => anchorZ (γ * ly + δ) (γ * hy + δ) (γ * y i + δ)) =
      slope S (fun i => anchorZ lx hx (x i)) (fun i => anchorZ ly hy (y i)) := by
  simp only [anchorZ_affine α β _ _ _ hα, anchorZ_affine γ δ _ _ _ hγ]

/-- Relative fitness score `RFS = 2z - 1` is zero exactly at the midpoint of the anchors. -/
theorem rfs_zero_iff_midpoint (lo hi s : ℝ) (h : lo ≠ hi) :
    2 * anchorZ lo hi s - 1 = 0 ↔ s = (lo + hi) / 2 := by
  unfold anchorZ
  have h' : hi - lo ≠ 0 := sub_ne_zero.mpr (Ne.symm h)
  rw [sub_eq_zero, mul_div_assoc', div_eq_one_iff_eq h']
  constructor
  · intro H
    linarith
  · intro H
    rw [H]
    ring

/-! ### Shared offsets and rank correlations (S1 Text, Section S4) -/

/-- With true ratio one-half, adding the same offset `c > 0` to both axes gives a ratio that
decreases strictly with the alone-state effect. -/
theorem offset_ratio_strictAntiOn (c : ℝ) (hc : 0 < c) :
    StrictAntiOn (fun x : ℝ => (x / 2 + c) / (x + c)) (Set.Ioi 0) := by
  intro a ha b hb hab
  simp only [Set.mem_Ioi] at ha hb
  show (b / 2 + c) / (b + c) < (a / 2 + c) / (a + c)
  rw [div_lt_div_iff₀ (by linarith) (by linarith)]
  nlinarith [mul_pos hc (sub_pos.mpr hab)]

/-- Rank of `v i` within `s`: the number of values strictly below it. -/
noncomputable def rankIn (s : Finset ι) (v : ι → ℝ) (i : ι) : ℝ :=
  ((s.filter fun j => v j < v i).card : ℝ)

/-- Pearson correlation. -/
noncomputable def pearson (s : Finset ι) (u v : ι → ℝ) : ℝ :=
  (∑ i ∈ s, (u i - mean s u) * (v i - mean s v)) /
    Real.sqrt ((∑ i ∈ s, (u i - mean s u) ^ 2) * ∑ i ∈ s, (v i - mean s v) ^ 2)

/-- Spearman correlation: the Pearson correlation of ranks. -/
noncomputable def spearman (s : Finset ι) (u v : ι → ℝ) : ℝ :=
  pearson s (rankIn s u) (rankIn s v)

/-- Pearson correlation only depends on the values of its arguments on `s`. -/
private lemma pearson_congr_right (s : Finset ι) (u v v' : ι → ℝ) (h : ∀ i ∈ s, v i = v' i) :
    pearson s u v = pearson s u v' := by
  have hm : mean s v = mean s v' := by
    unfold mean
    rw [Finset.sum_congr rfl h]
  have h1 : ∑ i ∈ s, (u i - mean s u) * (v i - mean s v)
      = ∑ i ∈ s, (u i - mean s u) * (v' i - mean s v') :=
    Finset.sum_congr rfl fun i hi => by rw [h i hi, hm]
  have h2 : ∑ i ∈ s, (v i - mean s v) ^ 2 = ∑ i ∈ s, (v' i - mean s v') ^ 2 :=
    Finset.sum_congr rfl fun i hi => by rw [h i hi, hm]
  unfold pearson
  rw [h1, h2]

/-- A variable that is not constant on `s` has Pearson self-correlation one. -/
private lemma pearson_self (s : Finset ι) (r : ι → ℝ) {i j : ι} (hi : i ∈ s) (hj : j ∈ s)
    (hij : r i ≠ r j) : pearson s r r = 1 := by
  have hsum : ∑ k ∈ s, (r k - mean s r) * (r k - mean s r) = ∑ k ∈ s, (r k - mean s r) ^ 2 :=
    Finset.sum_congr rfl fun k _ => (sq _).symm
  have hpos : 0 < ∑ k ∈ s, (r k - mean s r) ^ 2 := by
    rcases (Finset.sum_nonneg fun k (_ : k ∈ s) => sq_nonneg (r k - mean s r)).lt_or_eq with
      hlt | heq
    · exact hlt
    · exfalso
      have h0 := (Finset.sum_eq_zero_iff_of_nonneg
        fun k (_ : k ∈ s) => sq_nonneg (r k - mean s r)).mp heq.symm
      have hi0 : r i - mean s r = 0 := (pow_eq_zero_iff two_ne_zero).mp (h0 i hi)
      have hj0 : r j - mean s r = 0 := (pow_eq_zero_iff two_ne_zero).mp (h0 j hj)
      exact hij (by linarith)
  unfold pearson
  rw [hsum, Real.sqrt_mul_self hpos.le, div_self hpos.ne']

/-- A strictly larger value has a strictly larger rank. -/
private lemma rankIn_lt (s : Finset ι) (u : ι → ℝ) {i j : ι} (hi : i ∈ s) (huij : u i < u j) :
    rankIn s u i < rankIn s u j := by
  unfold rankIn
  have hsub : s.filter (fun k => u k < u i) ⊂ s.filter (fun k => u k < u j) := by
    rw [Finset.ssubset_iff_of_subset]
    · exact ⟨i, Finset.mem_filter.mpr ⟨hi, huij⟩,
        fun hmem => lt_irrefl _ (Finset.mem_filter.mp hmem).2⟩
    · intro k hk
      rw [Finset.mem_filter] at hk ⊢
      exact ⟨hk.1, hk.2.trans huij⟩
  exact_mod_cast Finset.card_lt_card hsub

/-- Two variables that order the units identically have Spearman correlation one. -/
theorem spearman_eq_one_of_same_order (s : Finset ι) (u v : ι → ℝ)
    (h : ∀ i ∈ s, ∀ j ∈ s, (u i < u j ↔ v i < v j)) (hne : ∃ i ∈ s, ∃ j ∈ s, u i < u j) :
    spearman s u v = 1 := by
  have hrank : ∀ i ∈ s, rankIn s v i = rankIn s u i := by
    intro i hi
    have hf : s.filter (fun j => v j < v i) = s.filter (fun j => u j < u i) :=
      Finset.filter_congr fun j hj => (h j hj i hi).symm
    unfold rankIn
    rw [hf]
  obtain ⟨i, hi, j, hj, hij⟩ := hne
  unfold spearman
  rw [pearson_congr_right s (rankIn s u) (rankIn s v) (rankIn s u) hrank]
  exact pearson_self s (rankIn s u) hi hj (rankIn_lt s u hi hij).ne

/-- S1 Text, Section S4: `x i = i` for `i = 1, …, 100`, every true ratio 0.5, offsets 1 and 2 in
two endpoints. The two ratio vectors have Spearman correlation one. -/
theorem shared_offset_counterexample :
    spearman (Icc 1 100) (fun i : ℕ => ((i : ℝ) / 2 + 1) / (i + 1))
      (fun i : ℕ => ((i : ℝ) / 2 + 2) / (i + 2)) = 1 := by
  have h1 := offset_ratio_strictAntiOn 1 one_pos
  have h2 := offset_ratio_strictAntiOn 2 two_pos
  have hpos : ∀ i ∈ Icc 1 100, ((i : ℕ) : ℝ) ∈ Set.Ioi (0 : ℝ) := by
    intro i hi
    have h1i := (Finset.mem_Icc.mp hi).1
    rw [Set.mem_Ioi, Nat.cast_pos]
    omega
  apply spearman_eq_one_of_same_order
  · intro i hi j hj
    exact (h1.lt_iff_gt (hpos i hi) (hpos j hj)).trans
      (h2.lt_iff_gt (hpos i hi) (hpos j hj)).symm
  · refine ⟨2, by norm_num, 1, by norm_num, ?_⟩
    norm_num

end Dominance
