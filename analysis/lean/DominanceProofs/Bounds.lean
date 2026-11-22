/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Bounds at known and interval-valued fractions under a monotone profile

S1 Text, Section S9, and main text, "Unknown subunit fractions and mutant-dose scaling have
distinct consequences". Under the additional restriction `0 = g 0 ≤ ⋯ ≤ g n = 1`, a single exact
ratio `d` at fraction `p` constrains each interior coordinate `g j` to a sharp interval
`[Lbound, Ubound]`; both ends are nonincreasing in `p`, so over a fraction scenario `[a, b]`
the sharp coordinatewise hull is `[Lbound(b), Ubound(a)]`.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

/-- Monotone normalized profile `0 = g 0 ≤ g 1 ≤ ⋯ ≤ g n = 1`. -/
def MonotoneProfile (n : ℕ) (g : ℕ → ℝ) : Prop :=
  g 0 = 0 ∧ g n = 1 ∧ ∀ i j, i ≤ j → j ≤ n → g i ≤ g j

/-- Smallest binomial expectation given `g j = t`: lower interior entries 0, higher ones `t`. -/
noncomputable def lowerExt (n j : ℕ) (p t : ℝ) : ℝ :=
  p ^ n + t * ∑ k ∈ Ico j n, w n k p

/-- Largest binomial expectation given `g j = t`: lower interior entries `t`, higher ones 1. -/
noncomputable def upperExt (n j : ℕ) (p t : ℝ) : ℝ :=
  ∑ k ∈ Ioc j n, w n k p + t * ∑ k ∈ Icc 1 j, w n k p

/-! ### Auxiliary lemmas -/

private lemma w_self (n : ℕ) (p : ℝ) : w n n p = p ^ n := by
  simp [w]

private lemma w_zero_right (n : ℕ) (p : ℝ) : w n 0 p = (1 - p) ^ n := by
  simp [w]

private lemma monotoneProfile_bounds {n : ℕ} {g : ℕ → ℝ} (hg : MonotoneProfile n g) {k : ℕ}
    (hk : k ≤ n) : 0 ≤ g k ∧ g k ≤ 1 := by
  obtain ⟨h0, hn, hm⟩ := hg
  have h1 := hm 0 k (Nat.zero_le k) hk
  have h2 := hm k n hk le_rfl
  constructor <;> linarith

/-- The transfer ratio is monotone in the profile (the weights are nonnegative). -/
private lemma transfer_le_transfer (n : ℕ) {f h : ℕ → ℝ} (hfh : ∀ k, k ≤ n → f k ≤ h k)
    {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : transfer n f p ≤ transfer n h p := by
  unfold transfer
  refine Finset.sum_le_sum fun k hk => ?_
  have hk' : k ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
  exact mul_le_mul_of_nonneg_left (hfh k hk') (w_nonneg hp0 hp1)

/-- The transfer ratio is linear in the profile. -/
private lemma transfer_affine (n : ℕ) (f h : ℕ → ℝ) (a b p : ℝ) :
    transfer n (fun k => a * f k + b * h k) p = a * transfer n f p + b * transfer n h p := by
  unfold transfer
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

private lemma monotoneProfile_combo {n : ℕ} {f h : ℕ → ℝ} (hf : MonotoneProfile n f)
    (hh : MonotoneProfile n h) {c : ℝ} (hc0 : 0 ≤ c) (hc1 : c ≤ 1) :
    MonotoneProfile n (fun k => (1 - c) * f k + c * h k) := by
  obtain ⟨hf0, hfn, hfm⟩ := hf
  obtain ⟨hh0, hhn, hhm⟩ := hh
  refine ⟨by simp [hf0, hh0], by simp [hfn, hhn], fun i j hij hjn => ?_⟩
  have h1 := mul_le_mul_of_nonneg_left (hfm i j hij hjn) (by linarith : (0 : ℝ) ≤ 1 - c)
  have h2 := mul_le_mul_of_nonneg_left (hhm i j hij hjn) hc0
  simp only
  linarith

/-- Any value between the transfer ratios of two monotone profiles is the transfer ratio of a
convex combination of them. -/
private lemma convex_attain (n : ℕ) (p : ℝ) {f h : ℕ → ℝ} (hf : MonotoneProfile n f)
    (hh : MonotoneProfile n h) (d : ℝ) (hd : transfer n f p ≤ d ∧ d ≤ transfer n h p) :
    ∃ c : ℝ, 0 ≤ c ∧ c ≤ 1 ∧ MonotoneProfile n (fun k => (1 - c) * f k + c * h k) ∧
      transfer n (fun k => (1 - c) * f k + c * h k) p = d := by
  rcases eq_or_lt_of_le (hd.1.trans hd.2) with heq | hlt
  · refine ⟨0, le_rfl, zero_le_one, monotoneProfile_combo hf hh le_rfl zero_le_one, ?_⟩
    rw [transfer_affine]
    have : transfer n f p = d := by linarith [hd.1, hd.2]
    simp [this]
  · have hpos : 0 < transfer n h p - transfer n f p := sub_pos.mpr hlt
    have hc0 : 0 ≤ (d - transfer n f p) / (transfer n h p - transfer n f p) :=
      div_nonneg (sub_nonneg.mpr hd.1) hpos.le
    have hc1 : (d - transfer n f p) / (transfer n h p - transfer n f p) ≤ 1 :=
      (div_le_one hpos).mpr (by linarith [hd.2])
    refine ⟨_, hc0, hc1, monotoneProfile_combo hf hh hc0 hc1, ?_⟩
    rw [transfer_affine]
    have hne : transfer n h p - transfer n f p ≠ 0 := hpos.ne'
    field_simp
    ring

/-- Lower extremal profile: `0` below `j`, `t` on `j, …, n-1`, `1` from `n` on. -/
private noncomputable def gLow (n j : ℕ) (t : ℝ) (k : ℕ) : ℝ :=
  if k < j then 0 else if k < n then t else 1

/-- Upper extremal profile: `0` at `0`, `t` on `1, …, j`, `1` above `j`. -/
private noncomputable def gUp (j : ℕ) (t : ℝ) (k : ℕ) : ℝ :=
  if k = 0 then 0 else if k ≤ j then t else 1

private lemma gLow_profile {n j : ℕ} (hj1 : 1 ≤ j) (hjn : j ≤ n) {t : ℝ} (ht0 : 0 ≤ t)
    (ht1 : t ≤ 1) : MonotoneProfile n (gLow n j t) := by
  refine ⟨?_, ?_, fun a b hab hbn => ?_⟩
  · simp [gLow, show 0 < j by omega]
  · simp [gLow, show ¬ n < j by omega]
  · unfold gLow
    split_ifs <;> linarith

private lemma gUp_profile {n j : ℕ} (hjn : j < n) {t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    MonotoneProfile n (gUp j t) := by
  refine ⟨?_, ?_, fun a b hab hbn => ?_⟩
  · simp [gUp]
  · simp [gUp, show n ≠ 0 by omega, show ¬ n ≤ j by omega]
  · unfold gUp
    split_ifs <;> first | linarith | (exfalso; omega)

private lemma transfer_gLow {n j : ℕ} (hjn : j ≤ n) (t p : ℝ) :
    transfer n (gLow n j t) p = lowerExt n j p t := by
  unfold transfer lowerExt
  rw [Finset.sum_range_succ, ← Finset.sum_range_add_sum_Ico _ hjn]
  have h1 : ∑ k ∈ range j, w n k p * gLow n j t k = 0 :=
    Finset.sum_eq_zero fun k hk => by
      simp [gLow, Finset.mem_range.mp hk]
  have h2 : ∑ k ∈ Ico j n, w n k p * gLow n j t k = t * ∑ k ∈ Ico j n, w n k p := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k hk => ?_
    obtain ⟨hk1, hk2⟩ := Finset.mem_Ico.mp hk
    rw [show gLow n j t k = t by simp [gLow, show ¬ k < j by omega, hk2]]
    ring
  have h3 : gLow n j t n = 1 := by simp [gLow, show ¬ n < j by omega]
  rw [h1, h2, h3, w_self]
  ring

private lemma sum_range_split_up (n j : ℕ) (hjn : j ≤ n) (f : ℕ → ℝ) :
    ∑ k ∈ range (n + 1), f k = f 0 + ∑ k ∈ Icc 1 j, f k + ∑ k ∈ Ioc j n, f k := by
  have hs : range (n + 1) = insert 0 (Icc 1 j ∪ Ioc j n) := by
    ext k
    simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_union, Finset.mem_Icc,
      Finset.mem_Ioc]
    omega
  have hd : Disjoint (Icc 1 j) (Ioc j n) := by
    rw [Finset.disjoint_left]
    intro k hk1 hk2
    simp only [Finset.mem_Icc, Finset.mem_Ioc] at hk1 hk2
    omega
  have h0 : 0 ∉ Icc 1 j ∪ Ioc j n := by
    simp only [Finset.mem_union, Finset.mem_Icc, Finset.mem_Ioc, not_or]
    omega
  rw [hs, Finset.sum_insert h0, Finset.sum_union hd, add_assoc]

private lemma transfer_gUp {n j : ℕ} (hjn : j ≤ n) (t p : ℝ) :
    transfer n (gUp j t) p = upperExt n j p t := by
  unfold transfer upperExt
  rw [sum_range_split_up n j hjn]
  have h1 : ∑ k ∈ Icc 1 j, w n k p * gUp j t k = t * ∑ k ∈ Icc 1 j, w n k p := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k hk => ?_
    obtain ⟨hk1, hk2⟩ := Finset.mem_Icc.mp hk
    rw [show gUp j t k = t by simp [gUp, show k ≠ 0 by omega, hk2]]
    ring
  have h2 : ∑ k ∈ Ioc j n, w n k p * gUp j t k = ∑ k ∈ Ioc j n, w n k p := by
    refine Finset.sum_congr rfl fun k hk => ?_
    obtain ⟨hk1, hk2⟩ := Finset.mem_Ioc.mp hk
    rw [show gUp j t k = 1 by simp [gUp, show k ≠ 0 by omega, show ¬ k ≤ j by omega]]
    ring
  rw [h1, h2, show gUp j t 0 = 0 by simp [gUp]]
  ring

private lemma gLow_le {n j : ℕ} {g : ℕ → ℝ} (hg : MonotoneProfile n g) :
    ∀ k, k ≤ n → gLow n j (g j) k ≤ g k := by
  obtain ⟨h0, hn, hm⟩ := hg
  intro k hk
  unfold gLow
  split_ifs with h1 h2
  · have := hm 0 k (Nat.zero_le _) hk
    linarith
  · exact hm j k (by omega) hk
  · have : k = n := by omega
    simp [this, hn]

private lemma le_gUp {n j : ℕ} (hjn : j ≤ n) {g : ℕ → ℝ} (hg : MonotoneProfile n g) :
    ∀ k, k ≤ n → g k ≤ gUp j (g j) k := by
  obtain ⟨h0, hn, hm⟩ := hg
  intro k hk
  unfold gUp
  split_ifs with h1 h2
  · simp [h1, h0]
  · exact hm k j h2 hjn
  · have := hm k n hk le_rfl
    linarith

/-- The extremal expectations bound every monotone profile with `g j = t`. -/
theorem ext_bounds (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1)
    (g : ℕ → ℝ) (hg : MonotoneProfile n g) :
    lowerExt n j p (g j) ≤ transfer n g p ∧ transfer n g p ≤ upperExt n j p (g j) := by
  constructor
  · calc lowerExt n j p (g j) = transfer n (gLow n j (g j)) p := (transfer_gLow hj.2.le _ _).symm
      _ ≤ transfer n g p := transfer_le_transfer n (gLow_le hg) hp.1 hp.2
  · calc transfer n g p ≤ transfer n (gUp j (g j)) p :=
          transfer_le_transfer n (le_gUp hj.2.le hg) hp.1 hp.2
      _ = upperExt n j p (g j) := transfer_gUp hj.2.le _ _

/-- Every expectation between the extrema is attained by a monotone profile with `g j = t`
(the feasible profile space at fixed `t` is convex). -/
theorem ext_attained (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1)
    (t : ℝ) (ht : 0 ≤ t ∧ t ≤ 1) (d : ℝ) (hd : lowerExt n j p t ≤ d ∧ d ≤ upperExt n j p t) :
    ∃ g, MonotoneProfile n g ∧ g j = t ∧ transfer n g p = d := by
  -- The range hypothesis on `p` is not needed: the transfer ratio is linear in the profile.
  obtain ⟨-, -⟩ := hp
  have hL := gLow_profile hj.1 hj.2.le ht.1 ht.2
  have hU := gUp_profile hj.2 ht.1 ht.2
  obtain ⟨c, -, -, hmono, htr⟩ := convex_attain n p hL hU d
    ⟨by rw [transfer_gLow hj.2.le]; exact hd.1, by rw [transfer_gUp hj.2.le]; exact hd.2⟩
  refine ⟨_, hmono, ?_, htr⟩
  have e1 : gLow n j t j = t := by simp [gLow, hj.2]
  have e2 : gUp j t j = t := by simp [gUp, show j ≠ 0 by omega]
  show (1 - c) * gLow n j t j + c * gUp j t j = t
  rw [e1, e2]
  ring

/-- A monotone profile is compatible with an exact ratio iff `p^n ≤ d ≤ 1 - (1-p)^n`. -/
theorem feasible_iff (n : ℕ) (hn : 1 ≤ n) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) (d : ℝ) :
    (∃ g, MonotoneProfile n g ∧ transfer n g p = d) ↔ (p ^ n ≤ d ∧ d ≤ 1 - (1 - p) ^ n) := by
  have hlow : transfer n (gLow n n 0) p = p ^ n := by
    rw [transfer_gLow le_rfl]
    simp [lowerExt]
  have hup : transfer n (gUp 0 0) p = 1 - (1 - p) ^ n := by
    rw [transfer_gUp (Nat.zero_le n)]
    have h := sum_range_split_up n 0 (Nat.zero_le n) (fun k => w n k p)
    rw [sum_w, w_zero_right] at h
    simp only [upperExt]
    simp only [Finset.Icc_eq_empty_of_lt (show 0 < 1 by norm_num), Finset.sum_empty] at h ⊢
    linarith
  have hPL := gLow_profile (n := n) (j := n) (t := 0) hn le_rfl le_rfl zero_le_one
  have hPU := gUp_profile (n := n) (j := 0) (t := 0) (by omega) le_rfl zero_le_one
  constructor
  · rintro ⟨g, hg, rfl⟩
    constructor
    · rw [← hlow]
      refine transfer_le_transfer n (fun k hk => ?_) hp.1 hp.2
      have hb := monotoneProfile_bounds hg hk
      unfold gLow
      split_ifs with h1
      · exact hb.1
      · have hkn : k = n := by omega
        rw [hkn, hg.2.1]
    · rw [← hup]
      refine transfer_le_transfer n (fun k hk => ?_) hp.1 hp.2
      have hb := monotoneProfile_bounds hg hk
      unfold gUp
      split_ifs with h1 h2
      · rw [h1, hg.1]
      · have hk0 : k = 0 := by omega
        rw [hk0, hg.1]
      · exact hb.2
  · rintro ⟨h1, h2⟩
    obtain ⟨c, -, -, hmono, htr⟩ := convex_attain n p hPL hPU d
      ⟨by rw [hlow]; exact h1, by rw [hup]; exact h2⟩
    exact ⟨_, hmono, htr⟩

/-- An observed ratio outside `[0,1]` has no compatible fraction under the monotone model. -/
theorem no_fraction_outside_unit (n : ℕ) (hn : 1 ≤ n) (d : ℝ) (hd : d < 0 ∨ 1 < d)
    (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) :
    ¬ ∃ g, MonotoneProfile n g ∧ transfer n g p = d := by
  intro hex
  have h := (feasible_iff n hn p hp d).mp hex
  have h1 : 0 ≤ p ^ n := pow_nonneg hp.1 n
  have h2 : 0 ≤ (1 - p) ^ n := pow_nonneg (by linarith [hp.2]) n
  rcases hd with hd | hd <;> linarith [h.1, h.2]

/-- Sharp lower marginal bound `L_j(p, d)`. -/
noncomputable def Lbound (n j : ℕ) (p d : ℝ) : ℝ :=
  max 0 ((d - ∑ k ∈ Ioc j n, w n k p) / ∑ k ∈ Icc 1 j, w n k p)

/-- Sharp upper marginal bound `U_j(p, d)`. -/
noncomputable def Ubound (n j : ℕ) (p d : ℝ) : ℝ :=
  min 1 ((d - p ^ n) / ∑ k ∈ Ico j n, w n k p)

private lemma sum_Icc_pos {n j : ℕ} (hj1 : 1 ≤ j) (hjn : j ≤ n) {p : ℝ} (hp0 : 0 < p)
    (hp1 : p < 1) : 0 < ∑ k ∈ Icc 1 j, w n k p :=
  Finset.sum_pos (fun k hk => w_pos (by have := (Finset.mem_Icc.mp hk).2; omega) hp0 hp1)
    ⟨1, Finset.mem_Icc.mpr ⟨le_rfl, hj1⟩⟩

private lemma sum_Ico_pos {n j : ℕ} (hjn : j < n) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    0 < ∑ k ∈ Ico j n, w n k p :=
  Finset.sum_pos (fun k hk => w_pos (by have := (Finset.mem_Ico.mp hk).2; omega) hp0 hp1)
    ⟨j, Finset.mem_Ico.mpr ⟨le_rfl, hjn⟩⟩

/-- Sharp marginal bounds: under feasibility, `g j` ranges exactly over `[Lbound, Ubound]`. -/
theorem sharp_marginal_bounds (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (p : ℝ) (hp : 0 < p ∧ p < 1)
    (d : ℝ) (hd : p ^ n ≤ d ∧ d ≤ 1 - (1 - p) ^ n) (t : ℝ) :
    (∃ g, MonotoneProfile n g ∧ transfer n g p = d ∧ g j = t) ↔
      (Lbound n j p d ≤ t ∧ t ≤ Ubound n j p d) := by
  -- The feasibility hypothesis is not needed: both sides fail together when `d` is infeasible.
  obtain ⟨-, -⟩ := hd
  have hB := sum_Icc_pos (n := n) hj.1 hj.2.le hp.1 hp.2
  have hC := sum_Ico_pos (n := n) hj.2 hp.1 hp.2
  unfold Lbound Ubound
  constructor
  · rintro ⟨g, hg, hT, rfl⟩
    have hb := ext_bounds n j hj p ⟨hp.1.le, hp.2.le⟩ g hg
    rw [hT] at hb
    unfold lowerExt upperExt at hb
    have h01 := monotoneProfile_bounds hg hj.2.le
    refine ⟨max_le h01.1 ?_, le_min h01.2 ?_⟩
    · rw [div_le_iff₀ hB]
      linarith [hb.2]
    · rw [le_div_iff₀ hC]
      linarith [hb.1]
  · rintro ⟨hL, hU⟩
    obtain ⟨ht0, hL'⟩ := max_le_iff.mp hL
    obtain ⟨ht1, hU'⟩ := le_min_iff.mp hU
    rw [div_le_iff₀ hB] at hL'
    rw [le_div_iff₀ hC] at hU'
    obtain ⟨g, hg, hgj, hT⟩ := ext_attained n j hj p ⟨hp.1.le, hp.2.le⟩ t ⟨ht0, ht1⟩ d
      ⟨by unfold lowerExt; linarith, by unfold upperExt; linarith⟩
    exact ⟨g, hg, hT, hgj⟩

private lemma w_succ_succ (n k : ℕ) (p : ℝ) :
    w (n + 1) (k + 1) p = p * w n k p + (1 - p) * w n (k + 1) p := by
  unfold w
  rw [Nat.choose_succ_succ', show n + 1 - (k + 1) = n - k by omega]
  push_cast
  rcases Nat.lt_or_ge k n with hk | hk
  · rw [show n - k = n - (k + 1) + 1 by omega]
    ring
  · rw [Nat.choose_eq_zero_of_lt (show n < k + 1 by omega), Nat.cast_zero]
    ring

private lemma w_succ_zero (n : ℕ) (p : ℝ) : w (n + 1) 0 p = (1 - p) * w n 0 p := by
  rw [w_zero_right, w_zero_right]
  ring

private lemma w_succ_self (n : ℕ) (p : ℝ) : w n (n + 1) p = 0 := by
  simp [w, Nat.choose_succ_self]

/-- One-step recursion of the binomial expectation (Pascal's rule). -/
private lemma transfer_succ (n : ℕ) (g : ℕ → ℝ) (p : ℝ) :
    transfer (n + 1) g p = p * transfer n (fun k => g (k + 1)) p + (1 - p) * transfer n g p := by
  unfold transfer
  have hA : ∑ k ∈ range (n + 1), w n (k + 1) p * g (k + 1) + w n 0 p * g 0 =
      ∑ k ∈ range (n + 1), w n k p * g k := by
    rw [Finset.sum_range_succ (fun k => w n (k + 1) p * g (k + 1)), w_succ_self,
      Finset.sum_range_succ' (fun k => w n k p * g k)]
    ring
  rw [Finset.sum_range_succ' _ (n + 1), ← hA]
  simp only [w_succ_succ, w_succ_zero, add_mul, Finset.sum_add_distrib, mul_assoc,
    ← Finset.mul_sum]
  ring

/-- Stochastic monotonicity: the binomial expectation of a nondecreasing profile is
nondecreasing in the mutant fraction. -/
theorem transfer_monotoneOn (n : ℕ) (g : ℕ → ℝ) (hg : ∀ i j, i ≤ j → j ≤ n → g i ≤ g j) :
    MonotoneOn (transfer n g) (Set.Icc 0 1) := by
  induction n generalizing g with
  | zero =>
    intro p _ q _ _
    simp [transfer, w]
  | succ n ih =>
    intro p hp q hq hpq
    have hA := ih (fun k => g (k + 1))
      (fun i j hij hjn => hg (i + 1) (j + 1) (by omega) (by omega)) hp hq hpq
    have hB := ih g (fun i j hij hjn => hg i j hij (by omega)) hp hq hpq
    have hC : transfer n g p ≤ transfer n (fun k => g (k + 1)) p :=
      transfer_le_transfer n (fun k hk => hg k (k + 1) (by omega) (by omega)) hp.1 hp.2
    rw [transfer_succ, transfer_succ]
    have h1 := mul_le_mul_of_nonneg_left hA hq.1
    have h2 := mul_le_mul_of_nonneg_left hB (sub_nonneg.mpr hq.2)
    have h3 := mul_le_mul_of_nonneg_left hC (sub_nonneg.mpr hpq)
    linarith

/-- Interior fractions compatible with `d` under the monotone model. -/
def Feasible (n : ℕ) (d p : ℝ) : Prop :=
  0 < p ∧ p < 1 ∧ p ^ n ≤ d ∧ d ≤ 1 - (1 - p) ^ n

private lemma upperExt_one {n j : ℕ} (hjn : j ≤ n) (p : ℝ) :
    upperExt n j p 1 = 1 - (1 - p) ^ n := by
  have h := sum_range_split_up n j hjn (fun k => w n k p)
  rw [sum_w, w_zero_right] at h
  unfold upperExt
  linarith

/-- The clipped lower bound is nonincreasing on the feasible fraction set. -/
theorem Lbound_antitoneOn (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (d : ℝ) :
    AntitoneOn (fun p => Lbound n j p d) {p | Feasible n d p} := by
  intro p hp q hq hpq
  obtain ⟨hp0, hp1, -, hpd2⟩ := hp
  obtain ⟨hq0, hq1, -, -⟩ := hq
  show Lbound n j q d ≤ Lbound n j p d
  have hBp := sum_Icc_pos (n := n) hj.1 hj.2.le hp0 hp1
  have hBq := sum_Icc_pos (n := n) hj.1 hj.2.le hq0 hq1
  have hs0 : 0 ≤ Lbound n j p d := le_max_left _ _
  have hs1 : Lbound n j p d ≤ 1 := by
    refine max_le zero_le_one ?_
    rw [div_le_one hBp]
    have := upperExt_one hj.2.le p
    unfold upperExt at this
    linarith
  have hsd : d ≤ upperExt n j p (Lbound n j p d) := by
    have h : (d - ∑ k ∈ Ioc j n, w n k p) / ∑ k ∈ Icc 1 j, w n k p ≤ Lbound n j p d :=
      le_max_right _ _
    rw [div_le_iff₀ hBp] at h
    unfold upperExt
    linarith
  have hmono : upperExt n j p (Lbound n j p d) ≤ upperExt n j q (Lbound n j p d) := by
    rw [← transfer_gUp hj.2.le (Lbound n j p d) p, ← transfer_gUp hj.2.le (Lbound n j p d) q]
    exact transfer_monotoneOn n _ (gUp_profile hj.2 hs0 hs1).2.2 ⟨hp0.le, hp1.le⟩
      ⟨hq0.le, hq1.le⟩ hpq
  show max 0 ((d - ∑ k ∈ Ioc j n, w n k q) / ∑ k ∈ Icc 1 j, w n k q) ≤ Lbound n j p d
  refine max_le hs0 ?_
  rw [div_le_iff₀ hBq]
  unfold upperExt at hmono hsd
  linarith

/-- The clipped upper bound is nonincreasing on the feasible fraction set. -/
theorem Ubound_antitoneOn (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (d : ℝ) :
    AntitoneOn (fun p => Ubound n j p d) {p | Feasible n d p} := by
  intro p hp q hq hpq
  obtain ⟨hp0, hp1, hpd1, -⟩ := hp
  obtain ⟨hq0, hq1, -, -⟩ := hq
  show Ubound n j q d ≤ Ubound n j p d
  have hCp := sum_Ico_pos (n := n) hj.2 hp0 hp1
  have hCq := sum_Ico_pos (n := n) hj.2 hq0 hq1
  by_cases hr : 1 ≤ (d - p ^ n) / ∑ k ∈ Ico j n, w n k p
  · have hU1 : Ubound n j p d = 1 := min_eq_left hr
    rw [hU1]
    exact min_le_left _ _
  · replace hr := not_le.mp hr
    have hUp : Ubound n j p d = (d - p ^ n) / ∑ k ∈ Ico j n, w n k p := min_eq_right hr.le
    have hr0 : 0 ≤ (d - p ^ n) / ∑ k ∈ Ico j n, w n k p := div_nonneg (by linarith) hCp.le
    have hrd : lowerExt n j p ((d - p ^ n) / ∑ k ∈ Ico j n, w n k p) = d := by
      unfold lowerExt
      rw [div_mul_cancel₀ _ hCp.ne']
      ring
    have hmono : lowerExt n j p ((d - p ^ n) / ∑ k ∈ Ico j n, w n k p) ≤
        lowerExt n j q ((d - p ^ n) / ∑ k ∈ Ico j n, w n k p) := by
      rw [← transfer_gLow hj.2.le _ p, ← transfer_gLow hj.2.le _ q]
      exact transfer_monotoneOn n _ (gLow_profile hj.1 hj.2.le hr0 hr.le).2.2 ⟨hp0.le, hp1.le⟩
        ⟨hq0.le, hq1.le⟩ hpq
    rw [hUp]
    show min 1 ((d - q ^ n) / ∑ k ∈ Ico j n, w n k q) ≤ _
    refine (min_le_right _ _).trans ?_
    rw [div_le_iff₀ hCq]
    rw [hrd] at hmono
    unfold lowerExt at hmono
    linarith

/-- Over a feasible fraction scenario `[a, b]` the sharp coordinatewise hull is
`[Lbound n j b d, Ubound n j a d]`: every compatible profile at every fraction in the scenario
lies inside it, and both ends are attained. -/
theorem hull_over_interval (n j : ℕ) (hj : 1 ≤ j ∧ j < n) (d a b : ℝ) (hab : a ≤ b)
    (ha : Feasible n d a) (hb : Feasible n d b) :
    (∀ p ∈ Set.Icc a b, ∀ g, MonotoneProfile n g → transfer n g p = d →
        Lbound n j b d ≤ g j ∧ g j ≤ Ubound n j a d) ∧
    (∃ g, MonotoneProfile n g ∧ transfer n g b = d ∧ g j = Lbound n j b d) ∧
    (∃ g, MonotoneProfile n g ∧ transfer n g a = d ∧ g j = Ubound n j a d) := by
  -- `hab` is not needed: for an empty scenario the first clause is vacuous.
  have _ := hab
  have hfeas : ∀ p ∈ Set.Icc a b, Feasible n d p := by
    rintro p ⟨hap, hpb⟩
    obtain ⟨ha0, ha1, had1, had2⟩ := ha
    obtain ⟨hb0, hb1, hbd1, hbd2⟩ := hb
    refine ⟨by linarith, by linarith, ?_, ?_⟩
    · exact (pow_le_pow_left₀ (by linarith) hpb n).trans hbd1
    · have : (1 - p) ^ n ≤ (1 - a) ^ n := pow_le_pow_left₀ (by linarith) (by linarith) n
      linarith
  have hLU : ∀ p, Feasible n d p → Lbound n j p d ≤ Ubound n j p d := by
    rintro p ⟨hp0, hp1, hpd1, hpd2⟩
    obtain ⟨g, hg, hT⟩ := (feasible_iff n (by omega) p ⟨hp0.le, hp1.le⟩ d).mpr ⟨hpd1, hpd2⟩
    have h := (sharp_marginal_bounds n j hj p ⟨hp0, hp1⟩ d ⟨hpd1, hpd2⟩ (g j)).mp
      ⟨g, hg, hT, rfl⟩
    exact h.1.trans h.2
  refine ⟨?_, ?_, ?_⟩
  · intro p hp g hg hT
    have hpf := hfeas p hp
    obtain ⟨hp0, hp1, hpd1, hpd2⟩ := hpf
    have h := (sharp_marginal_bounds n j hj p ⟨hp0, hp1⟩ d ⟨hpd1, hpd2⟩ (g j)).mp
      ⟨g, hg, hT, rfl⟩
    have hL : Lbound n j b d ≤ Lbound n j p d := Lbound_antitoneOn n j hj d (hfeas p hp) hb hp.2
    have hU : Ubound n j p d ≤ Ubound n j a d := Ubound_antitoneOn n j hj d ha (hfeas p hp) hp.1
    exact ⟨hL.trans h.1, h.2.trans hU⟩
  · obtain ⟨hb0, hb1, hbd1, hbd2⟩ := hb
    exact (sharp_marginal_bounds n j hj b ⟨hb0, hb1⟩ d ⟨hbd1, hbd2⟩ (Lbound n j b d)).mpr
      ⟨le_rfl, hLU b ⟨hb0, hb1, hbd1, hbd2⟩⟩
  · obtain ⟨ha0, ha1, had1, had2⟩ := ha
    exact (sharp_marginal_bounds n j hj a ⟨ha0, ha1⟩ d ⟨had1, had2⟩ (Ubound n j a d)).mpr
      ⟨hLU a ⟨ha0, ha1, had1, had2⟩, le_rfl⟩

/-- The compatible fractions form the interval `[1 - (1-d)^(1/n), d^(1/n)]`. -/
theorem feasible_fraction_iff (n : ℕ) (hn : 1 ≤ n) (d : ℝ) (hd : 0 ≤ d ∧ d ≤ 1)
    (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) :
    (p ^ n ≤ d ∧ d ≤ 1 - (1 - p) ^ n) ↔
      (1 - (1 - d) ^ ((n : ℝ)⁻¹) ≤ p ∧ p ≤ d ^ ((n : ℝ)⁻¹)) := by
  have hnpos : (0 : ℝ) < n := by exact_mod_cast (show 0 < n by omega)
  have key : ∀ x y : ℝ, 0 ≤ x → 0 ≤ y → (x ^ n ≤ y ↔ x ≤ y ^ ((n : ℝ)⁻¹)) := by
    intro x y hx hy
    rw [Real.le_rpow_inv_iff_of_pos hx hy hnpos, Real.rpow_natCast]
  have h1 := key p d hp.1 hd.1
  have h2 := key (1 - p) (1 - d) (by linarith [hp.2]) (by linarith [hd.2])
  constructor
  · rintro ⟨hl, hu⟩
    refine ⟨?_, h1.mp hl⟩
    have := h2.mp (by linarith)
    linarith
  · rintro ⟨hl, hu⟩
    refine ⟨h1.mpr hu, ?_⟩
    have := h2.mpr (by linarith)
    linarith

/-- A ratio of zero requires `p = 0`, and a ratio of one requires `p = 1`. -/
theorem endpoint_ratios_force_pure (n : ℕ) (hn : 1 ≤ n) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) :
    (p ^ n ≤ 0 → p = 0) ∧ (1 ≤ 1 - (1 - p) ^ n → p = 1) := by
  have hn0 : n ≠ 0 := by omega
  constructor
  · intro h
    have h0 : p ^ n = 0 := le_antisymm h (pow_nonneg hp.1 n)
    exact (pow_eq_zero_iff hn0).mp h0
  · intro h
    have h0 : (1 - p) ^ n = 0 := le_antisymm (by linarith) (pow_nonneg (by linarith [hp.2]) n)
    have := (pow_eq_zero_iff hn0).mp h0
    linarith

/-! ### The synthetic tetramer of S1 Text, Section S9 (`d = 0.5`) -/

/-- At the known fraction `p = 0.5`: `g₁ ∈ [0, 0.5]`, `g₂ ∈ [0.3, 0.7]`, `g₃ ∈ [0.5, 1]`. -/
theorem tetramer_example_point :
    Lbound 4 1 (1/2) (1/2) = 0 ∧ Ubound 4 1 (1/2) (1/2) = 1/2 ∧
    Lbound 4 2 (1/2) (1/2) = 3/10 ∧ Ubound 4 2 (1/2) (1/2) = 7/10 ∧
    Lbound 4 3 (1/2) (1/2) = 1/2 ∧ Ubound 4 3 (1/2) (1/2) = 1 := by
  simp only [Lbound, Ubound, Finset.Ioc_ofNat_ofNat, Finset.Icc_ofNat_ofNat,
    Finset.Ico_ofNat_ofNat, Finset.sum_singleton]
  norm_num [w, Nat.choose]

/-- The scenario `[0.4, 0.6]` is feasible throughout for `d = 0.5`. -/
theorem tetramer_example_feasible :
    Feasible 4 (1/2) (2/5) ∧ Feasible 4 (1/2) (3/5) := by
  unfold Feasible
  norm_num

/-- Over `[0.4, 0.6]` the hulls are `[0, 0.5616]`, `[0.0497, 0.9503]`, `[0.4384, 1]`
(four-decimal rounding stated as half-open intervals). -/
theorem tetramer_example_scenario :
    Lbound 4 1 (3/5) (1/2) = 0 ∧
    (0.56155 ≤ Ubound 4 1 (2/5) (1/2) ∧ Ubound 4 1 (2/5) (1/2) < 0.56165) ∧
    (0.04965 ≤ Lbound 4 2 (3/5) (1/2) ∧ Lbound 4 2 (3/5) (1/2) < 0.04975) ∧
    (0.95025 ≤ Ubound 4 2 (2/5) (1/2) ∧ Ubound 4 2 (2/5) (1/2) < 0.95035) ∧
    (0.43835 ≤ Lbound 4 3 (3/5) (1/2) ∧ Lbound 4 3 (3/5) (1/2) < 0.43845) ∧
    Ubound 4 3 (2/5) (1/2) = 1 := by
  simp only [Lbound, Ubound, Finset.Ioc_ofNat_ofNat, Finset.Icc_ofNat_ofNat,
    Finset.Ico_ofNat_ofNat, Finset.sum_singleton]
  norm_num [w, Nat.choose]

end Dominance
