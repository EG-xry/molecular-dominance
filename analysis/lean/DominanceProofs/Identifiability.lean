/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Identification of the count-indexed profile from mixed-ratio titrations

Main text, "A mixed measurement supplies one constraint on the activity profile", and S1 Text,
Section S1. With known stoichiometry `n` and fixed endpoints `g 0 = 0`, `g n = 1`, the interior
entries `g 1, …, g (n-1)` are identified by the fractions in a design `S ⊆ [0,1]` if and only if
`S` contains at least `n - 1` distinct interior fractions. The sufficiency proof follows S1 Text:
dividing row `i` of the interior design matrix by `pᵢ(1-pᵢ)^(n-1)` and column `k` by `C(n,k)`
leaves a Vandermonde matrix in `pᵢ/(1-pᵢ)`.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

/-- Interior design matrix for `n - 1` fractions: row `i`, column `k` holds `w n (k+1) (p i)`,
the weight of the interior entry `g (k+1)`. -/
noncomputable def design (n : ℕ) (p : Fin (n - 1) → ℝ) : Matrix (Fin (n - 1)) (Fin (n - 1)) ℝ :=
  Matrix.of fun i k => w n (k + 1) (p i)

/-- The rescaling used in S1 Text: each entry factors as a row factor `pᵢ(1-pᵢ)^(n-1)`, a
Vandermonde entry `(pᵢ/(1-pᵢ))^k`, and a column factor `C(n, k+1)`. -/
theorem design_entry_rescaled (n : ℕ) (hn : 2 ≤ n) (p : Fin (n - 1) → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (i k : Fin (n - 1)) :
    design n p i k =
      (p i * (1 - p i) ^ (n - 1)) * (p i / (1 - p i)) ^ (k : ℕ) * (n.choose (k + 1) : ℝ) := by
  have hk : (k : ℕ) < n - 1 := k.isLt
  have h1 : (1 - p i) ≠ 0 := (sub_pos.mpr (hp i).2).ne'
  have hpow : (1 - p i) ^ (n - 1) = (1 - p i) ^ (n - ((k : ℕ) + 1)) * (1 - p i) ^ (k : ℕ) := by
    rw [← pow_add]
    congr 1
    omega
  simp only [design, Matrix.of_apply, w]
  rw [hpow, div_pow]
  field_simp
  ring

/-- Distinct interior fractions give a nonsingular interior design matrix. -/
theorem design_det_ne_zero (n : ℕ) (hn : 2 ≤ n) (p : Fin (n - 1) → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (hinj : Function.Injective p) :
    (design n p).det ≠ 0 := by
  have hfac : design n p = Matrix.diagonal (fun i => p i * (1 - p i) ^ (n - 1)) *
      Matrix.vandermonde (fun i => p i / (1 - p i)) *
      Matrix.diagonal (fun k : Fin (n - 1) => (n.choose ((k : ℕ) + 1) : ℝ)) := by
    ext i k
    rw [Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.vandermonde_apply,
      design_entry_rescaled n hn p hp i k]
  rw [hfac, Matrix.det_mul, Matrix.det_mul, Matrix.det_diagonal, Matrix.det_diagonal]
  refine mul_ne_zero (mul_ne_zero ?_ ?_) ?_
  · rw [Finset.prod_ne_zero_iff]
    intro i _
    exact mul_ne_zero (hp i).1.ne' (pow_ne_zero _ (sub_pos.mpr (hp i).2).ne')
  · rw [Matrix.det_vandermonde_ne_zero_iff]
    intro i j hij
    apply hinj
    have hi1 : (1 - p i) ≠ 0 := (sub_pos.mpr (hp i).2).ne'
    have hj1 : (1 - p j) ≠ 0 := (sub_pos.mpr (hp j).2).ne'
    have hij' : p i / (1 - p i) = p j / (1 - p j) := hij
    rw [div_eq_div_iff hi1 hj1] at hij'
    linear_combination hij'
  · rw [Finset.prod_ne_zero_iff]
    intro k _
    have hk : (k : ℕ) + 1 ≤ n := by have := k.isLt; omega
    exact_mod_cast (Nat.choose_pos hk).ne'

/-- Sufficiency: `n - 1` distinct interior fractions identify the interior profile. -/
theorem interior_profile_identified (n : ℕ) (hn : 2 ≤ n) (p : Fin (n - 1) → ℝ)
    (hp : ∀ i, 0 < p i ∧ p i < 1) (hinj : Function.Injective p)
    (g h : ℕ → ℝ) (hg0 : g 0 = 0) (hgn : g n = 1) (hh0 : h 0 = 0) (hhn : h n = 1)
    (heq : ∀ i, transfer n g (p i) = transfer n h (p i)) :
    ∀ k ∈ Ico 1 n, g k = h k := by
  have hmul : (design n p).mulVec
      (fun k : Fin (n - 1) => g ((k : ℕ) + 1) - h ((k : ℕ) + 1)) = 0 := by
    funext i
    have hg := transfer_sub_endpoint n (by omega) g hg0 hgn (p i)
    have hh := transfer_sub_endpoint n (by omega) h hh0 hhn (p i)
    have hdiff : ∑ k ∈ Ico 1 n, w n k (p i) * (g k - h k) = 0 := by
      have hsplit : ∑ k ∈ Ico 1 n, w n k (p i) * (g k - h k) =
          ∑ k ∈ Ico 1 n, w n k (p i) * g k - ∑ k ∈ Ico 1 n, w n k (p i) * h k := by
        rw [← Finset.sum_sub_distrib]
        refine Finset.sum_congr rfl fun k _ => ?_
        ring
      rw [hsplit, ← hg, ← hh, heq i]
      ring
    rw [Finset.sum_Ico_eq_sum_range, ← Fin.sum_univ_eq_sum_range] at hdiff
    rw [Matrix.mulVec_apply_eq_sum, Pi.zero_apply]
    refine (Finset.sum_congr rfl fun k _ => ?_).trans hdiff
    simp only [design, Matrix.of_apply, add_comm 1 (k : ℕ)]
  have hc0 := Matrix.eq_zero_of_mulVec_eq_zero (design_det_ne_zero n hn p hp hinj) hmul
  intro k hk
  rw [Finset.mem_Ico] at hk
  have hk1 := congrFun hc0 ⟨k - 1, by omega⟩
  simp only [Pi.zero_apply] at hk1
  have hk' : k - 1 + 1 = k := by omega
  rw [hk'] at hk1
  linarith

/-- Number of distinct interior fractions in a design. -/
noncomputable def interiorCount (S : Finset ℝ) : ℕ :=
  (S.filter fun p => 0 < p ∧ p < 1).card

/-- The fractions in `S` identify every interior entry of profiles with fixed endpoints. -/
def Identifies (n : ℕ) (S : Finset ℝ) : Prop :=
  ∀ g h : ℕ → ℝ, g 0 = 0 → g n = 1 → h 0 = 0 → h n = 1 →
    (∀ p ∈ S, transfer n g p = transfer n h p) → ∀ k ∈ Ico 1 n, g k = h k

/-- Main identification result (necessary and sufficient). Fractions in `[0,1]` identify the
interior profile of an `n`-subunit complex, `n ≥ 2`, if and only if they include at least
`n - 1` distinct interior fractions. -/
theorem identifies_iff (n : ℕ) (hn : 2 ≤ n) (S : Finset ℝ) (hS : ∀ p ∈ S, 0 ≤ p ∧ p ≤ 1) :
    Identifies n S ↔ n - 1 ≤ interiorCount S := by
  constructor
  · intro hid
    rw [← not_lt]
    intro hlt
    set T := S.filter (fun p => 0 < p ∧ p < 1) with hT
    have hTcard : T.card < n - 1 := hlt
    obtain ⟨M, hM⟩ : ∃ M : Matrix T (Fin (n - 1)) ℝ, ∀ q k, M q k = w n ((k : ℕ) + 1) (q : ℝ) :=
      ⟨Matrix.of fun q k => w n ((k : ℕ) + 1) (q : ℝ), fun q k => rfl⟩
    have hrank : Module.finrank ℝ (T → ℝ) < Module.finrank ℝ (Fin (n - 1) → ℝ) := by
      rw [Module.finrank_fin_fun, Module.finrank_fintype_fun_eq_card, Fintype.card_coe]
      exact hTcard
    obtain ⟨c, hcker, hc0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot
      (LinearMap.ker_ne_bot_of_finrank_lt (f := M.mulVecLin) hrank)
    rw [LinearMap.mem_ker, Matrix.mulVecLin_apply] at hcker
    obtain ⟨j, hj⟩ := Function.ne_iff.mp hc0
    have hj' : c j ≠ 0 := hj
    have hjlt := j.isLt
    obtain ⟨d, hd⟩ : ∃ d : ℕ → ℝ,
        ∀ k, d k = ∑ i : Fin (n - 1), if k = (i : ℕ) + 1 then c i else 0 :=
      ⟨_, fun k => rfl⟩
    obtain ⟨h, hh⟩ : ∃ h : ℕ → ℝ, ∀ k, h k = if k = n then 1 else 0 := ⟨_, fun k => rfl⟩
    obtain ⟨g, hg⟩ : ∃ g : ℕ → ℝ, ∀ k, g k = h k + d k := ⟨_, fun k => rfl⟩
    have hd0 : d 0 = 0 := by
      rw [hd]
      apply Finset.sum_eq_zero
      intro i _
      exact ite_eq_right (by omega)
    have hdn : d n = 0 := by
      rw [hd]
      apply Finset.sum_eq_zero
      intro i _
      have := i.isLt
      exact ite_eq_right (by omega)
    have hdj : d ((j : ℕ) + 1) = c j := by
      rw [hd, Finset.sum_eq_single j]
      · exact ite_eq_left rfl
      · intro b _ hb
        refine ite_eq_right ?_
        intro hbj
        exact hb (Fin.ext (by omega))
      · intro hjs
        exact absurd (Finset.mem_univ j) hjs
    have htd : ∀ q, transfer n d q = ∑ k : Fin (n - 1), w n ((k : ℕ) + 1) q * c k := by
      intro q
      unfold transfer
      simp only [hd, Finset.mul_sum, mul_ite, mul_zero]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [Finset.sum_ite_eq' (range (n + 1)) ((k : ℕ) + 1)]
      have := k.isLt
      rw [ite_eq_left (Finset.mem_range.mpr (by omega))]
    have hgh : ∀ q, transfer n g q = transfer n h q + transfer n d q := by
      intro q
      unfold transfer
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [hg]
      ring
    have hd_zero : ∀ q ∈ S, transfer n d q = 0 := by
      intro q hq
      by_cases hqT : 0 < q ∧ q < 1
      · have hqT' : q ∈ T := Finset.mem_filter.mpr ⟨hq, hqT⟩
        rw [htd]
        have hrow := congrFun hcker ⟨q, hqT'⟩
        rw [Matrix.mulVec_apply_eq_sum, Pi.zero_apply] at hrow
        simp only [hM] at hrow
        exact hrow
      · have hq01 := hS q hq
        have hq' : q = 0 ∨ q = 1 := by
          rcases eq_or_lt_of_le hq01.1 with h0 | h0
          · exact Or.inl h0.symm
          · rcases eq_or_lt_of_le hq01.2 with h1 | h1
            · exact Or.inr h1
            · exact absurd ⟨h0, h1⟩ hqT
        rcases hq' with rfl | rfl
        · rw [transfer_at_zero]
          exact hd0
        · rw [transfer_at_one]
          exact hdn
    have hn0 : (0 : ℕ) ≠ n := by omega
    have key := hid g h
      (by rw [hg, hh, hd0, ite_eq_right hn0, add_zero])
      (by rw [hg, hh, hdn, ite_eq_left rfl, add_zero])
      (by rw [hh, ite_eq_right hn0])
      (by rw [hh, ite_eq_left rfl])
      (fun q hq => by rw [hgh, hd_zero q hq, add_zero])
    have hk := key ((j : ℕ) + 1) (Finset.mem_Ico.mpr ⟨by omega, by omega⟩)
    rw [hg, hdj] at hk
    exact hj' (by linarith)
  · intro hcount g h hg0 hgn hh0 hhn heq
    have hTcard : n - 1 ≤ (S.filter fun p => 0 < p ∧ p < 1).card := hcount
    have hpT : ∀ i : Fin (n - 1), (S.filter fun p => 0 < p ∧ p < 1).orderEmbOfFin rfl
        (Fin.castLE hTcard i) ∈ S.filter fun p => 0 < p ∧ p < 1 :=
      fun i => Finset.orderEmbOfFin_mem _ rfl _
    refine interior_profile_identified n hn
      (fun i => (S.filter fun p => 0 < p ∧ p < 1).orderEmbOfFin rfl (Fin.castLE hTcard i))
      (fun i => (Finset.mem_filter.mp (hpT i)).2) ?_ g h hg0 hgn hh0 hhn
      (fun i => heq _ (Finset.mem_filter.mp (hpT i)).1)
    intro i j hij
    exact Fin.castLE_injective hTcard
      (((S.filter fun p => 0 < p ∧ p < 1).orderEmbOfFin rfl).injective hij)

/-- Fractions of zero or one carry no information about the interior entries. -/
theorem endpoint_fractions_uninformative (n : ℕ) (g h : ℕ → ℝ)
    (hg0 : g 0 = 0) (hgn : g n = 1) (hh0 : h 0 = 0) (hhn : h n = 1) :
    transfer n g 0 = transfer n h 0 ∧ transfer n g 1 = transfer n h 1 := by
  refine ⟨?_, ?_⟩
  · rw [transfer_at_zero, transfer_at_zero, hg0, hh0]
  · rw [transfer_at_one, transfer_at_one, hgn, hhn]

/-- For `n > 2`, one fraction leaves the profile unidentified. -/
theorem one_fraction_insufficient (n : ℕ) (hn : 3 ≤ n) (p : ℝ) (hp : 0 ≤ p ∧ p ≤ 1) :
    ¬ Identifies n {p} := by
  rw [identifies_iff n (by omega) {p}
    (fun q hq => by rw [Finset.mem_singleton.mp hq]; exact hp)]
  have h1 : interiorCount {p} ≤ 1 :=
    (Finset.card_filter_le _ _).trans (Finset.card_singleton p).le
  omega

/-- A dimer has one interior entry, so one known interior fraction suffices. -/
theorem dimer_one_fraction (p : ℝ) (hp : 0 < p ∧ p < 1) : Identifies 2 {p} := by
  rw [identifies_iff 2 le_rfl {p}
    (fun q hq => by rw [Finset.mem_singleton.mp hq]; exact ⟨hp.1.le, hp.2.le⟩)]
  have h1 : interiorCount {p} = 1 := by
    unfold interiorCount
    rw [Finset.filter_singleton, ite_eq_left hp, Finset.card_singleton]
  omega

end Dominance
