/-
Copyright (c) 2026 Eric Gao. Released under the MIT license.

# Count-indexed profiles average over arrangements

Main text, "A mixed measurement supplies one constraint on the activity profile", and S1 Text,
Sections S1 and S12. Under independent random placement, the mean activity of a complex with
configuration-level activities `A` depends on `A` only through its averages over configurations
with the same mutant count. Arrangements sharing a count, such as adjacent and opposite mutant
pairs in a tetramer, are therefore not separated by any number of fractions.
-/
import DominanceProofs.Basic

open Finset

namespace Dominance

/-- Probability of one labelled configuration (the set `s` of mutant positions) under
independent placement with mutant fraction `p`. -/
noncomputable def configWeight (n : ℕ) (p : ℝ) (s : Finset (Fin n)) : ℝ :=
  p ^ s.card * (1 - p) ^ (n - s.card)

/-- Mean activity under independent placement, from configuration-level activities `A`. -/
noncomputable def configMean (n : ℕ) (A : Finset (Fin n) → ℝ) (p : ℝ) : ℝ :=
  ∑ s : Finset (Fin n), configWeight n p s * A s

/-- Average of `A` over the configurations with exactly `k` mutant subunits. -/
noncomputable def countAverage (n : ℕ) (A : Finset (Fin n) → ℝ) (k : ℕ) : ℝ :=
  (∑ s ∈ (univ : Finset (Fin n)).powersetCard k, A s) / (n.choose k : ℝ)

/-- The configuration probabilities sum to one. -/
theorem configWeight_sum (n : ℕ) (p : ℝ) : ∑ s : Finset (Fin n), configWeight n p s = 1 := by
  have h := Finset.sum_powerset_apply_card (fun m => p ^ m * (1 - p) ^ (n - m))
    (x := (univ : Finset (Fin n)))
  simp only [Finset.powerset_univ, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  unfold configWeight
  rw [h]
  refine Eq.trans ?_ (sum_w n p)
  refine Finset.sum_congr rfl fun k _ => ?_
  unfold w
  ring

/-- Count reduction: the configuration-level mean equals the binomial mixture of the
count-indexed averages, `a k = countAverage n A k` (main-text equation (3), activity form). -/
theorem configMean_eq_binomial (n : ℕ) (A : Finset (Fin n) → ℝ) (p : ℝ) :
    configMean n A p = meanActivity n (countAverage n A) p := by
  unfold configMean meanActivity countAverage
  have hsplit := Finset.sum_fiberwise_of_maps_to (s := (univ : Finset (Finset (Fin n))))
    (t := range (n + 1)) (g := Finset.card)
    (fun s _ => Finset.mem_range.mpr (Nat.lt_succ_of_le (by simpa using Finset.card_le_univ s)))
    (fun s => configWeight n p s * A s)
  rw [← hsplit]
  refine Finset.sum_congr rfl fun k hk => ?_
  have hkn : k ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
  have hc : (n.choose k : ℝ) ≠ 0 := by exact_mod_cast (Nat.choose_pos hkn).ne'
  have hfib : (univ : Finset (Finset (Fin n))).filter (fun s => s.card = k)
      = (univ : Finset (Fin n)).powersetCard k := by
    ext s
    simp
  rw [hfib]
  have hw : ∀ s ∈ (univ : Finset (Fin n)).powersetCard k,
      configWeight n p s * A s = p ^ k * (1 - p) ^ (n - k) * A s := by
    intro s hs
    unfold configWeight
    rw [Finset.mem_powersetCard_univ.mp hs]
  rw [Finset.sum_congr rfl hw, ← Finset.mul_sum]
  unfold w
  rw [mul_div_assoc', eq_div_iff hc]
  ring

/-- Configuration activities with the same count-indexed averages give the same mean activity
at every fraction, so no titration can separate them. -/
theorem configMean_congr (n : ℕ) (A B : Finset (Fin n) → ℝ)
    (h : ∀ k ≤ n, countAverage n A k = countAverage n B k) (p : ℝ) :
    configMean n A p = configMean n B p := by
  rw [configMean_eq_binomial, configMean_eq_binomial]
  unfold meanActivity
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [h k (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))]

/-! ### Tetramer: adjacent and opposite pairs -/

/-- Two positions of a tetramer (`Fin 4`, around the fourfold axis) are adjacent when they are
`{i, i + 1}` for some `i`. -/
def IsAdjacentPair (s : Finset (Fin 4)) : Prop := ∃ i : Fin 4, s = {i, i + 1}

instance : DecidablePred IsAdjacentPair := fun s => by
  unfold IsAdjacentPair; infer_instance

/-- Four of the six two-mutant configurations are adjacent pairs. -/
theorem tetramer_adjacent_count :
    (((univ : Finset (Fin 4)).powersetCard 2).filter IsAdjacentPair).card = 4 := by
  decide

/-- Two of the six are opposite pairs. -/
theorem tetramer_opposite_count :
    (((univ : Finset (Fin 4)).powersetCard 2).filter fun s => ¬ IsAdjacentPair s).card = 2 := by
  decide

/-- The two wild-type positions are adjacent exactly when the two mutant positions are. -/
theorem adjacent_iff_compl (s : Finset (Fin 4)) (hs : s.card = 2) :
    IsAdjacentPair s ↔ IsAdjacentPair sᶜ := by
  revert s
  decide

/-- Under independent placement, `a₂ = (4 a_adj + 2 a_opp) / 6`. -/
theorem tetramer_countAverage_two (A : Finset (Fin 4) → ℝ) (aAdj aOpp : ℝ)
    (hA : ∀ s : Finset (Fin 4), s.card = 2 → A s = if IsAdjacentPair s then aAdj else aOpp) :
    countAverage 4 A 2 = (4 * aAdj + 2 * aOpp) / 6 := by
  unfold countAverage
  rw [← Finset.sum_filter_add_sum_filter_not _ IsAdjacentPair]
  have h1 : ∑ s ∈ ((univ : Finset (Fin 4)).powersetCard 2).filter IsAdjacentPair, A s
      = 4 * aAdj := by
    rw [Finset.sum_congr rfl (g := fun _ => aAdj), Finset.sum_const, tetramer_adjacent_count,
      nsmul_eq_mul]
    · norm_num
    · intro s hs
      rw [Finset.mem_filter, Finset.mem_powersetCard_univ] at hs
      rw [hA s hs.1, ite_eq_left hs.2]
  have h2 : ∑ s ∈ ((univ : Finset (Fin 4)).powersetCard 2).filter (fun s => ¬ IsAdjacentPair s),
      A s = 2 * aOpp := by
    rw [Finset.sum_congr rfl (g := fun _ => aOpp), Finset.sum_const, tetramer_opposite_count,
      nsmul_eq_mul]
    · norm_num
    · intro s hs
      rw [Finset.mem_filter, Finset.mem_powersetCard_univ] at hs
      rw [hA s hs.1, ite_eq_right hs.2]
  rw [h1, h2]
  norm_num [Nat.choose]

/-- S1 Text example: `(a_adj, a_opp) = (0.25, 1)` and `(0.75, 0)` give identical mean activity at
every fraction when all other configurations agree. -/
theorem tetramer_arrangements_indistinguishable (A B : Finset (Fin 4) → ℝ)
    (hAB : ∀ s : Finset (Fin 4), s.card ≠ 2 → A s = B s)
    (hA : ∀ s : Finset (Fin 4), s.card = 2 → A s = if IsAdjacentPair s then 1/4 else 1)
    (hB : ∀ s : Finset (Fin 4), s.card = 2 → B s = if IsAdjacentPair s then 3/4 else 0)
    (p : ℝ) :
    configMean 4 A p = configMean 4 B p := by
  apply configMean_congr
  intro k _
  by_cases h2 : k = 2
  · subst h2
    rw [tetramer_countAverage_two A (1/4) 1 hA, tetramer_countAverage_two B (3/4) 0 hB]
    norm_num
  · unfold countAverage
    congr 1
    refine Finset.sum_congr rfl fun s hs => ?_
    rw [Finset.mem_powersetCard_univ] at hs
    exact hAB s (by rw [hs]; exact h2)

/-- S1 Text, Section S12: activation only with an adjacent wild-type pair gives `a₂ = 2/3`, hence
`g₂ = 1/3` when the all-mutant activity is zero. -/
theorem adjacent_only_g2 :
    (1 - countAverage 4 (fun s => if IsAdjacentPair s then 1 else 0) 2) / (1 - 0) = 1/3 := by
  rw [tetramer_countAverage_two _ 1 0 (fun s _ => rfl)]
  norm_num

/-- Activation only with an opposite wild-type pair gives `g₂ = 2/3`. -/
theorem opposite_only_g2 :
    (1 - countAverage 4 (fun s => if IsAdjacentPair s then 0 else 1) 2) / (1 - 0) = 2/3 := by
  rw [tetramer_countAverage_two _ 0 1 (fun s _ => rfl)]
  norm_num

/-! ### Non-random dimer assembly (S1 Text, Section S1) -/

/-- At overall mutant fraction one-half, dimer probabilities `(q, 1-2q, q)` are valid and keep the
marginal fraction at one-half for every `0 ≤ q ≤ 1/2`; `q = 1/4` is random assembly. -/
theorem dimer_family_marginal (q : ℝ) (hq : 0 ≤ q ∧ q ≤ 1/2) :
    0 ≤ q ∧ 0 ≤ 1 - 2 * q ∧ q + (1 - 2 * q) + q = 1 ∧ (1 - 2 * q) * (1/2) + q * 1 = 1/2 := by
  exact ⟨hq.1, by linarith [hq.2], by ring, by ring⟩

/-- With activities `(1, 0, 0.338)` the mean activity is `1.338 q`, which changes with `q` while
the fraction and composition-specific activities do not. -/
theorem dimer_family_mean (q : ℝ) : q * 1 + (1 - 2 * q) * 0 + q * 0.338 = 1.338 * q := by
  ring

/-- At `q = 0.0845` and a reference of 218.4%, the mean is about 24.7%, below the random-assembly
floor of 54.6%. -/
theorem dimer_family_crosses_floor :
    |218.4 * (1.338 * 0.0845) - 24.7| < (0.01 : ℝ) ∧
      218.4 * (1.338 * 0.0845) < (218.4 : ℝ) * (1 - 1/2) ^ 2 := by
  norm_num [abs_lt]

end Dominance
