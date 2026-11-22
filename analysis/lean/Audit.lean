import DominanceProofs

/-! Axiom audit: every public theorem should depend only on
`propext`, `Classical.choice`, and `Quot.sound`. Run with
`lake env lean Audit.lean`. -/

-- Basic.lean (15 theorems)
#print axioms Dominance.w_nonneg
#print axioms Dominance.w_pos
#print axioms Dominance.sum_w
#print axioms Dominance.w_at_zero
#print axioms Dominance.w_at_one
#print axioms Dominance.transfer_at_zero
#print axioms Dominance.transfer_at_one
#print axioms Dominance.normalize_zero
#print axioms Dominance.normalize_self
#print axioms Dominance.transfer_normalize
#print axioms Dominance.transfer_sub_endpoint
#print axioms Dominance.transfer_linear
#print axioms Dominance.dimer_transfer
#print axioms Dominance.dimer_half
#print axioms Dominance.dimer_inversion

-- Identifiability.lean (7 theorems)
#print axioms Dominance.design_entry_rescaled
#print axioms Dominance.design_det_ne_zero
#print axioms Dominance.interior_profile_identified
#print axioms Dominance.identifies_iff
#print axioms Dominance.endpoint_fractions_uninformative
#print axioms Dominance.one_fraction_insufficient
#print axioms Dominance.dimer_one_fraction

-- Bounds.lean (14 theorems)
#print axioms Dominance.ext_bounds
#print axioms Dominance.ext_attained
#print axioms Dominance.feasible_iff
#print axioms Dominance.no_fraction_outside_unit
#print axioms Dominance.sharp_marginal_bounds
#print axioms Dominance.transfer_monotoneOn
#print axioms Dominance.Lbound_antitoneOn
#print axioms Dominance.Ubound_antitoneOn
#print axioms Dominance.hull_over_interval
#print axioms Dominance.feasible_fraction_iff
#print axioms Dominance.endpoint_ratios_force_pure
#print axioms Dominance.tetramer_example_point
#print axioms Dominance.tetramer_example_feasible
#print axioms Dominance.tetramer_example_scenario

-- Floor.lean (24 theorems)
#print axioms Dominance.assembly_floor
#print axioms Dominance.complete_poisoning_attains_floor
#print axioms Dominance.below_floor_negative_heterodimer
#print axioms Dominance.negative_heterodimer_need_not_breach
#print axioms Dominance.gt_one_needs_entry_gt_one
#print axioms Dominance.gt_one_needs_interior_entry
#print axioms Dominance.g1_strictAntiOn
#print axioms Dominance.g1_compatible_iff
#print axioms Dominance.g1_at_lower_boundary
#print axioms Dominance.g1_gt_one_below
#print axioms Dominance.g1_neg_above
#print axioms Dominance.g1_strictMono_in_d
#print axioms Dominance.dimer_scenario_window
#print axioms Dominance.profile_order_m_lt_one
#print axioms Dominance.profile_order_m_gt_one
#print axioms Dominance.g_gt_one_iff
#print axioms Dominance.bounded_profile_half
#print axioms Dominance.nonneg_activity_not_enough
#print axioms Dominance.ratio_gt_iff_contrast
#print axioms Dominance.ratio_gt_one_iff
#print axioms Dominance.floor_relative_se
#print axioms Dominance.logratio_tipping
#print axioms Dominance.point_estimate_fraction
#print axioms Dominance.propagated_z_le_fixed

-- Arrangements.lean (13 theorems)
#print axioms Dominance.configWeight_sum
#print axioms Dominance.configMean_eq_binomial
#print axioms Dominance.configMean_congr
#print axioms Dominance.tetramer_adjacent_count
#print axioms Dominance.tetramer_opposite_count
#print axioms Dominance.adjacent_iff_compl
#print axioms Dominance.tetramer_countAverage_two
#print axioms Dominance.tetramer_arrangements_indistinguishable
#print axioms Dominance.adjacent_only_g2
#print axioms Dominance.opposite_only_g2
#print axioms Dominance.dimer_family_marginal
#print axioms Dominance.dimer_family_mean
#print axioms Dominance.dimer_family_crosses_floor

-- Slope.lean (14 theorems)
#print axioms Dominance.slope_eq_weighted
#print axioms Dominance.slope_weights_sum
#print axioms Dominance.slope_mem_Icc
#print axioms Dominance.slope_scale_left
#print axioms Dominance.slope_scale_right
#print axioms Dominance.raw_moment_decomposition
#print axioms Dominance.slope_without_association
#print axioms Dominance.moment_adjusted_recovers
#print axioms Dominance.anchorZ_affine
#print axioms Dominance.anchored_slope_invariant
#print axioms Dominance.rfs_zero_iff_midpoint
#print axioms Dominance.offset_ratio_strictAntiOn
#print axioms Dominance.spearman_eq_one_of_same_order
#print axioms Dominance.shared_offset_counterexample

-- Examples.lean (13 theorems)
#print axioms Dominance.nonmonotone_profile
#print axioms Dominance.nonmonotone_ratios
#print axioms Dominance.unconstrained_recovers
#print axioms Dominance.monotone_fit_minimum
#print axioms Dominance.monotone_fit_rmse
#print axioms Dominance.monotonicity_example
#print axioms Dominance.pal_mapping
#print axioms Dominance.falconer_mapping
#print axioms Dominance.floor_at_half
#print axioms Dominance.scn5a_dimer_floor_row
#print axioms Dominance.other_order_floors
#print axioms Dominance.reference_excess
#print axioms Dominance.d785n_contrast
