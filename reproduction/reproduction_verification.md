# Isolated reproduction verification — 2026-11-22

**Passed.** All requested commands ran from the actual S1 File archive in a fresh system temporary directory. Every manifest entry was checked before extraction. Excluded third-party inputs were supplied explicitly and verified against INPUT_RIGHTS.json before analysis. No undeclared source-tree fallback was used. Network connections were disabled during analysis and recorded if attempted. This verifies reproduction with the declared external inputs; it does not establish that the ZIP is a complete offline input bundle or that historical remote data remain retrievable.

The unchanged full profile simulation was deliberately not rerun. Its packaged outputs were retained as inputs to the bounds cross-check and excluded from reproduced-artifact counts.

- Archive: `S1_File.zip`
- Archive SHA-256: `873011ae39149d411f71ce299262ece3b258cf7a26f142070c8bcb32c793fe73`
- Manifest files verified: 244
- Isolated analysis directory: `<temporary directory>/archive/analysis`
- Interpreter: `python3`
- Total elapsed time: 32.61 s
- Network attempts: 0
- Matching artifact comparisons: 63/63
- CSV/JSON numerical tolerance: relative and absolute 1e-12.
- JSON generation dates and absolute package prefixes were normalized; scientific metadata and numerical fields were retained. PNG pixels and PDF bytes were compared, excluding PDF creation/modification dates.

## Commands and exit status

| Driver | Exit status | Elapsed seconds | Log |
| --- | --- | --- | --- |
| extract_scn5a | 0 | 0.32 | [extract_scn5a](reproduction_logs/extract_scn5a.log) |
| run_deconvolution | 0 | 6.72 | [run_deconvolution](reproduction_logs/run_deconvolution.log) |
| run_kv7_controls | 0 | 1.84 | [run_kv7_controls](reproduction_logs/run_kv7_controls.log) |
| run_tp53 | 0 | 8.84 | [run_tp53](reproduction_logs/run_tp53.log) |
| run_kv7 | 0 | 1.46 | [run_kv7](reproduction_logs/run_kv7.log) |
| run_tp53_hotspot | 0 | 1.23 | [run_tp53_hotspot](reproduction_logs/run_tp53_hotspot.log) |
| run_survey | 0 | 1.1 | [run_survey](reproduction_logs/run_survey.log) |
| run_survey_validation | 0 | 0.11 | [run_survey_validation](reproduction_logs/run_survey_validation.log) |
| run_selection | 0 | 0.57 | [run_selection](reproduction_logs/run_selection.log) |
| run_kv7_provenance | 0 | 2.73 | [run_kv7_provenance](reproduction_logs/run_kv7_provenance.log) |
| run_rank_null | 0 | 0.43 | [run_rank_null](reproduction_logs/run_rank_null.log) |
| run_kv7_rank | 0 | 15.82 | [run_kv7_rank](reproduction_logs/run_kv7_rank.log) |
| run_scn5a | 0 | 16.16 | [run_scn5a](reproduction_logs/run_scn5a.log) |
| run_figure2 | 0 | 8.04 | [run_figure2](reproduction_logs/run_figure2.log) |
| run_clinical | 0 | 0.48 | [run_clinical](reproduction_logs/run_clinical.log) |
| run_bounds | 0 | 2.34 | [run_bounds](reproduction_logs/run_bounds.log) |
| run_clinvar | 0 | 0.45 | [run_clinvar](reproduction_logs/run_clinvar.log) |
| run_titration | 0 | 3.28 | [run_titration](reproduction_logs/run_titration.log) |
| run_titration_main | 0 | 1.48 | [run_titration_main](reproduction_logs/run_titration_main.log) |
| run_uncertainty | 0 | 3.59 | [run_uncertainty](reproduction_logs/run_uncertainty.log) |
| run_schematic | 0 | 7.61 | [run_schematic](reproduction_logs/run_schematic.log) |
| pytest | 0 | 8.58 | [pytest](reproduction_logs/pytest.log) |

Exact commands are recorded in `reproduction_results.json`. Tests: `test_dominance.py test_profile_revision.py test_assay_revision.py test_bounds.py test_tp53_anchors.py test_survey_revision.py test_revision_nov2026.py test_lean_proofs.py test_kcnq2_archive_nov20.py test_normalization_titration_revision.py test_scn5a_round2.py test_archive_reading_round2.py test_titration_round2.py`.

## Artifact comparisons

| Artifact | Comparison | Result |
| --- | --- | --- |
| scn5a_peak.csv | CSV | match |
| scn5a_controls.csv | CSV | match |
| scn5a_specifications.csv | CSV | match |
| kv7_controls.csv | CSV | match |
| kv7_rank.csv | CSV | match |
| deconvolution.csv | CSV | match |
| tp53_portability.csv | CSV | match |
| tp53_controls.csv | CSV | match |
| tp53_anchors.csv | CSV | match |
| clinical_null.csv | CSV | match |
| clinvar_crosswalk.csv | CSV | match |
| survey_decisions.csv | CSV | match |
| bounds.csv | CSV | match |
| tp53_hotspot.csv | CSV | match |
| survey_validation.csv | CSV | match |
| kv7_row_map.csv | CSV | match |
| deconvolution_origin_differences.csv | CSV | match |
| deconvolution_control_means.csv | CSV | match |
| deconvolution_se_comparator.csv | CSV | match |
| scn5a_robustness.csv | CSV | match |
| scn5a_joint_uncertainty.csv | CSV | match |
| scn5a_reference_dose.csv | CSV | match |
| survey_second_reading_comparison.csv | CSV | match |
| survey_blinded_reading_comparison.csv | CSV | match |
| titration_profiles.csv | CSV | match |
| titration_candidate_tests.csv | CSV | match |
| titration_sensitivity_ranges.csv | CSV | match |
| titration_fit_sensitivity_ranges.csv | CSV | match |
| titration_replication.csv | CSV | match |
| scn5a_specifications.json | JSON | match |
| kv7_controls.json | JSON | match |
| kv7_rank.json | JSON | match |
| deconvolution.json | JSON | match |
| tp53_portability.json | JSON | match |
| tp53_anchors.json | JSON | match |
| clinical_null.json | JSON | match |
| clinvar_crosswalk.json | JSON | match |
| survey_summary.json | JSON | match |
| survey_zheng_case.json | JSON | match |
| bounds.json | JSON | match |
| tp53_hotspot.json | JSON | match |
| survey_validation.json | JSON | match |
| selection.json | JSON | match |
| kv7_panels.json | JSON | match |
| kv7_provenance.json | JSON | match |
| rank_null.json | JSON | match |
| uncertainty.json | JSON | match |
| titration.json | JSON | match |
| scn5a_floor_scenarios.json | JSON | match |
| Schematic.png | PNG pixels | match |
| Schematic.pdf | PDF bytes excluding dates | match |
| Figure1.png | PNG pixels | match |
| Figure1.pdf | PDF bytes excluding dates | match |
| Figure2.png | PNG pixels | match |
| Figure2.pdf | PDF bytes excluding dates | match |
| Figure3.png | PNG pixels | match |
| Figure3.pdf | PDF bytes excluding dates | match |
| Titration.png | PNG pixels | match |
| Titration.pdf | PDF bytes excluding dates | match |
| TitrationMain.png | PNG pixels | match |
| TitrationMain.pdf | PDF bytes excluding dates | match |
| SCN5A_robustness.png | PNG pixels | match |
| SCN5A_robustness.pdf | PDF bytes excluding dates | match |

## Scope and limitations

This verifies execution and artifact reproduction with the archived inputs and cached ClinVar records in the existing environment. It does not test a fresh dependency installation, validate live clinical records, or establish biological assumptions. Figure3 is the analysis filename for supplementary S1 Fig. TIFF exports and manuscript rendering were not regenerated in this check. Skipped or unrequested drivers and their retained outputs are not counted as reproduced.

Pytest result:

```text
........................................................................ [ 48%]
...............................................s........................ [ 97%]
....                                                                     [100%]
147 passed, 1 skipped in 8.35s
```
