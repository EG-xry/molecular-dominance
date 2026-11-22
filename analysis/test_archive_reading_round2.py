"""Check recorded second-pass evidence and its limited sampling interpretation."""
import csv
import json
from pathlib import Path

import pytest
import run_survey_validation as validation

ROOT = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def validated(tmp_path_factory):
    output = tmp_path_factory.mktemp("second-reading")
    old = Path.cwd()
    import os
    try:
        os.chdir(output)
        return validation.validate(decisions=ROOT / "survey_decisions.csv")
    finally:
        os.chdir(old)


def test_all_candidates_have_separately_recorded_evidence(validated):
    second = validated["second_reading"]
    assert second["n_read"] == validated["n_flagged_experiment_sets"] == 46
    assert second["exact_agreement"]["mutant_wt_design"] == [46, 46]
    assert second["exact_agreement"]["alone_arm_in_experiment_set"] == [43, 46]
    assert second["blinded"] is False
    assert second["human_independent_validation"] is False
    assert len(second["disagreements"]) == 5


def test_archived_pair_count_and_purposive_sample_scope(validated):
    assert validated["verified_designs_with_both_archived_arms"] == [
        "urn:mavedb:00000068", "urn:mavedb:00000094", "urn:mavedb:00000664"]
    assert validated["designs_stating_ratio_or_fraction"] == [0, 5]
    assert "designs_stating_ratio_or_fraction_ci_wilson" not in validated
    assert validated["recall_estimate"] is None
    assert validated["prevalence_estimate"] is None
    blind = validated["blinded_second_reading"]
    assert blind["blinded_to_prior_labels_and_counts"] is True
    assert blind["human_independent_validation"] is False
    assert blind["design_exact_agreement"] == [43, 46]
    assert blind["verified_designs_alone_arm_semantic_agreement"] == [5, 5]
    assert blind["verified_designs_ratio_label_agreement"] == [5, 5]
    assert blind["final_counts"] == {"yes": 5, "no": 36, "ambiguous": 5}
    assert {r["experiment_prefix"] for r in blind["design_disagreements"]} == {
        "urn:mavedb:00001204", "urn:mavedb:00001206", "urn:mavedb:00001240"}
    # The author's check (below) resolves the three disputed sets as other.
    assert list(validated["ambiguous_designs"]) == ["urn:mavedb:00001215", "urn:mavedb:00001226"]


def test_author_check_is_final_for_every_contested_set(validated):
    """The author checked each set any reading did not call other; those calls are final."""
    check = validated["author_check"]
    assert check["n_checked"] == 10
    assert check["agrees_with_first_design"] == 10
    assert check["agrees_with_blind_design"] == 7
    assert check["changed_from_blind_adjudication"] == [
        "urn:mavedb:00001204", "urn:mavedb:00001206", "urn:mavedb:00001240"]
    assert check["provenance"]["blinded"] is False
    assert validated["final_counts"] == {"yes": 5, "no": 39, "ambiguous": 2}
    assert validated["n_verified_designs"] == 5
    assert validated["designs_stating_ratio_or_fraction"] == [0, 5]


def test_changed_author_check_is_rejected(tmp_path, monkeypatch):
    labels = tmp_path / "survey_author_check.csv"
    text = validation.AUTHOR_LABELS.read_text().replace("urn:mavedb:00001240,no", "urn:mavedb:00001240,ambiguous")
    labels.write_text(text)
    monkeypatch.setattr(validation, "AUTHOR_LABELS", labels)
    with pytest.raises(ValueError, match="Author check labels differ"):
        validation.validate(decisions=ROOT / "survey_decisions.csv")


def test_blind_repeat_kappa_uses_raw_pre_adjudication_labels(validated):
    """The three disputed records must not be reconciled before agreement is scored."""
    agreement = validated["blinded_second_reading"]["design_agreement_statistics"]
    assert agreement["first_marginal_counts"] == {"ambiguous": 2, "no": 39, "yes": 5}
    assert agreement["repeat_marginal_counts"] == {"ambiguous": 5, "no": 36, "yes": 5}
    assert agreement["observed_agreement"] == pytest.approx(43 / 46)
    assert agreement["chance_agreement_from_marginals"] == pytest.approx(1439 / 2116)
    assert agreement["cohens_kappa_unweighted"] == pytest.approx(539 / 677)
    assert agreement["confusion_matrix_rows_first_columns_repeat"] == [
        [2, 0, 0], [3, 36, 0], [0, 0, 5]]
    first = validated["first_reading"]
    assert first["inherited_records"]["n"] == 36
    assert first["completion_records"]["n"] == 10
    assert "Not recorded" in first["inherited_records"]["reader_system_and_model"]
    assert validated["blinded_second_reading"]["reading_type"] == "blind repeat reading"


def test_nominal_agreement_handles_undefined_kappa_and_invalid_pairs():
    assert validation.nominal_agreement(["no"] * 3, ["no"] * 3)["cohens_kappa_unweighted"] is None
    with pytest.raises(ValueError, match="equal length"):
        validation.nominal_agreement(["yes"], ["yes", "no"])
    with pytest.raises(ValueError, match="nonempty"):
        validation.nominal_agreement([], [])


def test_publication_validation_does_not_promote_ambiguous_archive_designs():
    rows = list(csv.DictReader((ROOT / "data/archive/survey_publication_designs.csv").open()))
    assert len(rows) == 10
    assert sum(r["both_arms_archived"] == "yes" for r in rows) == 3
    card11 = [r for r in rows if r["gene"] == "CARD11"]
    assert len({r["doi"] for r in card11}) == 1
    assert all(r["archive_design"] == "ambiguous" for r in card11)
    assert all("50:50" in r["intended_input_ratio"] for r in card11)
    kcnq2 = next(r for r in rows if r["gene"] == "KCNQ2")
    assert "no primary publication identifier" in kcnq2["publication_link_status"]


def test_changed_second_reading_evidence_is_rejected(tmp_path):
    records = json.loads((ROOT / "data/archive/scoresets.json").read_text())
    rows = list(csv.DictReader(validation.SECOND_LABELS.open()))
    rows[0]["evidence_quote"] = "THIS IS NOT SOURCE METADATA"
    labels = tmp_path / "changed.csv"
    with labels.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    experiments = {row["experiment_prefix"]: {"sets": [row["source_urn"]]} for row in rows}
    with pytest.raises(ValueError, match="quote not found"):
        validation.second_reading(records, experiments, {}, set(experiments), labels=labels)
