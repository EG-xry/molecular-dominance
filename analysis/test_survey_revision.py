"""Regression checks for the archive screen and its scoped publication example."""
import json
import pytest

import run_survey as survey


@pytest.fixture(scope="module")
def screened_snapshot(tmp_path_factory):
    output = tmp_path_factory.mktemp("survey")
    return survey.run_survey(output_dir=output)


def test_original_keyword_counts_are_preserved(screened_snapshot):
    rows, summary, _ = screened_snapshot
    expected = {
        "n_score_sets": 2809,
        "n_experiments": 1288,
        "paired_candidates_scoresets_strict": 164,
        "paired_candidates_scoresets_loose": 189,
        "paired_candidate_experiments": 25,
        "physical_unit": 1699,
        "standardised_only": 113,
        "ratio_keyword_match": 33,
        "ratio_keyword_match_experiments": 13,
    }
    assert {key: summary[key] for key in expected} == expected
    assert len(rows) == len({r["urn"] for r in rows}) == 2809


def test_zheng_is_one_study_with_fourteen_metadata_records(screened_snapshot):
    _, _, case = screened_snapshot
    archive = case["archive_evidence"]
    assert archive["n_score_sets"] == 14
    assert archive["n_study_prefixes"] == 1
    assert archive["experiment_urns"] == ["urn:mavedb:00000094-a"]
    assert archive["paired_candidate_strict_count"] == 14
    assert archive["ratio_keyword_match_count"] == 0
    assert archive["literal_1_to_1_count_in_screened_fields"] == 0
    assert archive["equimolar_count_in_screened_fields"] == 0
    publication = case["publication_evidence"]
    assert publication["intended_dna_mass_ratio_WT_to_variant"] == "1:1"
    assert publication["WT_plasmid_micrograms_per_well"] == 0.8
    assert publication["variant_plasmid_micrograms_per_well"] == 0.8
    assert publication["achieved_assembly_competent_fraction"] is None
    assert case["supporting_document_verified"] is True


def test_explicit_input_ratio_is_detected_without_becoming_achieved_fraction(tmp_path):
    record = {
        "urn": "urn:mavedb:00000094-a-1",
        "title": "Example",
        "experiment": {"urn": "urn:mavedb:00000094-a",
                       "methodText": "Co-transfected variant with wild-type plasmid at 1:1."},
    }
    source = tmp_path / "snapshot.json"
    source.write_text(json.dumps([record]))
    _, summary, case = survey.run_survey(source, tmp_path)
    assert summary["ratio_keyword_match"] == 1
    assert case["archive_evidence"]["ratio_keyword_match_count"] == 1
    assert case["publication_evidence"]["achieved_assembly_competent_fraction"] is None
    assert case["conclusion"].startswith("Ratio-related language occurs")


def test_modified_supporting_source_fails_verification(tmp_path):
    original = json.loads((survey.DATA_DIR / "zheng2022_source.json").read_text())
    (tmp_path / original["source_file"]).write_bytes(b"not the source document")
    manifest = tmp_path / "source.json"
    manifest.write_text(json.dumps(original))
    with pytest.raises(ValueError, match="checksum"):
        survey.audit_zheng_case([], [], manifest)


def test_out_of_scope_archive_record_fields_do_not_enter_screen():
    record = {"title": "Example", "experiment": {"methodText": "co-transfected"},
              "unsearched_notes": "1:1"}
    assert survey.PAIRED_STRICT.search(survey.text_of(record))
    assert not survey.RATIO.search(survey.text_of(record))
