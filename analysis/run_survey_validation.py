#!/usr/bin/env python3
"""Documented AI-assisted metadata readings of the archive keyword screen.

    python3 run_survey.py              # writes survey_decisions.csv first
    python3 run_survey_validation.py

Every experiment set flagged by the strict or broad paired-design screen, by the ratio
screen, or by MaveDB's controlled "Phenotypic Assay Mechanism" keyword for a
dominant-negative effect was classified by reading its archived metadata.
The classifications, with a verbatim evidence quote for each, are stored in
data/archive/survey_validation_labels.csv.  This script

  * checks that every flagged experiment set has exactly one label;
  * checks that every evidence quote occurs verbatim in that experiment's
    screened metadata (score-set and experiment titles, descriptions,
    abstract and methods);
  * applies the author's individual check of the ten sets that any reading
    classified as verified or ambiguous (data/archive/survey_author_check.csv),
    whose calls are final;
  * reports the precision of each screen and which route detected each
    verified mutant/wild-type design.

The verified designs are those found by these routes plus a review of the
broader screen; they are not an independent gold standard, so detection
fractions are not unbiased recall estimates.  Writes survey_validation.json
and survey_validation.csv.
"""
import csv
from collections import Counter
import hashlib
import json
import math
from datetime import date
from pathlib import Path

from run_survey import DEFAULT_SOURCE, text_of

LABELS = Path(__file__).resolve().parent / "data" / "archive" / "survey_validation_labels.csv"
FIRST_PROVENANCE = LABELS.with_name("survey_first_reading_provenance.json")
SECOND_LABELS = LABELS.with_name("survey_second_reading.csv")
SECOND_PROVENANCE = LABELS.with_name("survey_second_reading_provenance.json")
BLIND_LABELS = LABELS.with_name("survey_blinded_reading.csv")
BLIND_PROVENANCE = LABELS.with_name("survey_blinded_reading_provenance.json")
AUTHOR_LABELS = LABELS.with_name("survey_author_check.csv")
AUTHOR_PROVENANCE = LABELS.with_name("survey_author_check_provenance.json")
DN_KEYWORD_KEY = "Phenotypic Assay Mechanism"


def wilson(k, n, z=1.959963984540054):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0.0, centre - half), min(1.0, centre + half)]


def dn_keyword(rec):
    for k in (rec.get("experiment") or {}).get("keywords") or []:
        kw = k.get("keyword") or {}
        if kw.get("key") == DN_KEYWORD_KEY and "ominant-negative" in (kw.get("label") or ""):
            return kw.get("label")
    return None


def nominal_agreement(first, repeat):
    """Unweighted agreement of two nominal labelings, without a gold standard."""
    if len(first) != len(repeat) or not first:
        raise ValueError("Agreement requires two nonempty labelings of equal length")
    n = len(first)
    a, b = Counter(first), Counter(repeat)
    labels = sorted(a.keys() | b.keys())
    observed = sum(x == y for x, y in zip(first, repeat)) / n
    expected = sum(a[label] * b[label] for label in labels) / n**2
    # A sole shared label gives no marginal variation with which to define kappa.
    kappa = (observed - expected) / (1 - expected) if expected < 1 else None
    return {
        "labels": labels,
        "first_marginal_counts": dict(sorted(a.items())),
        "repeat_marginal_counts": dict(sorted(b.items())),
        "confusion_matrix_rows_first_columns_repeat": [
            [sum(x == row and y == col for x, y in zip(first, repeat))
             for col in labels] for row in labels],
        "observed_agreement": observed,
        "chance_agreement_from_marginals": expected,
        "cohens_kappa_unweighted": kappa,
        "interpretation": "Repeatability of nominal design labels; not accuracy against an expert reference standard.",
    }


def second_reading(by_urn, exps, first, flagged, labels=SECOND_LABELS, provenance=SECOND_PROVENANCE):
    """Compare the recorded separate pass; do not call an exposed AI pass blinded."""
    rows = list(csv.DictReader(open(labels)))
    second = {row["experiment_prefix"]: row for row in rows}
    if len(second) != len(rows) or set(second) != flagged:
        raise ValueError("Second reading must label each flagged experiment set exactly once")
    for prefix, row in second.items():
        if row["source_urn"] not in exps[prefix]["sets"]:
            raise ValueError(f"Second-reading evidence belongs to another experiment: {prefix}")
        if row["evidence_quote"] not in text_of(by_urn[row["source_urn"]]):
            raise ValueError(f"Second-reading quote not found verbatim: {prefix}")
        if row["mutant_wt_design"] not in {"yes", "no", "ambiguous"}:
            raise ValueError(f"Invalid second-reading design label: {prefix}")
    provenance_data = json.loads(Path(provenance).read_text())
    if hashlib.sha256(Path(labels).read_bytes()).hexdigest() != provenance_data["distributed_labels_sha256"]:
        raise ValueError("Second-reading file differs from its recorded provenance")
    fields = ["mutant_wt_design", "alone_arm_in_experiment_set", "ratio_or_fraction_stated"]
    comparison = []
    for prefix in sorted(flagged):
        for field in fields:
            a, b = first[prefix][field], second[prefix][field]
            comparison.append({"experiment_prefix": prefix, "field": field,
                               "first_label": a, "second_label": b, "exact_agreement": a == b})
    return comparison, {
        "reader": provenance_data["reader"], "date": provenance_data["date"],
        "blinded": provenance_data["blinded"], "human_independent_validation": False,
        "prior_context": provenance_data["prior_context"],
        "n_read": len(second),
        "exact_agreement": {field: [sum(r["exact_agreement"] for r in comparison if r["field"] == field), len(second)] for field in fields},
        "disagreements": [r for r in comparison if not r["exact_agreement"]],
        "adjudication": [
            "TP53 00000068: 'yes (p53-null axis)' and 'yes' denote the same archived alone arm.",
            "CARD11 00001215/00001226: no separate archived alone arm is documented; edited-cell allele configurations remain unresolved. 'no' versus 'unknown' reflects these different questions.",
            "CARD11 ratio labels 'no' versus 'NA' agree that no numerical WT:mutant ratio is established for the ambiguous archived designs.",
        ],
    }


def blinded_reading(by_urn, exps, first, flagged):
    """Check the frozen blind repeat reading and retain unresolved disagreements."""
    provenance = json.loads(BLIND_PROVENANCE.read_text())
    if hashlib.sha256(BLIND_LABELS.read_bytes()).hexdigest() != provenance["output_sha256"]:
        raise ValueError("Blinded labels differ from the frozen reading")
    rows = list(csv.DictReader(BLIND_LABELS.open()))
    blind = {r["experiment_prefix"]: r for r in rows}
    if len(blind) != len(rows) or set(blind) != flagged:
        raise ValueError("Blinded reading must label each flagged experiment exactly once")
    comparison, final = [], {p: dict(row) for p, row in first.items()}
    for prefix, row in blind.items():
        matching = [by_urn[u] for u in exps[prefix]["sets"]
                    if u == row["source_urn"] or (by_urn[u].get("experiment") or {}).get("urn") == row["source_urn"]]
        if not matching or not any(row["evidence_quote"] in text_of(r) for r in matching):
            raise ValueError(f"Blinded evidence does not resolve verbatim within group: {prefix}")
        if row["mutant_wt_design"] not in {"yes", "no", "ambiguous"}:
            raise ValueError(f"Invalid blinded design label: {prefix}")
        a, b = first[prefix]["mutant_wt_design"], row["mutant_wt_design"]
        comparison.append({"experiment_prefix": prefix, "first_design": a,
                           "blinded_design": b, "exact_agreement": a == b,
                           "final_design": a if a == b else "ambiguous"})
        if a != b:
            final[prefix].update(mutant_wt_design="ambiguous", alone_arm_in_experiment_set="unknown",
                                 ratio_or_fraction_stated=row["ratio_or_fraction_stated"],
                                 evidence_quote=row["evidence_quote"],
                                 design_note="Unresolved first/blinded reading discrepancy; retained as ambiguous. " + row["rationale"])
    verified = [p for p in first if first[p]["mutant_wt_design"] == blind[p]["mutant_wt_design"] == "yes"]
    # These are the counts before the author's check, which author_check() applies next.
    arm_agree = sum(first[p]["alone_arm_in_experiment_set"].split(" ")[0] == blind[p]["alone_arm_in_experiment_set"] for p in verified)
    ratio_agree = sum(first[p]["ratio_or_fraction_stated"] == blind[p]["ratio_or_fraction_stated"] for p in verified)
    agreement = nominal_agreement(
        [r["first_design"] for r in comparison],
        [r["blinded_design"] for r in comparison],
    )
    return final, sorted(comparison, key=lambda r: r["experiment_prefix"]), {
        "reading_type": "blind repeat reading",
        "reader_identity": provenance["reader_identity"],
        "blinded_to_prior_labels_and_counts": not provenance["prior_archive_labels_or_counts_accessed"],
        "human_independent_validation": False,
        "frozen_csv_sha256": provenance["output_sha256"],
        "n_read": len(blind),
        "design_exact_agreement": [sum(r["exact_agreement"] for r in comparison), len(blind)],
        "design_agreement_statistics": agreement,
        "design_disagreements": [r for r in comparison if not r["exact_agreement"]],
        "verified_designs_alone_arm_semantic_agreement": [arm_agree, len(verified)],
        "verified_designs_ratio_label_agreement": [ratio_agree, len(verified)],
        "arm_comparison_scope": "Five mutually verified mixed designs only; first-pass NA labels outside this group do not answer the blinded codebook's broader standalone-arm question.",
        "final_adjudication": "Every design disagreement remains ambiguous; no forced consensus or post hoc replacement of frozen labels.",
        "final_counts": {label: sum(r["mutant_wt_design"] == label for r in final.values()) for label in ("yes", "no", "ambiguous")},
        "limits": "Blind repeat AI-assisted metadata reading; agreement measures repeatability, not accuracy against an expert standard or recall/prevalence. HAP1 labels in 00000097/00000658 infer haploid context without direct ploidy verification.",
    }


def author_check(final, first, flagged):
    """Apply the author's individual check of every set any reading did not call other.

    The author checked, on 22 November 2026, the ten flagged sets that the first,
    blind or second reading classified as verified or ambiguous; those calls are
    final. Every reading classified the remaining flagged sets as other.
    """
    provenance = json.loads(AUTHOR_PROVENANCE.read_text())
    if hashlib.sha256(AUTHOR_LABELS.read_bytes()).hexdigest() != provenance["labels_sha256"]:
        raise ValueError("Author check labels differ from their recorded provenance")
    rows = {r["experiment_prefix"]: r for r in csv.DictReader(AUTHOR_LABELS.open())}
    readings = [first] + [{r["experiment_prefix"]: r for r in csv.DictReader(path.open())}
                          for path in (BLIND_LABELS, SECOND_LABELS)]
    contested = {p for reading in readings for p, r in reading.items() if r["mutant_wt_design"] != "no"}
    if set(rows) != contested or not set(rows) <= flagged:
        raise ValueError("The author check must cover exactly the sets any reading did not call other")
    blind = readings[1]
    for prefix, row in rows.items():
        if row["mutant_wt_design"] not in {"yes", "no", "ambiguous"}:
            raise ValueError(f"Invalid author design label: {prefix}")
        if row["mutant_wt_design"] == first[prefix]["mutant_wt_design"]:
            final[prefix] = dict(first[prefix])
        final[prefix].update(mutant_wt_design=row["mutant_wt_design"],
                             alone_arm_in_experiment_set=row["alone_arm_in_experiment_set"],
                             ratio_or_fraction_stated=row["ratio_or_fraction_stated"],
                             design_note=f"Author check 2026-11-22 ({row['status']}). " + first[prefix]["design_note"])
    return final, {
        "reading_type": "author's individual check of the contested sets (not blinded)",
        "provenance": provenance,
        "n_checked": len(rows),
        "agrees_with_first_design": sum(rows[p]["mutant_wt_design"] == first[p]["mutant_wt_design"] for p in rows),
        "agrees_with_blind_design": sum(rows[p]["mutant_wt_design"] == blind[p]["mutant_wt_design"] for p in rows),
        "changed_from_blind_adjudication": sorted(p for p in rows if rows[p]["mutant_wt_design"] != blind[p]["mutant_wt_design"]),
        "final_counts": {label: sum(r["mutant_wt_design"] == label for r in final.values()) for label in ("yes", "no", "ambiguous")},
    }


def validate(source=DEFAULT_SOURCE, decisions="survey_decisions.csv", labels=LABELS):
    recs = json.loads(Path(source).read_text())
    recs = list(recs.values()) if isinstance(recs, dict) else recs
    by_urn = {r["urn"]: r for r in recs}
    rows = list(csv.DictReader(open(decisions)))
    label_rows = list(csv.DictReader(open(labels)))
    lab = {r["experiment_prefix"]: r for r in label_rows}
    if len(lab) != len(label_rows):
        raise ValueError("duplicate experiment-set validation labels")
    first_provenance = json.loads(FIRST_PROVENANCE.read_text())
    if Path(labels).resolve() == LABELS.resolve():
        if hashlib.sha256(Path(labels).read_bytes()).hexdigest() != first_provenance["labels_sha256"]:
            raise ValueError("First-reading labels differ from their recorded provenance")

    exps = {}
    for r in rows:
        e = exps.setdefault(r["experiment_prefix"], {"sets": [], "strict": 0, "loose": 0, "ratio": 0, "keyword": None})
        e["sets"].append(r["urn"])
        e["strict"] += r["paired_candidate_strict"] == "True"
        e["ratio"] += r["ratio_keyword_match"] == "True"
        e["loose"] += r["paired_candidate_loose"] == "True"
        kw = dn_keyword(by_urn[r["urn"]])
        if kw:
            e["keyword"] = kw

    flagged = {p for p, e in exps.items() if e["strict"] or e["loose"] or e["ratio"] or e["keyword"]}
    missing = sorted(flagged - set(lab))
    if missing:
        raise ValueError(f"flagged experiment sets without a validation label: {missing}")
    for p, l in lab.items():
        if p not in exps:
            raise ValueError(f"label for an experiment set absent from the screen: {p}")
        texts = [text_of(by_urn[u]) for u in exps[p]["sets"]]
        if not any(l["evidence_quote"] in t for t in texts):
            raise ValueError(f"evidence quote not found verbatim for {p}: {l['evidence_quote']!r}")
        if l["mutant_wt_design"] not in {"yes", "no", "ambiguous"}:
            raise ValueError(f"invalid design label for {p}")

    comparison, second = second_reading(by_urn, exps, lab, flagged)
    if hashlib.sha256(Path(source).read_bytes()).hexdigest() != json.loads(BLIND_PROVENANCE.read_text())["full_source_snapshot_sha256"]:
        raise ValueError("Blinded reading source snapshot differs from its frozen provenance")
    first_labels = {p: dict(r) for p, r in lab.items()}
    lab, blind_comparison, blind_result = blinded_reading(by_urn, exps, lab, flagged)
    lab, author_result = author_check(lab, first_labels, flagged)

    def count(pred):
        return sorted(p for p in lab if pred(p))

    strict_exp = count(lambda p: exps[p]["strict"] > 0)
    ratio_exp = count(lambda p: exps[p]["ratio"] > 0)
    loose_exp = count(lambda p: exps[p]["loose"] > 0)
    kw_exp = count(lambda p: exps[p]["keyword"] is not None)
    designs = count(lambda p: lab[p]["mutant_wt_design"] == "yes")
    ambiguous = count(lambda p: lab[p]["mutant_wt_design"] == "ambiguous")
    strict_tp = [p for p in strict_exp if p in designs]
    ratio_tp = [p for p in ratio_exp if p in designs and lab[p]["ratio_or_fraction_stated"] == "yes"]
    kw_tp = [p for p in kw_exp if p in designs]
    strict_sets = sum(exps[p]["strict"] for p in strict_exp)
    strict_sets_tp = sum(exps[p]["strict"] for p in strict_tp)
    stated = [p for p in designs if lab[p]["ratio_or_fraction_stated"] == "yes"]
    both_arms = [p for p in designs if lab[p]["alone_arm_in_experiment_set"].startswith("yes")]
    known_in_paper = ["urn:mavedb:00000068", "urn:mavedb:00000094", "urn:mavedb:00000664",
                      "urn:mavedb:00001235"]

    result = {
        "generated": date.today().isoformat(),
        "snapshot_retrieved": "2026-09-26",
        "labels_file": str(Path(labels).name),
        "final_classification_basis": "The author checked individually the ten sets that any reading classified as verified or ambiguous, and those calls are final; every reading classified the other 36 flagged sets as other.",
        "classifier_identity": "First set: 36 inherited labels with reader/model unrecorded plus 10 completed by OpenAI Codex (GPT-6); blind repeat and second reading: OpenAI Codex (GPT-6); final check of the ten contested sets: the author (not blinded).",
        "first_reading": first_provenance,
        "second_reading": second,
        "blinded_second_reading": blind_result,
        "author_check": author_result,
        "final_counts": author_result["final_counts"],
        "schema_scope": {
            "snapshot": "2026-09-26 short-record snapshot, 2819 records; 2809 with experiment abstract or methods screened",
            "short_record_definition": "Saved score-set metadata summaries and their nested experiment metadata; the snapshot does not include score-table columns or values.",
            "screened_fields": ["title", "shortDescription", "experiment.title", "experiment.shortDescription", "experiment.abstractText", "experiment.methodText"],
            "candidate_routes": ["strict paired-design text", "broader paired-design text", "ratio-related text", "controlled Phenotypic Assay Mechanism / Dominant-negative keyword"],
            "additional_candidate_route": "experiment.keywords: Phenotypic Assay Mechanism / Dominant-negative",
            "dedicated_achieved_fraction_field_in_snapshot": False,
            "limitation": "This is a statement about the archived record structure, not a live API schema or a validated protein-fraction measurement method.",
        },
        "n_labelled_experiment_sets": len(lab),
        "n_flagged_experiment_sets": len(flagged),
        "n_unclassified_flagged": len(flagged - set(lab)),
        "n_unflagged_experiment_sets": len(set(exps) - flagged),
        "n_unflagged_reviewed": 0,
        "recall_estimate": None,
        "prevalence_estimate": None,
        "broad_screen": {"experiment_sets": len(loose_exp),
                         "mutant_wt_designs": [p for p in loose_exp if p in designs]},
        "strict_screen": {
            "experiment_sets": len(strict_exp), "score_sets": strict_sets,
            "mutant_wt_designs": strict_tp, "n_true": len(strict_tp),
            "precision_experiment_sets": len(strict_tp) / len(strict_exp),
            "precision_ci_wilson": wilson(len(strict_tp), len(strict_exp)),
            "score_sets_in_true_designs": strict_sets_tp,
        },
        "ratio_screen": {
            "experiment_sets": len(ratio_exp),
            "n_stating_mutant_wt_ratio": len(ratio_tp),
            "precision_experiment_sets": len(ratio_tp) / len(ratio_exp),
            "precision_ci_wilson": wilson(len(ratio_tp), len(ratio_exp)),
            "match_contexts": {p: lab[p]["match_refers_to"] for p in ratio_exp},
        },
        "dominant_negative_keyword": {
            "experiment_sets": len(kw_exp), "mutant_wt_designs": kw_tp,
            "labels": {p: exps[p]["keyword"] for p in kw_exp},
        },
        "verified_mutant_wt_designs": {
            p: {"note": lab[p]["design_note"],
                "detected_by_strict_text": p in strict_exp,
                "detected_by_keyword": p in kw_exp,
                "detected_by_ratio_text": p in ratio_exp,
                "alone_arm_in_experiment_set": lab[p]["alone_arm_in_experiment_set"],
                "ratio_or_fraction_stated_in_screened_fields": lab[p]["ratio_or_fraction_stated"]}
            for p in designs},
        "n_verified_designs": len(designs),
        "manuscript_design_labels": {"yes": "verified", "no": "other", "ambiguous": "ambiguous"},
        "publication_table_compilation": {
            "file": "data/archive/survey_publication_designs.csv",
            "method": "A separate Codex-assisted source check followed archived publication identifiers to the relevant article methods, figures and supplements; verified designs and ambiguous records were then tabulated with evidence locators. This check did not alter the frozen metadata labels.",
            "source_studies_checked": ["Giacomelli 2018", "Zheng 2022", "Jiang 2022", "Boettcher 2019", "Meitlis 2020"],
            "limits": "KCNQ2 linkage uses deposit metadata and a candidate article abstract only; the three formerly disputed metadata-only records, classified as other by the author's check, have unverified publication links. Publications' separate validation assays are distinguished from archived assay arms.",
        },
        "verified_designs_with_both_archived_arms": both_arms,
        "n_verified_designs_with_both_archived_arms": len(both_arms),
        "detection_by_strict_text": [len(strict_tp), len(designs)],
        "detection_by_keyword": [len(kw_tp), len(designs)],
        "detection_by_either": [len(set(strict_tp) | set(kw_tp)), len(designs)],
        "designs_stating_ratio_or_fraction": [len(stated), len(designs)],
        "known_paired_studies_analysed_in_paper": {
            p: {"strict_text": p in strict_exp, "keyword": p in kw_exp} for p in known_in_paper},
        "ambiguous_designs": {p: lab[p]["design_note"] for p in ambiguous},
        "limits": [
            "Verified designs were found by the same screens and a review of the broader screen; detection fractions are not unbiased recall estimates.",
            "The two label sets use archived metadata. A separate publication/design table checks linked source methods where available; it is not part of the blinded reading. An other label means that a qualifying design is not documented in the screened fields.",
            "All strict, broad, ratio and controlled-keyword matches are classified; no probability sample of unflagged records was reviewed, so recall and archive-wide prevalence are unestimated.",
            "No interval is assigned to the five purposively selected designs: they are not a probability sample of the archive.",
            "A bicistronic construct fixes the DNA-level co-expression design but does not state an achieved protein fraction.",
        ],
    }
    Path("survey_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    with open("survey_blinded_reading_comparison.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(blind_comparison[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(blind_comparison)
    with open("survey_second_reading_comparison.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(comparison[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparison)
    with open("survey_validation.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["experiment_prefix", "n_score_sets", "strict_score_sets", "ratio_score_sets",
                    "dominant_negative_keyword", "mutant_wt_design", "match_refers_to",
                    "ratio_or_fraction_stated", "evidence_quote"])
        for p in sorted(lab):
            e, l = exps[p], lab[p]
            w.writerow([p, len(e["sets"]), e["strict"], e["ratio"], e["keyword"] or "",
                        l["mutant_wt_design"], l["match_refers_to"],
                        l["ratio_or_fraction_stated"], l["evidence_quote"]])
    return result


if __name__ == "__main__":
    r = validate()
    s, q = r["strict_screen"], r["ratio_screen"]
    print(f"strict screen: {s['n_true']} of {s['experiment_sets']} experiment sets are mutant/wild-type "
          f"designs (descriptive yield {s['precision_experiment_sets']:.2f}); "
          f"{s['score_sets_in_true_designs']} of {s['score_sets']} score sets")
    print(f"ratio screen: {q['n_stating_mutant_wt_ratio']} of {q['experiment_sets']} experiment sets state a "
          "mutant:wild-type ratio")
    print(f"verified designs: {r['n_verified_designs']}; strict text detects "
          f"{r['detection_by_strict_text'][0]}, keyword detects {r['detection_by_keyword'][0]}, "
          f"either {r['detection_by_either'][0]}; ratio or fraction stated in "
          f"{r['designs_stating_ratio_or_fraction'][0]}")
    print(f"ambiguous: {list(r['ambiguous_designs'])}")
    agreement = r["blinded_second_reading"]["design_agreement_statistics"]
    print(f"blind repeat reading: agreement {agreement['observed_agreement']:.6f}; "
          f"chance agreement {agreement['chance_agreement_from_marginals']:.6f}; "
          f"Cohen kappa {agreement['cohens_kappa_unweighted']:.6f}")
    print("wrote survey_validation.json, survey_validation.csv")
