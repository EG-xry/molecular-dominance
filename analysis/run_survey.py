#!/usr/bin/env python3
"""Screen the MaveDB archive for paired designs, physical units and a stated ratio.

    python3 run_survey.py [path/to/scoresets.json]

This is a **keyword screen over metadata**, not an inventory of designs. It
reports which records mention the things an interpretable paired analysis
needs; whether a record that mentions them actually has them can only be
settled by validating the underlying publication and experiment. No systematic
validation of the candidate pool is claimed here.
Every regular expression and every per-record decision is written out to
survey_decisions.csv so the screen can be audited and re-tuned.
The additional survey_zheng_case.json records a purposively verified example
of published intended DNA-ratio information absent from the fields searched.

Source: MaveDB score-set metadata, retrieved 2026-09-26 (n = 2,809 with
retrievable metadata).
"""
import csv
import json
import re
import sys
from datetime import date
from pathlib import Path
import hashlib
import zipfile
import xml.etree.ElementTree as ET

DATA_DIR = Path(__file__).resolve().parent / "data" / "archive"
DEFAULT_SOURCE = DATA_DIR / "scoresets.json"

# --- the screen -----------------------------------------------------------
# PAIRED: language that would appear if a variant were assayed both alone and
# with wild type present.  Deliberately broad: a screen that misses candidates
# is worse than one that returns some to be discarded by hand.
# Two tiers.  STRICT asks for language that a paired alone/mixed design would
# have to use; LOOSE adds vocabulary that merely raises the possibility.  Both
# over-return: "co-transfect" describes a transfection method, not a design.
# These screens have unestimated sensitivity and specificity; neither
# bounds design prevalence.
PAIRED_STRICT = re.compile(
    r"co-?express(ed|ion)? with|in the presence of wild[- ]?type|"
    r"wild[- ]?type background|WT background|alone and (in |with )?(co|the presence)|"
    r"both alone and|co-?transfect|heterozygous (state|condition|background)|"
    r"1:1 (ratio|mix)|equimolar", re.I)
PAIRED_LOOSE = re.compile(
    r"co-?express|co-?transfect|cotransfect|heterozygo|homozygo|"
    r"wild[- ]?type background|WT background|dominant[- ]negative|"
    r"alone and|1:1 ratio|equimolar|haplo", re.I)
PAIRED = PAIRED_STRICT

# PHYSICAL: physical-unit keywords. Their presence does not establish matched
# score origins, a linear activity readout, or biological comparability.
PHYSICAL = re.compile(
    r"\bmV\b|\bpA/pF\b|\bpA\b|\bnA\b|% of wild|%WT|% WT|percent of wild|"
    r"kcal/mol|\bdG\b|\bddG\b|\bKd\b|\bKm\b|\bkcat\b|\buM\b|\bnM\b|"
    r"\bms\b|fold[- ]change of activity", re.I)

# RATIO: ratio-related language; a match need not describe an achieved fraction.
RATIO = re.compile(
    r"ratio of (plasmid|DNA|construct)|1:1|equimolar|allele[- ]specific|"
    r"expression ratio|relative expression|transfect(ed|ion) ratio|"
    r"stoichiometr", re.I)

STANDARDISED = re.compile(r"\bz[- ]score|log2|enrichment|fitness score|"
                          r"normali[sz]ed score", re.I)


def text_of(rec):
    """Everything this screen can see: the score set's title and description,
    and the parent experiment's title, description, abstract and methods.
    Score-column names are not present in the archive's short-record dump, so
    the screen does not see them."""
    exp = rec.get("experiment") or {}
    parts = [rec.get("title"), rec.get("shortDescription"),
             exp.get("title"), exp.get("shortDescription"),
             exp.get("abstractText"), exp.get("methodText")]
    return " ".join(p or "" for p in parts)


def audit_zheng_case(records, decisions, provenance_path=None):
    """One purposively checked publication-to-metadata case, not screen validation.

    The supporting document establishes an intended DNA dose only. No keyword
    outcome is converted into evidence that an achieved protein fraction was
    measured, absent from the publication, or missing across other studies.
    """
    provenance_path = Path(provenance_path or DATA_DIR / "zheng2022_source.json")
    evidence = json.loads(provenance_path.read_text())
    source_path = provenance_path.parent / evidence["source_file"]
    digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
    if digest != evidence["sha256"]:
        raise ValueError("Zheng supporting source checksum does not match provenance")
    with zipfile.ZipFile(source_path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    # Use only text nodes: Word can have other text-bearing XML elements.
    paragraphs = ["".join(t.text or "" for t in p.findall(".//w:t", ns))
                  for p in root.findall(".//w:p", ns)]
    dose = next((p for p in paragraphs if "0.8" in p and
                 "KCNQ4-Myc-mCherry" in p and "variant plasmids" in p), None)
    if dose is None or len(re.findall(r"0\.8\s*[μµ]g", dose)) != 2:
        raise ValueError("Expected equal WT and variant plasmid doses not found in source")
    selected = [r for r in records if r.get("urn", "").startswith("urn:mavedb:00000094-")]
    urns = {r["urn"] for r in selected}
    matched = [r for r in decisions if r["urn"] in urns]
    texts = [text_of(r) for r in selected]
    ratio_count = sum(bool(r["ratio_keyword_match"]) for r in matched)

    def role(rec):
        """Role of each deposited score set, read from its archive title and target."""
        title = (rec.get("title") or "")
        genes = {t.get("name") for t in rec.get("targetGenes") or []}
        if title.startswith("try_dataset"):
            return "trial upload"
        if "KCNQ1" in genes:
            return "KCNQ1 variants alone"
        if "co-expressing" in title:
            return "KCNQ4 variant co-expressed with wild type"
        if "synonymous" in title:
            return "KCNQ4 synonymous variants alone"
        if "homozygous state" in title:
            return "KCNQ4 variants alone"
        return "other"

    roles = {r["urn"]: role(r) for r in selected}
    role_counts = {k: sum(v == k for v in roles.values()) for k in sorted(set(roles.values()))}
    coexpressed_current_sets = [r for r in selected if roles[r["urn"]] ==
                                "KCNQ4 variant co-expressed with wild type"
                                and (r.get("title") or "").startswith("Normalized currents")]
    linked = sorted({p.get("identifier") for r in selected
                     for p in ((r.get("experiment") or {}).get("primaryPublicationIdentifiers") or [])
                     if p.get("identifier")})
    if not selected:
        conclusion = "This study has no eligible records in the supplied metadata snapshot."
    elif ratio_count:
        conclusion = ("Ratio-related language occurs in this study's supplied metadata; "
                      "a keyword match does not establish a measured protein fraction.")
    else:
        conclusion = ("Published intended DNA-ratio information is not captured by the "
                      "ratio keyword screen in this study's cached metadata.")
    return {
        "case_id": "zheng2022_kcnq4_publication_to_archive",
        "publication": evidence,
        "supporting_document_verified": True,
        "publication_evidence": {
            "section": evidence["evidence_section"],
            "intended_dna_mass_ratio_WT_to_variant": "1:1",
            "WT_plasmid_micrograms_per_well": 0.8,
            "variant_plasmid_micrograms_per_well": 0.8,
            "achieved_assembly_competent_fraction": None,
            "achieved_fraction_status": "not established by the documented DNA dose",
            "interpretation": evidence["interpretation"],
        },
        "archive_evidence": {
            "snapshot_retrieved": "2026-09-26",
            "score_set_urns": sorted(urns),
            "n_score_sets": len(selected),
            "n_study_prefixes": len({u.split("-")[0] for u in urns}),
            "experiment_urns": sorted({r["experiment"]["urn"] for r in selected}),
            "paired_candidate_strict_count": sum(bool(r["paired_candidate_strict"]) for r in matched),
            "ratio_keyword_match_count": ratio_count,
            "literal_1_to_1_count_in_screened_fields": sum(t.count("1:1") for t in texts),
            "equimolar_count_in_screened_fields": sum(len(re.findall("equimolar", t, re.I)) for t in texts),
            "score_set_roles": roles,
            "score_set_role_counts": role_counts,
            "n_coexpression_score_sets": role_counts.get("KCNQ4 variant co-expressed with wild type", 0),
            "n_variants_coexpressed": int(sum(r.get("numVariants") or 0 for r in coexpressed_current_sets)),
            "coexpression_selection": ("co-expression panels contain variants selected on their alone-state "
                                       "current or V1/2 (current-reduced, current-enhanced, or large V1/2 shifts), "
                                       "as stated in the deposited score-set titles"),
            "experiment_primary_publication_identifiers": linked,
            "fields_checked": ["title", "shortDescription", "experiment.title",
                               "experiment.shortDescription", "experiment.abstractText",
                               "experiment.methodText"],
        },
        "scope": ("One study contributes multiple score sets; these are not independent missingness observations. "
                  "Only the co-expression score sets concern a wild-type:variant input; the record links the "
                  "publication, so the dose is available one step away but has no machine-actionable field."),
        "conclusion": conclusion,
        "limits": ["This is a purposively checked case, not an estimate of screen sensitivity or reporting prevalence.",
                   "It does not establish a measured assembly-competent fraction in Zheng et al.",
                   "It concerns the cached fields searched, not all archive files, linked publications, or current records.",
                   "An absent metadata match is not evidence that the experiment lacked a stated input ratio."],
    }


def run_survey(source=DEFAULT_SOURCE, output_dir="."):
    source = Path(source)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    recs = json.loads(source.read_text())
    if isinstance(recs, dict):
        recs = list(recs.values())
    n_total = len(recs)
    # "retrievable metadata" = the parent experiment carries an abstract or methods
    recs = [r for r in recs
            if (r.get("experiment") or {}).get("abstractText")
            or (r.get("experiment") or {}).get("methodText")]
    print(f"{n_total} score sets in the dump; {len(recs)} carry retrievable "
          f"experiment metadata and are screened\n")

    rows, n_paired, n_phys, n_ratio, n_std, n_loose = [], 0, 0, 0, 0, 0
    experiments, exp_paired = set(), set()
    for r in recs:
        t = text_of(r)
        urn = r.get("urn", "")
        paired = bool(PAIRED_STRICT.search(t))
        loose = bool(PAIRED_LOOSE.search(t))
        exp = urn.split("-")[0]
        experiments.add(exp)
        if paired:
            exp_paired.add(exp)
        n_loose += loose
        phys = bool(PHYSICAL.search(t))
        ratio = bool(RATIO.search(t))
        std = bool(STANDARDISED.search(t))
        n_paired += paired
        n_phys += phys
        n_ratio += ratio
        n_std += std
        rows.append({
            "urn": urn,
            "title": (r.get("title", "") or "")[:120],
            "experiment_prefix": exp,
            "paired_candidate_strict": paired,
            "paired_candidate_loose": loose,
            "paired_match": (PAIRED_STRICT.search(t).group(0) if paired else ""),
            "physical_unit": phys,
            "unit_match": (PHYSICAL.search(t).group(0) if phys else ""),
            "ratio_keyword_match": ratio,
            "ratio_match": (RATIO.search(t).group(0) if ratio else ""),
            "standardised_only": std and not phys,
        })

    n = len(rows)
    with (output_dir / "survey_decisions.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "generated": date.today().isoformat(),
        "source": str(source.resolve()),
        "retrieved": "2026-09-26",
        "n_score_sets": n,
        "n_experiments": len(experiments),
        "paired_candidates_scoresets_strict": n_paired,
        "paired_candidates_scoresets_loose": n_loose,
        "paired_candidate_experiments": len(exp_paired),
        "paired_candidate_experiments_share": round(len(exp_paired) / len(experiments), 4),
        "physical_unit": n_phys,
        "physical_unit_share": round(n_phys / n, 4),
        "standardised_only": n_std,
        "ratio_keyword_match": n_ratio,
        "ratio_keyword_match_share_scoresets": round(n_ratio / n, 4),
        "ratio_keyword_match_experiments": None,   # filled below
        "caveat": ("keyword screen over metadata; a mention is a candidate, not a "
                   "confirmed design. Sensitivity and specificity are unestimated."),
        "patterns": {"paired_strict": PAIRED_STRICT.pattern,
                     "paired_loose": PAIRED_LOOSE.pattern, "physical": PHYSICAL.pattern,
                     "ratio": RATIO.pattern, "standardised": STANDARDISED.pattern},
    }
    # Every rate must be quoted against its own denominator. Score sets and
    # experiments are different universes: one experiment may deposit several score
    # sets, so 33/2,809 score sets is not 33/1,288 experiments.
    exp_ratio = {r["experiment_prefix"] for r in rows if r["ratio_keyword_match"]}
    summary["ratio_keyword_match_experiments"] = len(exp_ratio)
    summary["ratio_keyword_match_share_experiments"] = round(len(exp_ratio) / len(experiments), 4)
    summary["caveat_prevalence"] = (
        "metadata matches are not estimates or bounds for design prevalence or "
        "the prevalence of measured assembly-competent subunit fractions")
    case = audit_zheng_case(recs, rows)
    summary["verified_case_file"] = "survey_zheng_case.json"
    (output_dir / "survey_zheng_case.json").write_text(json.dumps(case, indent=2) + "\n")
    (output_dir / "survey_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    ne = len(experiments)
    print(f"MaveDB metadata screen, {n} score sets in {ne} experiments, "
          f"retrieved 2026-09-26\n")
    print(f"  paired CANDIDATE experiments (strict screen)  {len(exp_paired):5d}  "
          f"{len(exp_paired)/ne:6.1%}   of {ne} experiments")
    print(f"  the same, counted as score sets               {n_paired:5d}  {n_paired/n:6.1%}")
    print(f"  loose screen, score sets                      {n_loose:5d}  {n_loose/n:6.1%}")
    print(f"  physical-unit keyword match                      {n_phys:5d}  {n_phys/n:6.1%}")
    print(f"  standardised / enrichment only               {n_std:5d}  {n_std/n:6.1%}")
    print(f"  ratio/stoichiometry keyword, score sets        {n_ratio:5d}  {n_ratio/n:6.1%}  "
          f"of {n} score sets")
    print(f"  the same, counted as experiments              "
          f"{summary['ratio_keyword_match_experiments']:5d}  "
          f"{summary['ratio_keyword_match_share_experiments']:6.1%}  of {len(experiments)} experiments")
    print("\nA mention is a candidate, not a design. wrote survey_decisions.csv, "
          "survey_summary.json, survey_zheng_case.json")
    return rows, summary, case


if __name__ == "__main__":
    run_survey(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SOURCE)
