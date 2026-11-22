#!/usr/bin/env python3
"""Allele crosswalk for the SCN5A candidates: source label -> MANE -> ClinVar.

    python3 run_clinvar.py            # uses the cached records if present
    python3 run_clinvar.py --refresh  # re-queries NCBI E-utilities

Why a crosswalk is needed.  O'Neill et al. number SCN5A residues on the
2016-residue isoform a (NP_932173.1, encoded by NM_198056.3).  The MANE Select
transcript is NM_000335.5, whose protein NP_000326.2 is the 2015-residue
isoform b; the two differ by the single-codon deletion of Gln1077.  Every
label above residue 1077 therefore shifts DOWN BY ONE between the source
paper and MANE.  A previous version of this analysis recorded zero offsets
for all nine candidates while its own saved ClinVar titles showed the shift,
because the offset column stored the position that had been QUERIED rather
than the position in the record returned.  That error is corrected here by
checking each label's wild-type residue against both protein sequences.

Outputs clinvar_crosswalk.csv and clinvar_crosswalk.json.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}
AA1 = {v: k for k, v in AA3.items()}

SOURCE_TRANSCRIPT = "NM_198056.3"      # isoform a, 2016 aa
SOURCE_PROTEIN = "NP_932173.1"
MANE_TRANSCRIPT = "NM_000335.5"        # MANE Select, isoform b, 2015 aa
MANE_PROTEIN = "NP_000326.2"
Q_DEL = 1077                           # Gln1077 present in a, absent from b

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = "clinvar_cache"
REFRESH = "--refresh" in sys.argv

# All 51 variants in the panel.  The candidate set alone cannot answer whether
# the breach test flags clinically unsettled variants preferentially -- that
# needs the base rate across the whole panel, which is what a reviewer of the
# previous revision correctly identified as missing.
import pandas as _pd
_panel = _pd.read_csv("scn5a_peak.csv")
_spec = _pd.read_csv("scn5a_specifications.csv")
CANDIDATES = list(_panel.variant)
INTERSECTION = set(_spec.loc[_spec.in_intersection, "variant"])
BREACH_NUM = set(_spec.loc[_spec.het < 54.6, "variant"])          # numerical
BREACH_SIG = set(_spec.loc[_spec.p_fixed < 0.05, "variant"])      # spec (i)


def load_mane_protein():
    """NP_000326.2 sequence, from the copy saved during review."""
    for p in (os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "data", "refseq", "NP_000326.2.fasta"),
              "journal-review-2026-09-27/assay_NP_000326.2.fasta",
              "assay_NP_000326.2.fasta"):
        if os.path.exists(p):
            return "".join(l.strip() for l in open(p) if not l.startswith(">"))
    raise SystemExit("NP_000326.2 FASTA not found; fetch it from RefSeq first")


def mane_position(label, mane_seq):
    """Map a source label to its MANE position, verifying the residue.

    Returns (position, status).  Status is 'same' below the Gln1077 deletion,
    'shift-1' above it, or 'unverified' if neither position carries the
    expected wild-type residue."""
    wt, pos, mut = label[0], int(label[1:-1]), label[-1]
    at_same = mane_seq[pos - 1] if pos <= len(mane_seq) else None
    at_shift = mane_seq[pos - 2] if 1 < pos <= len(mane_seq) + 1 else None
    if at_same == wt and pos <= Q_DEL:
        return pos, "same"
    if at_shift == wt and pos > Q_DEL:
        return pos - 1, "shift-1"
    if at_same == wt:
        return pos, "same"
    if at_shift == wt:
        return pos - 1, "shift-1"
    return None, "unverified"


def eutils(path, **params):
    params.setdefault("tool", "dominance-crosswalk")
    url = EUTILS + path + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()


def is_exact(title, hgvs_p):
    """True if a ClinVar title is the queried substitution on the MANE transcript.

    A text search also returns other substitutions at the same codon and large
    deletions spanning it, and the first hit is not always the queried allele."""
    return title.startswith(f"{MANE_TRANSCRIPT}({'SCN5A'})") and title.endswith(f"({hgvs_p})")


def clinvar_lookup(hgvs_p, gene="SCN5A"):
    """Search ClinVar by MANE protein HGVS and summarise the exact match, if any."""
    os.makedirs(CACHE, exist_ok=True)
    key = os.path.join(CACHE, re.sub(r"[^A-Za-z0-9]", "_", hgvs_p) + ".json")
    if os.path.exists(key) and not REFRESH:
        cached = json.load(open(key))
        # Records cached before the exact match check are reused only if they
        # are the queried substitution; otherwise the search is repeated.
        if not cached.get("variation_id") or is_exact(cached.get("title", ""), hgvs_p):
            return cached
    term = f'{gene}[gene] AND "{hgvs_p}"'
    rec = {"term": term, "retrieved": date.today().isoformat()}
    try:
        es = json.loads(eutils("esearch.fcgi", db="clinvar", term=term,
                               retmode="json", retmax=20))
        ids = es["esearchresult"].get("idlist", [])
        rec["ids"] = ids
        if ids:
            time.sleep(0.35)
            su = json.loads(eutils("esummary.fcgi", db="clinvar", id=",".join(ids),
                                   retmode="json"))["result"]
            rec["titles"] = {i: su[i].get("title", "") for i in ids}
            exact = [i for i in ids if is_exact(su[i].get("title", ""), hgvs_p)]
            rec["exact_matches"] = exact
        if ids and exact:
            r = su[exact[0]]
            germ = r.get("germline_classification", {})
            rec.update({
                "variation_id": exact[0],
                "title": r.get("title", ""),
                "classification": germ.get("description", "not provided"),
                "review_status": germ.get("review_status", ""),
                "last_evaluated": germ.get("last_evaluated", ""),
                "conditions": sorted({t.get("trait_name", "")
                                      for s in germ.get("trait_set", [])
                                      for t in [s]} - {""}),
                "canonical_spdi": r.get("variation_set", [{}])[0].get("canonical_spdi", ""),
                "nucleotide_hgvs": next(
                    (v.get("cdna_change", "") for v in r.get("variation_set", [])), ""),
            })
    except Exception as e:                                   # network optional
        rec["error"] = f"{type(e).__name__}: {e}"
    json.dump(rec, open(key, "w"), indent=2)
    time.sleep(0.35)
    return rec


def codon_neighbours(mpos, gene="SCN5A"):
    """Which other substitutions at this codon does ClinVar hold?

    Absence of a protein-HGVS match is not evidence that a codon is
    unrepresented, so a candidate reported as absent is checked at the codon
    level before the absence is stated."""
    os.makedirs(CACHE, exist_ok=True)
    key = os.path.join(CACHE, f"codon_{mpos}.json")
    if os.path.exists(key) and not REFRESH:
        return json.load(open(key))
    out = {"queried": f'{gene}[gene] AND "{AA1[mane[mpos-1]]}{mpos}"',
           "retrieved": date.today().isoformat(), "alleles": []}
    try:
        es = json.loads(eutils("esearch.fcgi", db="clinvar", term=out["queried"],
                               retmode="json", retmax=20))
        ids = es["esearchresult"].get("idlist", [])
        if ids:
            time.sleep(0.35)
            su = json.loads(eutils("esummary.fcgi", db="clinvar",
                                   id=",".join(ids), retmode="json"))["result"]
            out["alleles"] = [
                {"variation_id": i, "title": su[i].get("title", ""),
                 "classification": su[i].get("germline_classification", {})
                                        .get("description", "")}
                for i in ids]
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    json.dump(out, open(key, "w"), indent=2)
    time.sleep(0.35)
    return out


mane = load_mane_protein()
print(f"{MANE_PROTEIN}: {len(mane)} aa   {SOURCE_PROTEIN}: 2016 aa "
      f"(differ by Gln{Q_DEL} deletion)\n")

rows = []
for lab in CANDIDATES:
    wt, pos, mut = lab[0], int(lab[1:-1]), lab[-1]
    if wt not in AA1 or mut not in AA1:
        # W822X is a nonsense control, not a missense variant; it has no
        # p.XnnnY designation and is excluded from the crosswalk.
        continue
    mpos, status = mane_position(lab, mane)
    hgvs_p = f"p.{AA1[wt]}{mpos}{AA1[mut]}" if mpos else ""
    cv = clinvar_lookup(hgvs_p) if hgvs_p else {}
    rows.append({
        "source_label": lab,
        "source_transcript": SOURCE_TRANSCRIPT,
        "source_protein": SOURCE_PROTEIN,
        "source_position": pos,
        "mane_transcript": MANE_TRANSCRIPT,
        "mane_protein": MANE_PROTEIN,
        "mane_position": mpos,
        "offset": (pos - mpos) if mpos else None,
        "offset_status": status,
        "mane_hgvs_p": f"{MANE_PROTEIN}:{hgvs_p}" if hgvs_p else "",
        "clinvar_variation_id": cv.get("variation_id", ""),
        "clinvar_title": cv.get("title", "not found"),
        "clinvar_hgvs_c": cv.get("nucleotide_hgvs", ""),
        "canonical_spdi": cv.get("canonical_spdi", ""),
        "classification": cv.get("classification", "absent from ClinVar"),
        "review_status": cv.get("review_status", ""),
        "last_evaluated": cv.get("last_evaluated", ""),
        "conditions": "; ".join(cv.get("conditions", [])),
        "retrieved": cv.get("retrieved", date.today().isoformat()),
        "in_intersection": lab in INTERSECTION,
        "breach_numerical": lab in BREACH_NUM,
        "breach_significant": lab in BREACH_SIG,
        "query_error": cv.get("error", ""),
    })
    # A candidate with no protein-HGVS hit gets a codon-level check, so the
    # absence claim is about one substitution rather than the whole codon.
    if mpos and not cv.get("variation_id"):
        nb = codon_neighbours(mpos)
        rows[-1]["codon_alleles_in_clinvar"] = "; ".join(
            f"{a['title'].split('(')[-1].strip(')')} [{a['classification']}]"
            for a in nb.get("alleles", [])) or "none"
        rows[-1]["classification"] = (
            "this substitution absent from ClinVar; codon represented"
            if nb.get("alleles") else "codon absent from ClinVar")

import pandas as pd
df = pd.DataFrame(rows)
df.to_csv("clinvar_crosswalk.csv", index=False)
json.dump({"generated": date.today().isoformat(),
           "source": {"transcript": SOURCE_TRANSCRIPT, "protein": SOURCE_PROTEIN},
           "mane": {"transcript": MANE_TRANSCRIPT, "protein": MANE_PROTEIN},
           "rows": rows}, open("clinvar_crosswalk.json", "w"), indent=2)

UNSETTLED = ("Uncertain", "Conflicting", "not provided", "absent",
             "no classification", "this substitution absent", "codon absent")
df["unsettled"] = df.classification.str.startswith(UNSETTLED)
print(df[df.in_intersection | df.source_label.isin(
    ["C335R", "G386R", "W879R", "N1380K", "G1743E"])][
    ["source_label", "mane_position", "offset", "offset_status",
     "clinvar_variation_id", "classification"]].to_string(index=False))

print("\n=== BASE RATE: does breaching flag unsettled variants preferentially? ===")
from scipy import stats as _st
for name, s_ in (("numerical (het < 54.6)", "breach_numerical"),
                 ("significant (spec i)", "breach_significant"),
                 ("intersection of five", "in_intersection")):
    b = df[df[s_]]
    n = df[~df[s_]]
    if len(b) == 0:
        continue
    tab = [[int(b.unsettled.sum()), int((~b.unsettled).sum())],
           [int(n.unsettled.sum()), int((~n.unsettled).sum())]]
    orr, p = _st.fisher_exact(tab)
    print(f"  {name:24s} breaching {b.unsettled.sum():2d}/{len(b):2d} "
          f"({b.unsettled.mean():5.1%})  non-breaching {n.unsettled.sum():2d}/{len(n):2d} "
          f"({n.unsettled.mean():5.1%})  OR {orr:5.2f}  p {p:.3f}")
print(f"  panel base rate: {df.unsettled.sum()}/{len(df)} = {df.unsettled.mean():.1%} unsettled")
for r in rows:
    if r.get("codon_alleles_in_clinvar"):
        print(f"\ncodon check for {r['source_label']} "
              f"({MANE_PROTEIN} codon {r['mane_position']}): "
              f"{r['codon_alleles_in_clinvar']}")
n_shift = int((df.offset == 1).sum())
df.to_csv("clinvar_crosswalk.csv", index=False)
print(f"\n{n_shift}/{len(df)} labels shift by one residue between "
      f"{SOURCE_PROTEIN} and {MANE_PROTEIN}; "
      f"{int((df.offset == 0).sum())} are unchanged.")
print("wrote clinvar_crosswalk.csv and clinvar_crosswalk.json")
