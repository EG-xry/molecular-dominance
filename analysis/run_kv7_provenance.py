#!/usr/bin/env python3
"""run_kv7_provenance.py -- re-extraction of the two historical Kv7 extract tables
from the vendored supplementary workbooks (revision of 8 November 2026).

Outputs: kv7_provenance.json, kv7_provenance.log, kv7_row_map.csv, kv7_source_rule_log.csv.
The script exits non-zero unless every extract row is reproduced exactly.

Original description:

Rebuilds, from the original supplementary workbooks only,

  dominance/data/kv7/fit_KCNQ4_primary.csv        (Zheng et al. 2022, Genome Res 32:1573)
  dominance/data/kv7/fit_KCNQ2_KCNQ3_primary.csv  (Vanoye et al. 2022, JCI Insight 7:e156314)

using the explicit inclusion rules K1-K4 and V1-V5 written out below; compares every
numeric field (x, sx, y, sy, cell counts, residue positions, source-row numbers) with the
vendored extracts at a relative tolerance of 1e-9; checks exhaustiveness in both
directions; and re-applies run_kv7.py's own filters plus dominance.estimate() to
reproduce the reported panel summaries.

Inclusion rules (applied in order; each candidate is logged with the first rule it fails)

  KCNQ4, unit = one co-expression row of Supplemental_Dataset_C/D/E/F (KCNQ4_DMS.xlsx)
    K1  mixed arm from Data Set C (current-reduced) or E (current-enhanced); D and F are dropped
        because those variants were co-expressed because of a >10 mV homozygous V1/2 shift
    K2  variant is a missense HGVS name p.Xxx#Yyy with ref != alt
    K3  the variant has a homozygous V1/2 row in Supplemental_Dataset_B (exact-name join)
    K4  x = B!D, sx = B!E, y = C|E!E, sy = C|E!F are all numeric
    cells: mixed n = C|E!D; alone n = Supplemental_Dataset_A!D (current n; B reports no n)

  KCNQ2/KCNQ3, unit = one column-D label of sheet 'Relative to WT channel' (SDT3 alone, SDT4 mixed)
    V1  label present in both tables (paired by the column-D label, not by columns A-C)
    V2  label is a single missense substitution, e.g. R201H
    V3  x = SDT3!N, sx = SDT3!O, y = SDT4!O, sy = SDT4!P are all numeric (not 'nd')
    V4  phenotype (column F) is BFNE, DEE or BFNE/DEE in both tables (drops PV and 'unknown')
    V5  alone arm recorded in CHO-Q3 cells (drops A193D, R201C, P335L: CHO-K1 + KCNQ3)
    cells: alone n = SDT3!V, mixed n = SDT4!W (voltage-dependence block)

Repository is read-only: repository files are only opened for reading, no bytecode is
written (sys.dont_write_bytecode, and dominance.py is executed from its source text rather
than imported), and every output goes to the folder that contains this script, which is
refused if it lies inside the repository.

Run:   cd <folder containing this file> && python3 -I rederive_kv7.py
Needs: Python 3.10, openpyxl, pandas, numpy, scipy (dominance.py imports scipy).
"""
import sys

sys.dont_write_bytecode = True  # never write __pycache__ next to a repository file

import csv
import hashlib
import html
import json
import math
import re
import types
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

# --------------------------------------------------------------------------- paths
HERE = Path(__file__).resolve().parent
SOURCE = HERE / "data" / "kv7" / "source"
ZHENG = SOURCE / "zheng2022_supplemental_dataset.xlsx"
SDT3 = SOURCE / "vanoye2022_sdt3.xlsx"
SDT4 = SOURCE / "vanoye2022_sdt4.xlsx"
VANOYE_TEXT = SOURCE / "vanoye2022_article_text.txt"
EXT_K4 = HERE / "data/kv7/fit_KCNQ4_primary.csv"
EXT_VQ = HERE / "data/kv7/fit_KCNQ2_KCNQ3_primary.csv"
DOMINANCE = HERE / "dominance.py"
RUN_KV7 = HERE / "run_kv7.py"
INPUTS = [ZHENG, SDT3, SDT4, VANOYE_TEXT, EXT_K4, EXT_VQ, DOMINANCE, RUN_KV7]
# The historical extracts record the locations of the workbooks when they were
# made; the vendored copies are verified against the published files' checksums.
RECORDED = {
    "ZHENG": "/Users/eric/Downloads/PIEZO1-main/kcnq4-dms/KCNQ4_DMS.xlsx",
    "SDT3": "/Users/eric/Downloads/PIEZO1-main/deepen/codex/scratch/literature/vanoye_sdt3.xlsx",
    "SDT4": "/Users/eric/Downloads/PIEZO1-main/deepen/codex/scratch/literature/vanoye_sdt4.xlsx",
}
EXPECTED_SHA256 = {
    ZHENG: "bdfe8de2bedbeea7a4e64d152ef97ffb77dabcfded58a744517a2f1a105d7942",
    SDT3: "bfee4246fa5214db0e375ecf340adf051ee774922168021a7b97f3e1a300d41a",
    SDT4: "70569b45a9efb6d17698a9c453079b3e6808722e8c9c6e31dfa786d68dd989b6",
}
OUT = HERE

RTOL = 1e-9
REPORTED = {  # MANUSCRIPT.md line 158
    "KCNQ4": dict(b=0.36, n=190, residues=141, ci=(0.28, 0.44)),
    "KCNQ2_KCNQ3": dict(b=0.45, n=34, residues=27, ci=(0.36, 0.56)),
}

AA3 = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
       "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
HGVS3 = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
ONE_LETTER_SUB = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")

LOG = []


def say(*parts):
    line = " ".join(str(p) for p in parts)
    print(line, flush=True)
    LOG.append(line)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def nan_like(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def same_num(a, b):
    if nan_like(a) or nan_like(b):
        return nan_like(a) and nan_like(b)
    return math.isclose(float(a), float(b), rel_tol=RTOL, abs_tol=0.0)


def rel_diff(a, b):
    if nan_like(a) or nan_like(b):
        return 0.0 if (nan_like(a) and nan_like(b)) else math.inf
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


class Book:
    """A workbook opened twice: cached values, and stored formulas (to prove inputs are literals)."""

    def __init__(self, path):
        self.path = Path(path)
        self.v = openpyxl.load_workbook(self.path, data_only=True)
        self.f = openpyxl.load_workbook(self.path, data_only=False)

    def cell(self, sheet, row, col):
        value = self.v[sheet].cell(row, col).value
        stored = self.f[sheet].cell(row, col).value
        if isinstance(stored, str) and stored.startswith("="):
            raise ValueError(f"{self.path.name}!{sheet}!{openpyxl.utils.get_column_letter(col)}{row} "
                             f"is a formula ({stored}); only literal cells are accepted")
        return value


def rel(path):
    """Package-relative path for portable outputs."""
    return str(Path(path).resolve().relative_to(HERE))


def col(letter):
    return openpyxl.utils.column_index_from_string(letter)


# =========================================================================== KCNQ4
K4_SHEETS = {s: f"Supplemental_Dataset_{s}" for s in "ABCDEF"}
K4_HEADERS = {  # (sheet letter) -> {column letter: header text in row 3}
    "A": {"A": "Variant", "B": "Normalized current (Imut/Iwt, at +40 mV)", "D": "Number of cells"},
    "B": {"A": "Variant", "D": "Δ V1/2 (mV)", "E": "SEM (ΔV1/2)"},
    "C": {"A": "Variant", "D": "Number of cells", "E": "Δ V1/2 (mV)", "F": "SEM (ΔV1/2)"},
    "D": {"A": "Variant", "D": "Number of cells", "E": "Δ V1/2 (mV)", "F": "SEM (ΔV1/2)"},
    "E": {"A": "Variant", "D": "Number of cells", "E": "Δ V1/2 (mV)", "F": "SEM (ΔV1/2)"},
    "F": {"A": "Variant", "D": "Number of cells", "E": "Δ V1/2 (mV)", "F": "SEM (ΔV1/2)"},
}
K4_INCLUDED_MIXED = ("C", "E")  # rule K1


def k4_parse(book):
    """Return {letter: {"rows": {variant: row}, "title":..., "condition":..., "non_variant": [...]}}."""
    out = {}
    for letter, name in K4_SHEETS.items():
        ws = book.v[name]
        for c, text in K4_HEADERS[letter].items():
            got = book.cell(name, 3, col(c))
            assert got == text, (name, c, got, text)
        rows, non_variant = {}, []
        for r in range(4, ws.max_row + 1):
            a = book.cell(name, r, 1)
            if isinstance(a, str) and a.strip().startswith("p."):
                v = a.strip()
                assert v not in rows, f"duplicate variant {v} in {name}"
                rows[v] = r
            elif a is not None:
                non_variant.append((r, a))
        out[letter] = dict(rows=rows, title=book.cell(name, 1, 1), condition=book.cell(name, 2, 2),
                           non_variant=non_variant, sheet=name)
    return out


def k4_hgvs(v):
    m = HGVS3.match(v)
    if not m:
        return None
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if ref in AA3 and alt in AA3 and ref != alt:
        return ref, pos, alt
    return None


def k4_rederive(book, parsed):
    """Apply K1-K4 to every co-expression row (sheets C, D, E, F). Return (included df, rule log)."""
    A, B = parsed["A"], parsed["B"]
    log, keep = [], []
    for letter in "CDEF":
        sh = parsed[letter]
        for variant, row in sh["rows"].items():
            rec = dict(panel="KCNQ4", unit="co-expression row", variant=variant,
                       source_file=RECORDED["ZHENG"], sheet=sh["sheet"], source_row=row,
                       alone_sheet=B["sheet"], alone_source_row=B["rows"].get(variant, ""))
            if letter not in K4_INCLUDED_MIXED:
                rec["decision"] = "excluded:K1_gating_selected_sheet"
            elif k4_hgvs(variant) is None:
                rec["decision"] = "excluded:K2_not_missense"
            elif variant not in B["rows"]:
                rec["decision"] = "excluded:K3_no_alone_V1/2_in_Dataset_B"
            else:
                br = B["rows"][variant]
                x = book.cell(B["sheet"], br, col("D"))
                sx = book.cell(B["sheet"], br, col("E"))
                y = book.cell(sh["sheet"], row, col("E"))
                sy = book.cell(sh["sheet"], row, col("F"))
                if not all(is_num(t) for t in (x, sx, y, sy)):
                    rec["decision"] = "excluded:K4_nonnumeric_x_sx_y_sy"
                else:
                    rec["decision"] = "included"
                    ar = A["rows"].get(variant)
                    keep.append(dict(
                        variant=variant, y=y, sy=sy,
                        mixed_reported_cells=book.cell(sh["sheet"], row, col("D")),
                        sheet=letter, mixed_source_row=row, x=x, sx=sx, alone_source_row=br,
                        alone_current_cells=(book.cell(A["sheet"], ar, col("D")) if ar else np.nan),
                        alone_current_source_row=(ar if ar else np.nan),
                        pos=k4_hgvs(variant)[1], alone_gating_cells=np.nan,
                        source_file=RECORDED["ZHENG"], alone_sheet=B["sheet"], mixed_sheet=sh["sheet"],
                        x_column=K4_HEADERS["B"]["D"], sx_column=K4_HEADERS["B"]["E"],
                        y_column=K4_HEADERS[letter]["E"], sy_column=K4_HEADERS[letter]["F"],
                        primary=True))
            log.append(rec)
    return pd.DataFrame(keep), log


# =========================================================================== Vanoye
V_SHEET = "Relative to WT channel"
V_LAYOUT = {
    # header row, data-start row, residue, wt aa, var aa, label, phenotype, mean, sem, n(gating)
    "SDT3": dict(path=SDT3, header_row=12, first=13, A="A", B="B", C="C", D="D", F="F",
                 mean="N", sem="O", n="V",
                 checks={("N", 11): "ΔV½ from WT", ("N", 12): "mean", ("O", 12): "sem",
                         ("V", 12): "n", ("N", 9): "Voltage-dependence of Activation",
                         ("A", 12): "Residue #", ("D", 12): "Variant", ("F", 12): "Phenotype",
                         ("A", 4): "CHO-Q3 cells electroporated with KCNQ2-variant"}),
    "SDT4": dict(path=SDT4, header_row=12, first=13, A="A", B="B", C="C", D="D", F="F",
                 mean="O", sem="P", n="W",
                 checks={("O", 9): "ΔV½ from WT+WT", ("O", 10): "mean", ("P", 10): "sem",
                         ("W", 10): "n", ("O", 7): "Voltage-dependence of Activation",
                         ("A", 12): "Residue #", ("D", 12): "Variant", ("F", 12): "Phenotype",
                         ("A", 4): "CHO-Q3 cells co-electroporated with KCNQ2-variant + KCNQ2-WT"}),
}
EPILEPSY = {"BFNE", "DEE", "BFNE/DEE"}  # rule V4
# Rule V5. Vanoye et al. 2022, Methods ("Plasmids and transfection"): "Of the KCNQ2 variants, 3
# (A193D, R201C, and P335L) yielded too few positively transfected cells for automated patch clamp
# recording when expressed as homozygous channels in CHO-Q3 cells. We overcame this limitation by
# transiently cotransfecting CHO-K1 cells with each of these variants along with KCNQ3 ..."
CONSTRUCT_EXCEPTIONS = ("A193D", "R201C", "P335L")
V5_SENTENCE = re.compile(r"Of the KCNQ2 variants, 3 \(A193D, R201C, and P335L\) yielded too few "
                         r"positively transfected cells.{0,200}?transiently cotransfecting CHO-K1 cells")


def v_parse(book, key):
    L = V_LAYOUT[key]
    for (c, r), text in L["checks"].items():
        got = book.cell(V_SHEET, r, col(c))
        assert isinstance(got, str) and got.strip() == text, (key, c, r, got, text)
    ws = book.v[V_SHEET]
    rows, wt_rows, anomalies = {}, [], []
    for r in range(L["first"], ws.max_row + 1):
        label = book.cell(V_SHEET, r, col(L["D"]))
        if not (isinstance(label, str) and label.strip()):
            continue
        label = label.strip()
        if label.lower().startswith("wt"):
            wt_rows.append(dict(row=r, label=label,
                                mean=book.cell(V_SHEET, r, col(L["mean"])),
                                sem=book.cell(V_SHEET, r, col(L["sem"])),
                                n=book.cell(V_SHEET, r, col(L["n"]))))
            continue
        assert label not in rows, f"duplicate label {label} in {key}"
        rows[label] = r
        res = book.cell(V_SHEET, r, col(L["A"]))
        wt = book.cell(V_SHEET, r, col(L["B"]))
        var = book.cell(V_SHEET, r, col(L["C"]))
        m = ONE_LETTER_SUB.match(label)
        if m:  # record where the residue / aa columns disagree with the label itself
            if not (is_num(res) and float(res) == int(m.group(2))):
                anomalies.append(f"{key} row {r} {label}: residue column A = {res!r}")
            if wt != m.group(1):
                anomalies.append(f"{key} row {r} {label}: WT aa column B = {wt!r}")
            if var != m.group(3):
                anomalies.append(f"{key} row {r} {label}: Var aa column C = {var!r}")
        else:
            anomalies.append(f"{key} row {r} {label}: not a single substitution "
                             f"(A={res!r}, B={wt!r}, C={var!r})")
    title = book.cell(V_SHEET, 1, 1)
    return dict(rows=rows, wt_rows=wt_rows, anomalies=anomalies, title=title,
                condition=book.cell(V_SHEET, 4, 1))


def v_value(book, key, row, what):
    return book.cell(V_SHEET, row, col(V_LAYOUT[key][what]))


def v_rederive(b3, b4, p3, p4):
    log, keep = [], []
    labels = list(p3["rows"]) + [v for v in p4["rows"] if v not in p3["rows"]]
    for label in labels:
        r3, r4 = p3["rows"].get(label), p4["rows"].get(label)
        rec = dict(panel="KCNQ2_KCNQ3", unit="variant label", variant=label,
                   source_file=f"{RECORDED['SDT3']} | {RECORDED['SDT4']}", sheet=V_SHEET,
                   source_row=f"SDT3:{r3 or ''} SDT4:{r4 or ''}",
                   alone_sheet=V_SHEET, alone_source_row=r3 or "")
        ph3 = v_value(b3, "SDT3", r3, "F") if r3 else None
        ph4 = v_value(b4, "SDT4", r4, "F") if r4 else None
        if r3 is None or r4 is None:
            rec["decision"] = "excluded:V1_not_in_both_tables" + ("(alone_only)" if r4 is None else "(mixed_only)")
        elif not ONE_LETTER_SUB.match(label):
            rec["decision"] = "excluded:V2_not_single_missense"
        else:
            x, sx = v_value(b3, "SDT3", r3, "mean"), v_value(b3, "SDT3", r3, "sem")
            y, sy = v_value(b4, "SDT4", r4, "mean"), v_value(b4, "SDT4", r4, "sem")
            if not all(is_num(t) for t in (x, sx, y, sy)):
                rec["decision"] = "excluded:V3_nonnumeric_dV_or_SEM(nd)"
            elif not (isinstance(ph3, str) and ph3.strip() in EPILEPSY
                      and isinstance(ph4, str) and ph4.strip() in EPILEPSY):
                rec["decision"] = f"excluded:V4_phenotype_not_epilepsy({ph3}/{ph4})"
            elif label in CONSTRUCT_EXCEPTIONS:
                rec["decision"] = "excluded:V5_alone_arm_in_CHO-K1+KCNQ3"
            else:
                rec["decision"] = "included"
                keep.append(dict(
                    variant=label, pos=v_value(b3, "SDT3", r3, "A"), phenotype_alone=ph3.strip(),
                    alone_source_row=r3, x=x, sx=sx, alone_gating_cells=v_value(b3, "SDT3", r3, "n"),
                    pos_mixed=v_value(b4, "SDT4", r4, "A"), phenotype_mixed=ph4.strip(),
                    mixed_source_row=r4, y=y, sy=sy,
                    mixed_reported_cells=v_value(b4, "SDT4", r4, "n"),
                    alone_source_file=RECORDED["SDT3"], mixed_source_file=RECORDED["SDT4"], source_sheet=V_SHEET,
                    x_column=V_LAYOUT["SDT3"]["mean"], sx_column=V_LAYOUT["SDT3"]["sem"],
                    y_column=V_LAYOUT["SDT4"]["mean"], sy_column=V_LAYOUT["SDT4"]["sem"]))
        rec["phenotype_alone"], rec["phenotype_mixed"] = ph3, ph4
        log.append(rec)
    return pd.DataFrame(keep), log


# =========================================================================== comparison
def compare(panel, rebuilt, extract, num_fields, str_fields, trace):
    """Field-by-field comparison keyed by variant. `trace(row)` checks recorded coordinates."""
    reb = rebuilt.set_index("variant", drop=False)
    rows, mismatches = [], []
    stats = dict(fields_compared=0, numeric_compared=0, numeric_exact=0, max_rel_diff=0.0)
    for i, e in extract.reset_index(drop=True).iterrows():
        v = e["variant"]
        coord_ok, coord_detail = trace(e)
        if v not in reb.index:
            rows.append((i, e, None, coord_ok, coord_detail, 0, ["not produced by the inclusion rules"]))
            continue
        r = reb.loc[v]
        bad = []
        for f in num_fields:
            stats["fields_compared"] += 1
            stats["numeric_compared"] += 1
            a, b = e[f], r[f]
            d = rel_diff(a, b)
            stats["max_rel_diff"] = max(stats["max_rel_diff"], d)
            if (nan_like(a) and nan_like(b)) or (not nan_like(a) and not nan_like(b) and float(a) == float(b)):
                stats["numeric_exact"] += 1
            if not same_num(a, b):
                bad.append(f"{f}: extract={a!r} rederived={b!r} rel_diff={d:.3g}")
        for f in str_fields:
            stats["fields_compared"] += 1
            if str(e[f]).strip() != str(r[f]).strip():
                bad.append(f"{f}: extract={e[f]!r} rederived={r[f]!r}")
        if not coord_ok:
            bad.append("recorded coordinates: " + coord_detail)
        mismatches += [f"{panel} {v}: {b}" for b in bad]
        rows.append((i, e, r, coord_ok, coord_detail, len(num_fields) + len(str_fields), bad))
    return rows, mismatches, stats


def main():
    say("REDERIVE_KV7 python", sys.version.split()[0], "openpyxl", openpyxl.__version__,
        "pandas", pd.__version__, "numpy", np.__version__)
    for path, expected in EXPECTED_SHA256.items():
        if sha256(path) != expected:
            sys.exit(f"vendored source {path.name} does not match the published checksum")
    hashes_before = {str(p.relative_to(HERE)): sha256(p) for p in INPUTS}
    for p, h in hashes_before.items():
        say("INPUT sha256", h, p)

    results = dict(inputs_sha256=hashes_before, rtol=RTOL)

    # ------------------------------------------------------------------ KCNQ4
    zb = Book(ZHENG)
    kp = k4_parse(zb)
    results["kcnq4_source"] = {}
    for letter, sh in kp.items():
        say(f"SOURCE KCNQ4 {sh['sheet']}: variant_rows={len(sh['rows'])} condition={sh['condition']!r} "
            f"non_variant_rows={sh['non_variant']} title={sh['title']!r}")
        results["kcnq4_source"][sh["sheet"]] = dict(variant_rows=len(sh["rows"]), title=sh["title"],
                                                     condition=sh["condition"],
                                                     non_variant_rows=sh["non_variant"])
    overlap = {f"{a}&{b}": sorted(set(kp[a]["rows"]) & set(kp[b]["rows"]))
               for a in "CDEF" for b in "CDEF" if a < b}
    say("SOURCE KCNQ4 overlap between co-expression sheets:",
        {k: len(v) for k, v in overlap.items()})
    b_only = set(kp["B"]["rows"]) - set().union(*(kp[s]["rows"] for s in "CDEF"))
    say(f"SOURCE KCNQ4 Dataset_B rows with no co-expression row in C-F (never co-expressed): {len(b_only)}")
    k4, k4_log = k4_rederive(zb, kp)
    k4_counts = pd.Series([r["decision"] for r in k4_log]).value_counts().to_dict()
    per_sheet = pd.DataFrame(k4_log).groupby(["sheet", "decision"]).size().to_dict()
    say("RULES KCNQ4 (sequential, first failing rule):", k4_counts)
    for (s, d), n in sorted(per_sheet.items()):
        say(f"RULES KCNQ4   {s:24s} {d:45s} {n}")
    # what rule K3 selects on: homozygous current (Dataset A, column B) of the Dataset C rows
    cur = {v: zb.cell(kp["A"]["sheet"], r, col("B")) for v, r in kp["A"]["rows"].items()}
    c_in = [cur[r["variant"]] for r in k4_log if r["sheet"].endswith("_C") and r["decision"] == "included"]
    c_out = [cur[r["variant"]] for r in k4_log if r["sheet"].endswith("_C") and r["decision"].startswith("excluded:K3")]
    k3_mech = dict(included_min_current=min(c_in), included_median_current=float(np.median(c_in)),
                   excluded_max_current=max(c_out), excluded_median_current=float(np.median(c_out)))
    say("RULES KCNQ4 K3 mechanism (Dataset A Imut/Iwt of Dataset C rows):", k3_mech)
    results["kcnq4_rules"] = dict(counts=k4_counts,
                                  per_sheet={f"{s}|{d}": n for (s, d), n in per_sheet.items()},
                                  dataset_B_never_coexpressed=len(b_only), k3_current_mechanism=k3_mech)

    # ------------------------------------------------------------------ Vanoye
    b3, b4 = Book(SDT3), Book(SDT4)
    p3, p4 = v_parse(b3, "SDT3"), v_parse(b4, "SDT4")
    results["vanoye_source"] = {}
    for key, p in (("SDT3", p3), ("SDT4", p4)):
        say(f"SOURCE Vanoye {key} '{V_SHEET}': variant_rows={len(p['rows'])} title={p['title']!r} "
            f"condition={p['condition']!r} wt_rows={p['wt_rows']}")
        for a in p["anomalies"]:
            say("SOURCE Vanoye label/column anomaly:", a)
        results["vanoye_source"][key] = dict(variant_rows=len(p["rows"]), title=p["title"],
                                              condition=p["condition"], wt_rows=p["wt_rows"],
                                              anomalies=p["anomalies"])
    v5_found = {}
    texts = {"repo_cached_paper_text": VANOYE_TEXT,
             "fresh_publisher_download": OUT / "downloads/jci/article_156314.html"}
    for k, path in texts.items():
        if path.exists():  # read as text only; downloaded files are data, never executed
            text = path.read_text(errors="replace")
            if path.suffix == ".html":
                text = html.unescape(re.sub(r"<[^>]+>", " ", re.sub(r"<script.*?</script>", " ", text, flags=re.S)))
            v5_found[k] = bool(V5_SENTENCE.search(re.sub(r"\s+", " ", text)))
    say("RULE V5 source sentence found in:", v5_found)
    vq, vq_log = v_rederive(b3, b4, p3, p4)
    vq_counts = pd.Series([r["decision"] for r in vq_log]).value_counts().to_dict()
    say("RULES Vanoye (sequential, first failing rule):", vq_counts)
    for r in vq_log:
        if r["decision"] != "included":
            say(f"RULES Vanoye   {r['variant']:12s} {r['source_row']:18s} "
                f"{str(r['phenotype_alone']) + ' | ' + str(r['phenotype_mixed']):22s} {r['decision']}")
    discordant = [(r["variant"], r["phenotype_alone"], r["phenotype_mixed"]) for r in vq_log
                  if r["phenotype_alone"] is not None and r["phenotype_mixed"] is not None
                  and str(r["phenotype_alone"]).strip() != str(r["phenotype_mixed"]).strip()]
    say("RULES Vanoye phenotype label differs between SDT3 and SDT4:", discordant)
    # what rule V3 selects on: homozygous current (SDT3 column I, fraction of WT) of the rows
    i3 = {r["variant"]: b3.cell(V_SHEET, p3["rows"][r["variant"]], col("I")) for r in vq_log
          if r["variant"] in p3["rows"]}
    v_in = [i3[r["variant"]] for r in vq_log if r["decision"] in ("included", "excluded:V5_alone_arm_in_CHO-K1+KCNQ3")]
    v_out = [i3[r["variant"]] for r in vq_log if r["decision"].startswith("excluded:V3")]
    v3_mech = dict(measured_n=len(v_in), measured_median_current=float(np.median(v_in)),
                   measured_min_current=min(v_in), nd_n=len(v_out),
                   nd_median_current=float(np.median(v_out)), nd_max_current=max(v_out))
    say("RULES Vanoye V3 mechanism (SDT3 col I, epilepsy rows measured vs nd):", v3_mech)
    results["vanoye_rules"] = dict(counts=vq_counts, v5_sentence_found=v5_found,
                                   phenotype_discordant=discordant, v3_current_mechanism=v3_mech)

    # ------------------------------------------------------------------ write rebuilt tables
    # round_trip parsing turns each decimal string back into the exact double it was written
    # from; pandas' default fast parser can be off by an ulp or so (run_kv7.py uses the default,
    # so task 4 below re-reads the extracts exactly as run_kv7.py does).
    e4 = pd.read_csv(EXT_K4, float_precision="round_trip")
    ev = pd.read_csv(EXT_VQ, float_precision="round_trip")
    k4 = k4[[c for c in e4.columns if c in k4.columns]]
    vq = vq[[c for c in ev.columns if c in vq.columns]]
    pd.DataFrame(k4_log + vq_log).to_csv(OUT / "kv7_source_rule_log.csv", index=False)

    # ------------------------------------------------------------------ compare
    def trace_k4(e):
        problems = []
        checks = [(e["alone_sheet"], e["alone_source_row"]), (e["mixed_sheet"], e["mixed_source_row"]),
                  (K4_SHEETS["A"], e["alone_current_source_row"])]
        for sheet, row in checks:
            got = zb.cell(sheet, int(row), 1)
            if (got or "").strip() != e["variant"]:
                problems.append(f"{sheet}!A{int(row)}={got!r}")
        if e["mixed_sheet"] != K4_SHEETS.get(e["sheet"]):
            problems.append(f"sheet letter {e['sheet']} vs {e['mixed_sheet']}")
        if e["source_file"] != RECORDED["ZHENG"]:
            problems.append(f"source_file {e['source_file']}")
        return (not problems), "; ".join(problems) or "recorded rows hold this variant"

    def trace_vq(e):
        problems = []
        for book, row, path, key in ((b3, e["alone_source_row"], e["alone_source_file"], SDT3),
                                     (b4, e["mixed_source_row"], e["mixed_source_file"], SDT4)):
            got = book.cell(e["source_sheet"], int(row), col("D"))
            if (got or "").strip() != e["variant"]:
                problems.append(f"{key.name}!D{int(row)}={got!r}")
            if path != RECORDED["SDT3" if key == SDT3 else "SDT4"]:
                problems.append(f"file {path}")
        return (not problems), "; ".join(problems) or "recorded rows hold this variant"

    k4_num = ["x", "sx", "y", "sy", "mixed_reported_cells", "alone_current_cells", "alone_gating_cells",
              "pos", "mixed_source_row", "alone_source_row", "alone_current_source_row"]
    k4_str = ["sheet", "alone_sheet", "mixed_sheet", "source_file", "x_column", "sx_column",
              "y_column", "sy_column", "primary"]
    vq_num = ["x", "sx", "y", "sy", "alone_gating_cells", "mixed_reported_cells", "pos", "pos_mixed",
              "alone_source_row", "mixed_source_row"]
    vq_str = ["phenotype_alone", "phenotype_mixed", "alone_source_file", "mixed_source_file",
              "source_sheet", "x_column", "sx_column", "y_column", "sy_column"]
    k4_rows, k4_mm, k4_stats = compare("KCNQ4", k4, e4, k4_num, k4_str, trace_k4)
    vq_rows, vq_mm, vq_stats = compare("KCNQ2_KCNQ3", vq, ev, vq_num, vq_str, trace_vq)
    # the same comparison with pandas' default float parser, i.e. the values run_kv7.py sees
    _, k4_mm_d, k4_stats_d = compare("KCNQ4", k4, pd.read_csv(EXT_K4), k4_num, k4_str, trace_k4)
    _, vq_mm_d, vq_stats_d = compare("KCNQ2_KCNQ3", vq, pd.read_csv(EXT_VQ), vq_num, vq_str, trace_vq)
    for panel, st, mm in (("KCNQ4", k4_stats_d, k4_mm_d), ("KCNQ2_KCNQ3", vq_stats_d, vq_mm_d)):
        say(f"COMPARE {panel} (pandas default float parser, as in run_kv7.py): "
            f"numeric_bitwise_equal={st['numeric_exact']}/{st['numeric_compared']} "
            f"max_rel_diff={st['max_rel_diff']:.3g} mismatches_at_rtol={len(mm)}")

    exhaust = {}
    for panel, reb, ext, log in (("KCNQ4", k4, e4, k4_log), ("KCNQ2_KCNQ3", vq, ev, vq_log)):
        missing = sorted(set(reb.variant) - set(ext.variant))
        extra = sorted(set(ext.variant) - set(reb.variant))
        why = {r["variant"]: r["decision"] for r in log}
        exhaust[panel] = dict(source_rows_satisfying_rules=len(reb), extract_rows=len(ext),
                              extract_duplicates=int(ext.variant.duplicated().sum()),
                              missing_from_extract=missing,
                              extract_rows_not_produced_by_rules={v: why.get(v, "not found in any source sheet")
                                                                  for v in extra})
        say(f"EXHAUSTIVE {panel}: rule-satisfying source rows={len(reb)} extract rows={len(ext)} "
            f"missing_from_extract={missing} extract_not_from_rules={exhaust[panel]['extract_rows_not_produced_by_rules']}")
    for panel, st, mm in (("KCNQ4", k4_stats, k4_mm), ("KCNQ2_KCNQ3", vq_stats, vq_mm)):
        say(f"COMPARE {panel}: fields_compared={st['fields_compared']} numeric={st['numeric_compared']} "
            f"numeric_bitwise_equal={st['numeric_exact']} max_rel_diff={st['max_rel_diff']:.3g} "
            f"mismatches={len(mm)}")
        for m in mm:
            say("MISMATCH", m)
    results["compare"] = dict(KCNQ4=dict(stats=k4_stats, mismatches=k4_mm,
                                         default_parser_stats=k4_stats_d, default_parser_mismatches=k4_mm_d),
                              KCNQ2_KCNQ3=dict(stats=vq_stats, mismatches=vq_mm,
                                               default_parser_stats=vq_stats_d,
                                               default_parser_mismatches=vq_mm_d))
    results["exhaustiveness"] = exhaust

    # ------------------------------------------------------------------ row map
    with open(OUT / "kv7_row_map.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["extract_file", "extract_line", "variant", "alone_source_file", "alone_sheet",
                    "alone_source_row", "mixed_source_file", "mixed_sheet", "mixed_source_row",
                    "cells_source_rows", "coordinates_check", "fields_compared", "fields_mismatched",
                    "match_status", "detail"])
        for i, e, r, ok, cdet, nf, bad in k4_rows:
            w.writerow([rel(EXT_K4), i + 2, e["variant"], rel(ZHENG), e["alone_sheet"], int(e["alone_source_row"]),
                        rel(ZHENG), e["mixed_sheet"], int(e["mixed_source_row"]),
                        f"mixed n: {e['mixed_sheet']}!D{int(e['mixed_source_row'])}; alone n (current): "
                        f"Supplemental_Dataset_A!D{int(e['alone_current_source_row'])}; alone gating n: not reported",
                        "OK" if ok else "FAIL: " + cdet, nf, len(bad),
                        "MATCH" if (r is not None and not bad) else ("NOT_REDERIVED" if r is None else "MISMATCH"),
                        "; ".join(bad)])
        for i, e, r, ok, cdet, nf, bad in vq_rows:
            w.writerow([rel(EXT_VQ), i + 2, e["variant"], rel(SDT3), V_SHEET, int(e["alone_source_row"]),
                        rel(SDT4), V_SHEET, int(e["mixed_source_row"]),
                        f"alone n: SDT3!V{int(e['alone_source_row'])}; mixed n: SDT4!W{int(e['mixed_source_row'])}",
                        "OK" if ok else "FAIL: " + cdet, nf, len(bad),
                        "MATCH" if (r is not None and not bad) else ("NOT_REDERIVED" if r is None else "MISMATCH"),
                        "; ".join(bad)])

    # ------------------------------------------------------------------ task 4: run_kv7 panel
    src = RUN_KV7.read_text()
    required = ['d = d[(d.sx > 0) & (d.sy > 0)].copy()', 'cls = d.variant.map(D.classify_variant)',
                'd["klass"] = [c[3] for c in cls]', 'mis = d[d.klass == "missense"]',
                'r = D.estimate(mis, units="physical")', 'kv7/fit_KCNQ4_primary.csv',
                'kv7/fit_KCNQ2_KCNQ3_primary.csv']
    missing_stmts = [s for s in required if s not in src]
    say("RUN_KV7 filter statements present in run_kv7.py:", not missing_stmts, missing_stmts)
    D = types.ModuleType("dominance_readonly")
    D.__file__ = str(DOMINANCE)
    sys.modules["dominance_readonly"] = D  # dataclasses resolves annotations through sys.modules
    exec(compile(DOMINANCE.read_text(), str(DOMINANCE), "exec"), D.__dict__)

    def run_kv7_panel(df):
        d = df[(df.sx > 0) & (df.sy > 0)].copy()
        dropped = sorted(set(df.variant) - set(d.variant))
        cls = d.variant.map(D.classify_variant)
        d["klass"] = [c[3] for c in cls]
        mis = d[d.klass == "missense"]
        r = D.estimate(mis, units="physical")
        return dict(extract_rows=len(df), dropped_by_sx_sy=dropped,
                    classes={k: int(v) for k, v in d.klass.value_counts().items()},
                    n=r.n, residues=r.n_residues, b=r.b, ci=list(r.ci), pearson_r=r.pearson_r,
                    interpretable=r.interpretable, warnings=r.warnings_)

    results["panel"] = {}
    for key, ext, reb in (("KCNQ4", pd.read_csv(EXT_K4), k4), ("KCNQ2_KCNQ3", pd.read_csv(EXT_VQ), vq)):
        out_e = run_kv7_panel(ext)
        out_r = run_kv7_panel(reb)
        rep = REPORTED[key]
        agree = (out_e["n"] == rep["n"] and out_e["residues"] == rep["residues"]
                 and round(out_e["b"], 2) == rep["b"]
                 and tuple(round(c, 2) for c in out_e["ci"]) == rep["ci"])
        say(f"PANEL {key} (vendored extract): rows={out_e['extract_rows']} dropped_by_sx/sy>0={out_e['dropped_by_sx_sy']} "
            f"classes={out_e['classes']} n={out_e['n']} residues={out_e['residues']} b={out_e['b']:.6f} "
            f"CI=[{out_e['ci'][0]:.6f}, {out_e['ci'][1]:.6f}] r={out_e['pearson_r']:.4f}")
        say(f"PANEL {key} (rederived table): n={out_r['n']} residues={out_r['residues']} b={out_r['b']:.6f} "
            f"CI=[{out_r['ci'][0]:.6f}, {out_r['ci'][1]:.6f}]")
        say(f"PANEL {key} manuscript: b={rep['b']} n={rep['n']} residues={rep['residues']} CI={rep['ci']} "
            f"-> reproduced at reported precision: {agree}")
        results["panel"][key] = dict(vendored=out_e, rederived=out_r, manuscript=rep, reproduced=agree)

    # ------------------------------------------------------------------ within-workbook consistency
    # Not part of the extraction: checks whether the sheet each value was taken from agrees with
    # the other sheets of the same workbook for the included variants.
    consistency = {}
    # KCNQ4: Dataset B repeats the homozygous current of Dataset A; equal values show that the B
    # row (x, sx) and the A row (alone_current_cells) describe the same variant series.
    same_cur = sum(1 for v in k4.variant
                   if zb.cell(kp["B"]["sheet"], kp["B"]["rows"][v], col("B"))
                   == zb.cell(kp["A"]["sheet"], kp["A"]["rows"][v], col("B")))
    consistency["KCNQ4_B_vs_A_current_equal"] = f"{same_cur}/{len(k4)}"
    say(f"CONSISTENCY KCNQ4 Dataset_B current == Dataset_A current for included variants: {same_cur}/{len(k4)}")
    # Vanoye: 'Relative to WT channel' (dV1/2 vs batch WT) against 'Variant mean values ' (absolute V1/2).
    for key, book, p, mean_sheet, cols in (("SDT3", b3, p3, "Variant mean values ", ("N", "O", "V")),
                                           ("SDT4", b4, p4, "Variant mean values ", ("O", "P", "W"))):
        ms = book.v[mean_sheet]
        mrows = {}
        for r in range(13, ms.max_row + 1):
            lab = book.cell(mean_sheet, r, 4)
            if isinstance(lab, str) and lab.strip():
                mrows[lab.strip()] = r
        nd_abs, sem_diff, n_diff, implied_wt = [], [], [], []
        for v in vq.variant:
            r, m = p["rows"][v], mrows[v]
            dv, dsem, dn = (book.cell(V_SHEET, r, col(c)) for c in cols)
            av, asem, an = (book.cell(mean_sheet, m, col(c)) for c in ("N", "O", "T"))
            if not is_num(av):
                nd_abs.append(f"{v} (abs V1/2 {av!r})")
                continue
            implied_wt.append(av - dv)
            if not math.isclose(dsem, asem, rel_tol=1e-6):
                sem_diff.append(f"{v} ({dsem:.4f} vs {asem:.4f})")
            if dn != an:
                n_diff.append(f"{v} ({dn} vs {an})")
        wt_avg = [book.cell(mean_sheet, r, col("N")) for r in range(13, 17)
                  if str(book.cell(mean_sheet, r, 4) or "").lower().startswith("wt")]
        c = dict(abs_V_nd_but_dV_given=nd_abs, sem_differs=sem_diff, n_differs=n_diff,
                 implied_batch_WT_V12_min=min(implied_wt), implied_batch_WT_V12_median=float(np.median(implied_wt)),
                 implied_batch_WT_V12_max=max(implied_wt), study_average_WT_V12=wt_avg)
        consistency[f"Vanoye_{key}"] = c
        say(f"CONSISTENCY Vanoye {key} included variants, 'Relative to WT channel' vs 'Variant mean values': "
            f"abs V1/2 'nd' but dV1/2 given={nd_abs}; SEM differs={len(sem_diff)} {sem_diff}; n differs={n_diff}; "
            f"implied batch-WT V1/2 (abs - dV) min/median/max = {min(implied_wt):.2f}/{np.median(implied_wt):.2f}/"
            f"{max(implied_wt):.2f} mV vs study-average WT {wt_avg}")
    results["within_workbook_consistency"] = consistency

    # ------------------------------------------------------------------ context: rule sensitivity
    def point(df):
        d = df[(df.sx > 0) & (df.sy > 0)].copy()
        d = d[[D.classify_variant(v)[3] == "missense" for v in d.variant]]
        return dict(n=len(d), residues=int(d.pos.nunique()), b=D.slope_through_origin(d.x, d.y))

    sens = {}
    allk = []  # KCNQ4 with all four co-expression sheets C, D, E, F (K1 dropped)
    for letter in "CDEF":
        sh = kp[letter]
        for v, row in sh["rows"].items():
            br = kp["B"]["rows"].get(v)
            if br is None or k4_hgvs(v) is None:
                continue
            vals = [zb.cell(kp["B"]["sheet"], br, col("D")), zb.cell(kp["B"]["sheet"], br, col("E")),
                    zb.cell(sh["sheet"], row, col("E")), zb.cell(sh["sheet"], row, col("F"))]
            if all(is_num(t) for t in vals):
                allk.append(dict(variant=v, pos=k4_hgvs(v)[1], x=vals[0], sx=vals[1], y=vals[2], sy=vals[3]))
    sens["KCNQ4 C+E (rules K1-K4)"] = point(k4)
    sens["KCNQ4 C+D+E+F (K1 dropped)"] = point(pd.DataFrame(allk))
    vall = []
    for label, r3 in p3["rows"].items():
        r4 = p4["rows"].get(label)
        if r4 is None or not ONE_LETTER_SUB.match(label):
            continue
        vals = [v_value(b3, "SDT3", r3, "mean"), v_value(b3, "SDT3", r3, "sem"),
                v_value(b4, "SDT4", r4, "mean"), v_value(b4, "SDT4", r4, "sem")]
        if all(is_num(t) for t in vals):
            vall.append(dict(variant=label, pos=int(v_value(b3, "SDT3", r3, "A")), x=vals[0], sx=vals[1],
                             y=vals[2], sy=vals[3], ph=f"{v_value(b3, 'SDT3', r3, 'F')}",
                             ph4=f"{v_value(b4, 'SDT4', r4, 'F')}"))
    vall = pd.DataFrame(vall)
    epi = vall[vall.ph.isin(EPILEPSY) & vall.ph4.isin(EPILEPSY)]
    sens["Vanoye V1-V5"] = point(vq)
    sens["Vanoye epilepsy incl. A193D/R201C/P335L (V5 dropped)"] = point(epi)
    sens["Vanoye all measurable pairs (V4, V5 dropped)"] = point(vall)
    nd3 = [s.split(" ")[0] for s in consistency["Vanoye_SDT3"]["abs_V_nd_but_dV_given"]]
    sens[f"Vanoye V1-V5 without {'/'.join(nd3)} (SDT3 absolute V1/2 'nd')"] = point(vq[~vq.variant.isin(nd3)])
    for k, v in sens.items():
        say(f"CONTEXT point slope {k}: n={v['n']} residues={v['residues']} b={v['b']:.4f}")
    results["context_point_slopes"] = sens

    # ------------------------------------------------------------------ finish
    hashes_after = {str(p.relative_to(HERE)): sha256(p) for p in INPUTS}
    results["inputs_unchanged_during_run"] = hashes_after == hashes_before
    say("INPUTS unchanged during run:", results["inputs_unchanged_during_run"])
    (OUT / "kv7_provenance.json").write_text(
        json.dumps(results, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    (OUT / "kv7_provenance.log").write_text("\n".join(LOG) + "\n")
    ok = (not k4_mm and not vq_mm and all(not v["missing_from_extract"]
                                           and not v["extract_rows_not_produced_by_rules"]
                                           for v in exhaust.values())
          and all(p["reproduced"] for p in results["panel"].values()))
    say("OVERALL", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
