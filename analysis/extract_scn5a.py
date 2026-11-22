

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Extract the SCN5A peak-current panel from the published supplement.

    python3 extract_scn5a.py /path/to/NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx

Source: O'Neill MJ, Muhammad A, Li B, et al. "Dominant negative effects of
SCN5A missense variants." Genet Med. 2022;24(6):1238-1248.  PMC9262418.
Supplementary workbook, sheet "Peak Current": peak sodium current density as
a percentage of wild type, measured alone ("Homozygous") and co-expressed with
wild type ("Heterozygous"), with standard errors and cell counts.

Writes scn5a_peak.csv (50 missense variants and W822X) and scn5a_controls.csv (WT and
WT+WT).  Both are committed, so this script is a provenance record and a
re-extraction check rather than a required step.
"""
import hashlib
import sys

import pandas as pd

SHEET = "Peak Current"
COLS = ["variant", "homo", "homo_se", "homo_n", "het", "het_se", "het_n"]

src = sys.argv[1] if len(sys.argv) > 1 else \
    _data("NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx", "../NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx")

print(f"source: {src}")
print(f"sha256: {hashlib.sha256(open(src,'rb').read()).hexdigest()}")

raw = pd.ExcelFile(src).parse(SHEET, header=0)
raw.columns = COLS
raw = raw[raw.variant.notna()].copy()

# The WT+WT row carries "-" in the alone columns: it has no alone state.
for c in COLS[1:]:
    raw[c] = pd.to_numeric(raw[c], errors="coerce")

ctl = raw[raw.variant.isin(["WT", "WT+WT"])]
var = raw[~raw.variant.isin(["WT", "WT+WT"])]
var = var[var.het.notna() & var.homo.notna()]

ctl.to_csv("scn5a_controls.csv", index=False)
var.to_csv("scn5a_peak.csv", index=False)
print(f"\nwrote scn5a_controls.csv ({len(ctl)} rows) and "
      f"scn5a_peak.csv ({len(var)} variants)")
print(ctl.to_string(index=False))
