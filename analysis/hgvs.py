"""HGVS protein names for SCN5A variants (revision of 8 November 2026).

Source labels in the O'Neill et al. workbook follow NP_932173.1, which carries
one extra glutamine at position 1077 relative to the MANE Select reference
NP_000326.2. Positions above 1077 therefore shift down by one. Every converted
name is checked against the vendored NP_000326.2 sequence.
"""
import re
from functools import lru_cache
from pathlib import Path

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
       "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
       "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "X": "Ter", "*": "Ter"}
REFSEQ = Path(__file__).resolve().parent / "data" / "refseq" / "NP_000326.2.fasta"


@lru_cache(maxsize=1)
def _sequence():
    return "".join(l.strip() for l in REFSEQ.read_text().splitlines() if not l.startswith(">"))


def scn5a_hgvs(label, accession=False):
    """Convert a source label such as G1712C or W822X to p.Gly1711Cys or p.Trp822Ter."""
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z*])", label)
    if not m:
        raise ValueError(f"unrecognised source label {label!r}")
    wt, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if pos == 1077:
        raise ValueError("position 1077 has no NP_000326.2 counterpart")
    mane = pos - 1 if pos > 1077 else pos
    if _sequence()[mane - 1] != wt:
        raise ValueError(f"{label}: NP_000326.2 residue {mane} is {_sequence()[mane - 1]}, not {wt}")
    name = f"p.{AA3[wt]}{mane}{AA3[alt]}"
    return f"NP_000326.2:{name}" if accession else name
