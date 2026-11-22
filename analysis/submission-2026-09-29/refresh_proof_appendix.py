#!/usr/bin/env python3
"""Refresh the proof inventory and headline signatures from the checked Lean sources."""
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
START = "<!-- BEGIN GENERATED LEAN STATEMENTS -->"
END = "<!-- END GENERATED LEAN STATEMENTS -->"
HEADLINES = ("identifies_iff", "configMean_congr", "assembly_floor", "sharp_marginal_bounds")


def main():
    files = sorted((ROOT / "lean" / "DominanceProofs").glob("*.lean"))
    sources = {p.name: p.read_text() for p in files}
    inventory = {name: re.findall(r"^theorem\s+(\w+)", source, re.M)
                 for name, source in sources.items()}
    count = sum(map(len, inventory.values()))
    statements = {}
    for name, source in sources.items():
        for theorem in HEADLINES:
            match = re.search(r"^theorem " + theorem + r"\b.*?(?= := by)", source, re.M | re.S)
            if match:
                statements[theorem] = (name, match.group())
    if set(statements) != set(HEADLINES):
        raise ValueError("A headline theorem was not found")
    text = (HERE / "supplement_source.md").read_text()
    text = re.sub(r"formalization comprises \d+ theorems in seven files",
                  f"formalization comprises {count} theorems in seven files", text)
    # Keep the numerical example in Table O synchronized with its formal statement.
    example = re.search(r"theorem dimer_family_crosses_floor.*?\(1\.338 \* ([\d.]+)\)",
                        sources["Arrangements.lean"], re.S)
    if not example:
        raise ValueError("Cannot locate the nonrandom-assembly example")
    text = re.sub(r"reproducing 24.7%, below the floor, at \$h=[\d.]+\$",
                  f"reproducing 24.7%, below the floor, at $h={example.group(1)}$", text)
    generated = [START, "### Source-derived theorem inventory and headline statements", "",
                 f"The source inventory contains {count} public theorems. The following counts and formal statements are regenerated from the Lean source files at each manuscript build; Table O maps the scientific claims to those names.", "",
                 "**Table O1. Public theorem inventory by Lean source file.**", "",
                 "| Lean source file | Public theorems |", "|:--|--:|"]
    generated += [f"| {name} | {len(names)} |" for name, names in inventory.items()]
    generated += ["", "The signatures below are exact source statements, with proof bodies omitted. In these statements, Identifies fixes the endpoints and asks whether the interior count profile is unique, interiorCount counts distinct fractions strictly between zero and one, and meanActivity is the binomial activity mixture. countAverage averages activities over configurations with the same mutant count. MonotoneProfile, Lbound, and Ubound encode the profile restriction and bounds of Section S9. The biological and measurement assumptions needed to apply these mathematical objects remain those in the main text.", ""]
    for theorem in HEADLINES:
        filename, signature = statements[theorem]
        generated += [f"{theorem} ({filename}):", "", "```lean", signature, "```", ""]
    generated += [END]
    block = "\n".join(generated)
    if START in text:
        text = re.sub(re.escape(START) + r".*?" + re.escape(END), lambda _: block, text, flags=re.S)
    else:
        text = text.replace("\n## References\n", "\n" + block + "\n\n## References\n")
    (HERE / "supplement_source.md").write_text(text)
    return count


if __name__ == "__main__":
    print(f"Refreshed {main()} public theorem declarations")
