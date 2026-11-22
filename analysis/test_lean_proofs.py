"""Checks on the Lean proofs in lean/ (revision of 8 November 2026).

The source checks always run. The build and axiom audit need the Lean
toolchain and the downloaded Mathlib library; without them those checks skip,
so the offline reproduction run is unaffected.
"""
from pathlib import Path
import os
import re
import shutil
import subprocess

import pytest

HERE = Path(__file__).resolve().parent
LEAN = HERE / "lean"
SOURCES = sorted((LEAN / "DominanceProofs").glob("*.lean"))
STANDARD = "[propext, Classical.choice, Quot.sound]"


def supplement_text():
    for candidate in (HERE / "submission-2026-09-29" / "supplement_source.md",
                      HERE.parent / "manuscript" / "supplement_source.md"):
        if candidate.is_file():
            return candidate.read_text()
    pytest.skip("supplement source not found")


def public_theorems():
    names = []
    for path in SOURCES:
        names += re.findall(r"^theorem\s+(\S+)", path.read_text(), re.M)
    return names


def test_lean_sources_contain_no_unproven_steps():
    assert len(SOURCES) == 7
    for path in SOURCES:
        text = re.sub(r"/-.*?-/", "", path.read_text(), flags=re.S)
        text = re.sub(r"--.*", "", text)
        for token in ("sorry", "admit", "native_decide", "implemented_by", "unsafe", "debug.skipKernelTC"):
            assert not re.search(r"(?<![A-Za-z_])" + re.escape(token) + r"(?![A-Za-z_])", text), (path.name, token)
        assert not re.search(r"^\s*axiom\s", text, re.M), path.name


def test_audit_lists_every_public_theorem():
    audited = re.findall(r"#print axioms Dominance\.(\S+)", (LEAN / "Audit.lean").read_text())
    assert sorted(audited) == sorted(public_theorems())
    assert len(audited) == 100


def test_supplement_theorem_map_names_real_theorems():
    text = supplement_text()
    section = text.split("## S13.", 1)[1].split("\n## ", 1)[0]
    cited = set(re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", section))
    theorems = set(public_theorems())
    files = {p.stem for p in SOURCES}
    unknown = sorted(c for c in cited if c not in theorems and c not in files and not c.endswith("lean"))
    assert not unknown, unknown
    assert len(cited & theorems) >= 40


def _lake():
    found = shutil.which("lake") or str(Path.home() / ".elan" / "bin" / "lake")
    return found if Path(found).exists() else None


@pytest.mark.skipif(not (LEAN / ".lake" / "packages" / "mathlib").is_dir(),
                    reason="Mathlib is not downloaded; run `lake exe cache get` in lean/")
def test_lean_build_and_axioms():
    lake = _lake()
    if lake is None:
        pytest.skip("Lean toolchain not installed")
    env = dict(os.environ, PATH=str(Path(lake).parent) + os.pathsep + os.environ.get("PATH", ""))
    build = subprocess.run([lake, "build"], cwd=LEAN, env=env, capture_output=True, text=True, timeout=3600)
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]
    audit = subprocess.run([lake, "env", "lean", "Audit.lean"], cwd=LEAN, env=env,
                           capture_output=True, text=True, timeout=3600)
    assert audit.returncode == 0, audit.stderr[-3000:]
    lines = [line for line in audit.stdout.splitlines() if "depends on axioms" in line]
    assert len(lines) == len(public_theorems())
    assert all(line.endswith(STANDARD) for line in lines), [l for l in lines if not l.endswith(STANDARD)]
