#!/usr/bin/env python3
"""Verify S1 File, restore declared external inputs, then reproduce offline.

Use --skip-profile for a focused revision whose full profile simulation is unchanged.
That mode retains packaged profile outputs and explicitly excludes them from the
reproduced-artifact count. Excluded third-party inputs must be explicitly supplied
with --input-cache (matching paths/checksums) or fetched with --download-inputs.
No undeclared source-tree fallback is used.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import date
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent
TEST_FILES = ["test_dominance.py", "test_profile_revision.py", "test_assay_revision.py",
              "test_bounds.py", "test_tp53_anchors.py", "test_survey_revision.py",
              "test_revision_nov2026.py", "test_lean_proofs.py", "test_kcnq2_archive_nov20.py",
              "test_normalization_titration_revision.py", "test_scn5a_round2.py",
              "test_archive_reading_round2.py", "test_titration_round2.py"]
CSV_OUTPUTS = ["scn5a_peak.csv", "scn5a_controls.csv", "scn5a_specifications.csv",
               "kv7_controls.csv", "kv7_rank.csv", "deconvolution.csv",
               "tp53_portability.csv", "tp53_controls.csv", "tp53_anchors.csv",
               "profile_scn5a.csv", "clinical_null.csv", "clinvar_crosswalk.csv",
               "survey_decisions.csv", "bounds.csv", "tp53_hotspot.csv", "survey_validation.csv",
               "kv7_row_map.csv", "deconvolution_origin_differences.csv",
               "deconvolution_control_means.csv", "deconvolution_se_comparator.csv",
               "scn5a_robustness.csv", "scn5a_joint_uncertainty.csv", "scn5a_reference_dose.csv",
               "survey_second_reading_comparison.csv", "survey_blinded_reading_comparison.csv", "titration_profiles.csv", "titration_candidate_tests.csv",
               "titration_sensitivity_ranges.csv", "titration_fit_sensitivity_ranges.csv", "titration_replication.csv"]
JSON_OUTPUTS = ["scn5a_specifications.json", "kv7_controls.json", "kv7_rank.json",
                "deconvolution.json", "tp53_portability.json", "tp53_anchors.json",
                "profile_identity.json", "clinical_null.json", "clinvar_crosswalk.json",
                "survey_summary.json", "survey_zheng_case.json", "bounds.json",
                "tp53_hotspot.json", "survey_validation.json", "selection.json", "kv7_panels.json",
                "kv7_provenance.json", "rank_null.json", "uncertainty.json", "titration.json",
                "scn5a_floor_scenarios.json"]
FIGURES = [f"{stem}.{suffix}" for stem in ("Schematic", "Figure1", "Figure2", "Figure3", "Titration", "TitrationMain", "SCN5A_robustness") for suffix in ("png", "pdf")]
PROFILE_OUTPUTS = {"profile_identity.json", "profile_scn5a.csv"}


def normalized(value, prefixes):
    if isinstance(value, dict):
        return {k: normalized(v, prefixes) for k, v in value.items() if k != "generated"}
    if isinstance(value, list):
        return [normalized(v, prefixes) for v in value]
    if isinstance(value, str):
        for prefix in sorted(prefixes, key=len, reverse=True):
            value = value.replace(prefix, "<PACKAGE>")
    return value


def differences(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return [f"{path}: keys differ {set(a) ^ set(b)}"]
        return [d for k in a for d in differences(a[k], b[k], f"{path}.{k}")]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: lengths {len(a)} vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, f"{path}[{i}]")]
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return [] if np.isclose(a, b, rtol=1e-12, atol=1e-12, equal_nan=True) else [f"{path}: {a} vs {b}"]
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=HERE / "S1_File.zip")
    parser.add_argument("--report-dir", type=Path, default=HERE)
    parser.add_argument("--input-cache", type=Path,
                        help="Explicit external-input cache mirroring paths in INPUT_RIGHTS.json")
    parser.add_argument("--download-inputs", action="store_true",
                        help="Retrieve declared external inputs before disabling network access")
    parser.add_argument("--metadata-bundle", type=Path,
                        help="Separate original-license historical metadata ZIP, verified against S1 hashes")
    parser.add_argument("--skip-profile", action="store_true",
                        help="Retain packaged profile outputs; do not claim they were regenerated")
    args = parser.parse_args()
    archive_path = args.package.resolve()
    report = args.report_dir.resolve()
    report.mkdir(parents=True, exist_ok=True)
    logs = report / "reproduction_logs"
    logs.mkdir(exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="piezo1-dominance-repro-"))
    unpacked, baseline = stage / "archive", stage / "baseline"
    baseline.mkdir()
    started = time.time()
    archive_sha = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive_path) as archive:
        manifest = json.loads(archive.read("FILE_MANIFEST.json"))["files"]
        for entry in manifest:
            payload = archive.read(entry["path"])
            if len(payload) != entry["bytes"] or hashlib.sha256(payload).hexdigest() != entry["sha256"]:
                raise RuntimeError(f"Archive checksum mismatch: {entry['path']}")
        for name in archive.namelist():
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError(f"Unsafe archive path: {name}")
        archive.extractall(unpacked)
    package = unpacked / "analysis"
    if not package.is_dir():
        raise ValueError("S1 File must contain an analysis directory")
    external_inputs = None
    if (package / "INPUT_RIGHTS.json").exists():
        spec = importlib.util.spec_from_file_location("packaged_input_restore", package / "restore_external_inputs.py")
        restore_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(restore_module)
        external_inputs = restore_module.restore(package, args.input_cache, args.download_inputs, args.metadata_bundle)
    input_manifest = {str(p.relative_to(package)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in sorted(package.rglob("*")) if p.is_file()}
    (report / "reproduction_input_manifest.json").write_text(json.dumps(input_manifest, indent=2) + "\n")
    all_outputs = CSV_OUTPUTS + JSON_OUTPUTS + FIGURES
    skipped_outputs = PROFILE_OUTPUTS if args.skip_profile else set()
    for filename in all_outputs:
        source = package / filename
        if not source.exists():
            raise FileNotFoundError(f"Required baseline is missing from S1 File: {filename}")
        shutil.copy2(source, baseline / filename)
        if filename not in skipped_outputs:
            source.unlink()

    guard = stage / "offline_guard"
    guard.mkdir()
    network_log = stage / "network_attempts.log"
    (guard / "sitecustomize.py").write_text(
        'import os, socket\nfrom pathlib import Path\n'
        'def offline_connect(self, *args, **kwargs):\n'
        '    with Path(os.environ["REPRO_NETWORK_LOG"]).open("a") as fh:\n'
        '        fh.write(repr(args) + "\\n")\n'
        '    raise OSError("Network disabled for self-contained reproduction")\n'
        'socket.socket.connect = offline_connect\n'
        'socket.socket.connect_ex = offline_connect\n')
    env = os.environ.copy()
    env.update(PYTHONPATH=str(guard), REPRO_NETWORK_LOG=str(network_log),
               MPLCONFIGDIR=str(stage / "matplotlib"), PYTHONDONTWRITEBYTECODE="1")
    python = Path(sys.executable)
    commands = {}

    def execute(name, command_args):
        log = logs / f"{name}.log"
        command_started = time.time()
        with log.open("w") as stream:
            result = subprocess.run([str(python), *command_args], cwd=package, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT)
        record = {"command": [str(python), *command_args], "exit_code": result.returncode,
                  "seconds": round(time.time() - command_started, 2), "log": str(log.relative_to(report))}
        print(f"{name}: exit {result.returncode}, {record['seconds']:.2f} s", flush=True)
        if result.returncode:
            print(log.read_text()[-3500:], flush=True)
        return record

    tasks = {
        "extract_scn5a": ([], ["extract_scn5a.py"]),
        "run_scn5a": (["extract_scn5a"], ["run_scn5a.py"]),
        "run_profile": (["run_scn5a"], ["run_profile.py"]),
        "run_bounds": (["run_scn5a", "run_profile"], ["run_bounds.py"]),
        "run_kv7_rank": ([], ["run_kv7_rank.py"]),
        "run_deconvolution": ([], ["run_deconvolution.py"]),
        "run_tp53": ([], ["run_tp53.py"]),
        "run_kv7_controls": ([], ["run_kv7_controls.py"]),
        "run_tp53_hotspot": (["run_tp53"], ["run_tp53_hotspot.py"]),
        "run_figure2": (["run_kv7_controls", "run_deconvolution", "run_tp53", "run_tp53_hotspot"], ["run_figure2.py"]),
        "run_kv7": ([], ["run_kv7.py"]),
        "run_clinical": (["run_scn5a"], ["run_clinical.py"]),
        "run_clinvar": (["run_scn5a"], ["run_clinvar.py"]),
        "run_survey": ([], ["run_survey.py"]),
        "run_survey_validation": (["run_survey"], ["run_survey_validation.py"]),
        "run_selection": ([], ["run_selection.py"]),
        "run_kv7_provenance": ([], ["run_kv7_provenance.py"]),
        "run_rank_null": ([], ["run_rank_null.py"]),
        "run_uncertainty": (["run_scn5a"], ["run_uncertainty.py"]),
        "run_schematic": ([], ["run_schematic.py"]),
        "run_titration": ([], ["run_titration.py"]),
        "run_titration_main": (["run_titration"], ["run_titration_main.py"]),
    }
    if args.skip_profile:
        del tasks["run_profile"]
        tasks["run_bounds"] = (["run_scn5a"], ["run_bounds.py"])
    pending, active = dict(tasks), {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        while pending or active:
            for name, (dependencies, command_args) in list(pending.items()):
                if len(active) == 4:
                    break
                if not all(dep in commands for dep in dependencies):
                    continue
                if any(commands[dep]["exit_code"] for dep in dependencies):
                    commands[name] = {"exit_code": -1, "error": "dependency failed"}
                    del pending[name]
                    continue
                active[executor.submit(execute, name, command_args)] = name
                del pending[name]
            if active:
                complete, _ = wait(active, return_when=FIRST_COMPLETED)
                for future in complete:
                    commands[active.pop(future)] = future.result()
            elif pending:
                raise RuntimeError(f"Unresolved dependencies: {list(pending)}")
    commands["pytest"] = execute("pytest", ["-m", "pytest", "-q", "-p", "no:cacheprovider", *TEST_FILES])

    comparisons = []
    prefixes = {str(SOURCE), str(SOURCE.resolve()), str(package), str(package.resolve())}
    archived_source = json.loads((baseline / "survey_summary.json").read_text()).get("source", "")
    source_suffix = "/data/archive/scoresets.json"
    if archived_source.endswith(source_suffix):
        prefixes.add(archived_source[:-len(source_suffix)])
    from matplotlib.image import imread
    for filename in all_outputs:
        if filename in skipped_outputs:
            continue
        original, regenerated = baseline / filename, package / filename
        kind = {".csv": "CSV", ".json": "JSON", ".png": "PNG pixels", ".pdf": "PDF bytes excluding dates"}[original.suffix]
        item = {"file": filename, "kind": kind, "match": False}
        try:
            if original.suffix == ".csv":
                a, b = pd.read_csv(original), pd.read_csv(regenerated)
                pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
                item.update(match=True, rows=len(b))
            elif original.suffix == ".json":
                a, b = [normalized(json.loads(p.read_text()), prefixes) for p in (original, regenerated)]
                diffs = differences(a, b)
                item.update(match=not diffs, differences=diffs[:30])
            elif original.suffix == ".png":
                a, b = imread(original), imread(regenerated)
                item.update(match=bool(a.shape == b.shape and np.array_equal(a, b)))
            else:
                a, b = original.read_bytes(), regenerated.read_bytes()
                clean = lambda value: re.sub(rb"/(?:CreationDate|ModDate)\s*\([^)]*\)", b"/Date <ignored>", value)
                item.update(match=bool(a.startswith(b"%PDF") and b.startswith(b"%PDF") and clean(a) == clean(b)))
        except (AssertionError, FileNotFoundError, ValueError) as exc:
            item["difference"] = str(exc)[:4000]
        comparisons.append(item)
    attempts = network_log.read_text().splitlines() if network_log.exists() else []
    (logs / "network_attempts.log").write_text("\n".join(attempts) + ("\n" if attempts else ""))
    state = {"archive": str(archive_path), "archive_sha256": archive_sha,
             "manifest_files_verified": len(manifest), "package": str(package), "python": str(python),
             "elapsed_seconds": round(time.time() - started, 2), "commands": commands,
             "comparisons": comparisons, "network_attempts": attempts,
             "external_inputs": external_inputs,
             "skipped": {"full_profile_simulation": args.skip_profile,
                         "retained_artifacts_not_reproduced": sorted(skipped_outputs)}}
    (report / "reproduction_results.json").write_text(json.dumps(state, indent=2) + "\n")
    matched = sum(item["match"] for item in comparisons)
    passed = all(command["exit_code"] == 0 for command in commands.values()) and matched == len(comparisons) and not attempts
    result_word = "Passed" if passed else "Differences or failures require review"
    profile_note = ("The unchanged full profile simulation was deliberately not rerun. Its packaged outputs were retained as inputs to the bounds cross-check and excluded from reproduced-artifact counts."
                    if args.skip_profile else "The profile driver ran in full mode; no --quick simulation was used.")
    text = f"""# Isolated reproduction verification — {date.today().isoformat()}

**{result_word}.** All requested commands ran from the actual S1 File archive in a fresh system temporary directory. Every manifest entry was checked before extraction. Excluded third-party inputs were supplied explicitly and verified against INPUT_RIGHTS.json before analysis. No undeclared source-tree fallback was used. Network connections were disabled during analysis and recorded if attempted. This verifies reproduction with the declared external inputs; it does not establish that the ZIP is a complete offline input bundle or that historical remote data remain retrievable.

{profile_note}

- Archive: `{archive_path}`
- Archive SHA-256: `{archive_sha}`
- Manifest files verified: {len(manifest)}
- Isolated analysis directory: `{package}`
- Interpreter: `{python}`
- Total elapsed time: {state['elapsed_seconds']:.2f} s
- Network attempts: {len(attempts)}
- Matching artifact comparisons: {matched}/{len(comparisons)}
- CSV/JSON numerical tolerance: relative and absolute 1e-12.
- JSON generation dates and absolute package prefixes were normalized; scientific metadata and numerical fields were retained. PNG pixels and PDF bytes were compared, excluding PDF creation/modification dates.

## Commands and exit status

| Driver | Exit status | Elapsed seconds | Log |
| --- | --- | --- | --- |
"""
    for name, command in commands.items():
        text += f"| {name} | {command['exit_code']} | {command.get('seconds', '—')} | [{name}](reproduction_logs/{name}.log) |\n"
    text += "\nExact commands are recorded in `reproduction_results.json`. Tests: `" + " ".join(TEST_FILES) + "`.\n\n## Artifact comparisons\n\n| Artifact | Comparison | Result |\n| --- | --- | --- |\n"
    for item in comparisons:
        text += f"| {item['file']} | {item['kind']} | {'match' if item['match'] else 'DIFFERENT'} |\n"
        if not item["match"]:
            print(json.dumps(item), flush=True)
    text += "\n## Scope and limitations\n\nThis verifies execution and artifact reproduction with the archived inputs and cached ClinVar records in the existing environment. It does not test a fresh dependency installation, validate live clinical records, or establish biological assumptions. Figure3 is the analysis filename for supplementary S1 Fig. TIFF exports and manuscript rendering were not regenerated in this check. Skipped or unrequested drivers and their retained outputs are not counted as reproduced.\n\nPytest result:\n\n```text\n" + (logs / "pytest.log").read_text().strip() + "\n```\n"
    (report / "reproduction_verification.md").write_text(text)
    print(json.dumps({"result": result_word, "package": str(package), "matched": f"{matched}/{len(comparisons)}",
                      "elapsed_seconds": state["elapsed_seconds"], "skipped": state["skipped"]}, indent=2), flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
