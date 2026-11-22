#!/usr/bin/env python3
"""Restore excluded third-party inputs for private reproduction, verifying hashes.

Run from the unpacked analysis directory:
  python3 restore_external_inputs.py --download
or use independently obtained, checksum-matching files in a local cache:
  python3 restore_external_inputs.py --input-cache /path/to/analysis-input-cache

The cache must mirror paths in INPUT_RIGHTS.json. A changed remote version is
never silently accepted; failed retrieval reports the primary access route.
The 18 historical MaveDB records that S1 File omits are archived at
https://doi.org/10.5281/zenodo.24438963; pass that ZIP with --metadata-bundle.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request
import zipfile


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def project(value, shape):
    if isinstance(shape, dict):
        if set(shape) == {"$constant"}:
            return shape["$constant"]
        return {key: project(value.get(key), child) for key, child in shape.items()}
    if isinstance(shape, list):
        if len(value) != len(shape):
            raise ValueError("Remote list length changed since the archived snapshot")
        return [project(item, child) for item, child in zip(value, shape)]
    return value


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "Dominance-reproduction/1.0 (source data retrieval)"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def read_historical_bundle(path, archive_manifest):
    """Read without extracting; pin every record to the S1 historical checksum."""
    if path is None:
        return {}
    expected = {row["urn"]: row for row in archive_manifest["omitted_records"]}
    with zipfile.ZipFile(path) as bundle:
        manifest = json.loads(bundle.read("HISTORICAL_METADATA_MANIFEST.json"))
        if manifest["full_snapshot_sha256"] != archive_manifest["full_sha256"]:
            raise ValueError("Historical bundle is for a different snapshot")
        rows = manifest["records"]
        if len(rows) != len(expected) or {r["urn"] for r in rows} != set(expected):
            raise ValueError("Historical bundle does not contain the expected records")
        records = {}
        for row in rows:
            payload = bundle.read(row["path"])
            reference = expected[row["urn"]]
            if sha(payload) != reference["sha256_json"]:
                raise ValueError(f"Historical bundle record checksum mismatch: {row['urn']}")
            value = json.loads(payload)
            if value["license"]["shortName"] != reference["license"]:
                raise ValueError("Historical source license differs")
            records[row["urn"]] = value
    return records


def restore(root, cache=None, download=False, metadata_bundle=None):
    root = Path(root).resolve()
    cache = Path(cache).resolve() if cache else None
    manifest = json.loads((root / "INPUT_RIGHTS.json").read_text())
    historical = read_historical_bundle(metadata_bundle, manifest["archive_metadata"])
    restored, failures = [], []
    for entry in manifest["files"]:
        if entry["included"] or entry["path"] == "data/archive/scoresets.json":
            continue
        target = root / entry["path"]
        if target.exists() and sha(target.read_bytes()) == entry["sha256"]:
            continue
        payload = None
        if cache and (cache / entry["path"]).exists():
            payload = (cache / entry["path"]).read_bytes()
        elif download:
            for url in entry.get("download_urls", []):
                try:
                    candidate = fetch(url)
                    if sha(candidate) == entry["sha256"]:
                        payload = candidate
                        break
                except Exception:
                    continue
        if payload is None or sha(payload) != entry["sha256"]:
            failures.append(f"{entry['path']}: obtain original at {entry['source_url']} (SHA-256 {entry['sha256']})")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        restored.append(entry["path"])
    archive = manifest["archive_metadata"]
    target = root / "data/archive/scoresets.json"
    if not target.exists() or sha(target.read_bytes()) != archive["full_sha256"]:
        cached = cache / "data/archive/scoresets.json" if cache else None
        if not historical and cached and cached.exists() and sha(cached.read_bytes()) == archive["full_sha256"]:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cached, target)
            restored.append("data/archive/scoresets.json")
        else:
            records = json.loads((root / "data/archive/scoresets_compatible.json").read_text())
            records.update(historical)
            for record in archive["omitted_records"]:
                if record["urn"] in historical:
                    continue
                try:
                    if not download:
                        raise ValueError("retrieval disabled")
                    value = project(json.loads(fetch(record["download_url"])), record["schema"])
                    if sha(json.dumps(value).encode()) != record["sha256_json"]:
                        raise ValueError("remote metadata differs from the archived snapshot")
                    records[record["urn"]] = value
                except Exception as error:
                    failures.append(f"{record['urn']}: {error}; {record['source_url']}")
            if all(urn in records for urn in archive["record_order"]):
                payload = json.dumps({urn: records[urn] for urn in archive["record_order"]}).encode()
                if sha(payload) != archive["full_sha256"]:
                    failures.append("Reconstructed archive checksum differs; obtain the 26 September 2026 snapshot")
                else:
                    target.write_bytes(payload)
                    restored.append("data/archive/scoresets.json")
    if failures:
        if any(item.startswith("urn:mavedb:") for item in failures):
            failures.append("The 18 historical records are archived at https://doi.org/10.5281/zenodo.24438963; "
                            "download MaveDB_2026-09-26_original_licenses.zip and pass it with --metadata-bundle.")
        raise RuntimeError("Excluded inputs not restored; no substitute accepted:\n" + "\n".join(failures))
    return {"restored": restored, "input_cache": str(cache) if cache else None,
            "historical_metadata_bundle": str(Path(metadata_bundle).resolve()) if metadata_bundle else None,
            "historical_metadata_bundle_sha256": sha(Path(metadata_bundle).read_bytes()) if metadata_bundle else None,
            "network_retrieval_enabled": download,
            "all_external_input_hashes_verified": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--input-cache", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--metadata-bundle", type=Path,
                        help="Separate historical MaveDB metadata ZIP under its original licenses")
    args = parser.parse_args()
    print(json.dumps(restore(args.root, args.input_cache, args.download, args.metadata_bundle), indent=2))


if __name__ == "__main__":
    main()
