# Reproduction report, 22 November 2026

`analysis/submission-2026-09-29/reproduction_driver.py` wrote these files when it verified the S1 File archive with SHA-256 `6da63d1597fb28705cb95c4d3dc2cb975dbd8a257fa654b44bdcc47c9ff74198`. That archive holds the same `analysis/` files as this repository. The 18 historical MaveDB records came from `MaveDB_2026-09-26_original_licenses.zip` (SHA-256 `5c83450ba7744c59dc5c567dabe67661d4e9444f33156feef33c4be895f79e74`), the file archived at <https://doi.org/10.5281/zenodo.24438963>.

- `reproduction_verification.md`: the summary. It records 246 checksum-verified files, 21 drivers and the regression tests run, 63 of 63 matching outputs, and no network requests.
- `reproduction_results.json`: the same results in machine-readable form.
- `reproduction_input_manifest.json`: the excluded inputs that were supplied, with their checksums.
- `reproduction_logs/`: the output of each command.

Absolute paths on the verifying computer were replaced by placeholders such as `<temporary directory>`. Nothing else was changed.
