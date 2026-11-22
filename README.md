# Molecular dominance from paired variant assays

Code, data and machine-checked proofs for the paper *What paired variant assays can and cannot identify about molecular dominance* (Zhou Quan Eric Gao, University of California San Diego).

The paper asks what assays that measure a variant alone and mixed with wild type can identify about molecular dominance. It applies the answer to published data on the cardiac sodium channel (*SCN5A*), *KCNQ2*, *TP53* and ATP-sensitive potassium channels, and to the MaveDB archive. This repository holds everything needed to regenerate the reported numbers, tables and figures. It contains the paper's S1 File, minus that file's copy of the manuscript.

## Contents

| Path | What it holds |
|:--|:--|
| `analysis/` | Analysis scripts, redistributable inputs, derived outputs, figures and regression tests. [`analysis/README.md`](analysis/README.md) lists every command, its outputs and the order to run them in. |
| `analysis/lean/` | The Lean 4 formalization: 100 theorems in seven files, checked by the Lean kernel. [`analysis/lean/README.md`](analysis/lean/README.md) explains how to build it. |
| `analysis/INPUT_RIGHTS.json` | Licence, attribution, source and SHA-256 checksum for every input file. |
| `FILE_MANIFEST.json` | SHA-256 checksum of every file in S1 File. The eleven `manuscript/` entries refer to the manuscript copy, which is not in this repository. |
| `reproduction/` | Report of the verification run of 22 November 2026, described below. |

## Reproducing the analyses

The verified environment is Python 3.10.19. No GPU is needed.

```sh
cd analysis
python3 -m pip install -r requirements.txt
python3 restore_external_inputs.py --download      # or: --input-cache DIR (see below)
python3 -m pytest -q test_*.py
```

With every input restored, 146 tests pass and 2 skip. One skipped test checks the theorem table in the supplementary text, which is not in this repository. The other builds the Lean proofs and runs once Mathlib is downloaded. Without the restored inputs, 12 tests fail because the files they read are missing.

Then run the drivers in the order given in `analysis/README.md`, for example `python3 run_scn5a.py` and `python3 run_titration.py`.

To check the proofs, install Lean through [elan](https://github.com/leanprover/elan), then run:

```sh
cd analysis/lean
lake exe cache get
lake build
lake env lean Audit.lean
```

Every theorem should depend only on `propext`, `Classical.choice` and `Quot.sound`.

## Inputs that are not in this repository

Four source inputs keep licences that do not allow them to be redistributed here. `restore_external_inputs.py` accepts each one only if its SHA-256 checksum matches the value recorded in `analysis/INPUT_RIGHTS.json`.

- **O'Neill et al. (2022) supplementary workbook** ([PMC9262418](https://pmc.ncbi.nlm.nih.gov/articles/PMC9262418/)). No redistribution grant has been established for it. It is needed for the *SCN5A* analyses.
- **Zheng et al. (2022) supplementary dataset and supplementary methods** (*Genome Research*), licensed CC BY-NC 4.0. The dataset is needed to re-derive the *KCNQ4* extract.
- **18 historical MaveDB metadata records** from the snapshot of 26 September 2026, licensed CC BY-SA 4.0 or CC BY-NC-SA 4.0. The current MaveDB API no longer returns all of them in this form. They will be archived separately under their original licences, with a link added here. Once you have that archive, pass it with `--metadata-bundle FILE.zip`. The other 2,801 records of the snapshot are in `analysis/data/archive/scoresets_compatible.json`.

`restore_external_inputs.py --download` fetches what is still online. `--input-cache DIR` uses copies you already have, laid out under the paths listed in `INPUT_RIGHTS.json`. Restored files are listed in `.gitignore`, so they cannot be committed back by accident.

## Verification

On 22 November 2026, `analysis/submission-2026-09-29/reproduction_driver.py` reran all 21 analysis drivers except the long full-profile simulation, then the regression tests. It ran from the S1 File archive in a fresh temporary directory and made no network requests, after restoring the excluded inputs from checksum-matched local copies. All 63 compared outputs matched the packaged ones. CSV and JSON values matched to a relative and absolute tolerance of 10⁻¹², PNG files matched pixel for pixel, and PDF files matched byte for byte apart from their dates. `reproduction/` holds the report, its logs and the input manifest. Local file paths in it are replaced by placeholders such as `<temporary directory>`.

## Licence

The code is under the [MIT licence](LICENSE). Data files keep the licences of their sources, which `analysis/INPUT_RIGHTS.json` lists file by file. For example, the MaveDB score sets are CC0 and the Vanoye et al. (2022) workbooks are CC BY 4.0. The MIT licence does not cover third-party data.

## Citation

Please cite the paper once it is published. [`CITATION.cff`](CITATION.cff) gives the citation for this repository.

## AI use

The paper's "Use of AI tools" section describes how AI tools were used in this work.
