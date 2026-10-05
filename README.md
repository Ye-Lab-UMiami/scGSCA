# scGSCA

Gene-set differential co-expression analysis for single-cell RNA sequencing data.

## Creator and contributors

- Creator: Zhuoli Jin
- Contributors: Zhuoli Jin, Rebecca Irlmeier, and Fei Ye
- Repository: https://github.com/Ye-Lab-UMiami/scGSCA

## Associated manuscript

**scGSCA: gene-set differential co-expression analysis for single-cell RNA sequencing data**

Zhuoli Jin, Rebecca Irlmeier, Jie Ping, Yuan Wang, Yan Guo, and Fei Ye.

Zhuoli Jin and Rebecca Irlmeier
contributed equally to this work. Fei Ye is the contact corresponding author.
Publication details and the article DOI will be added when available.

## Overview

scGSCA analyzes differences in gene-set co-expression in single-cell RNA
sequencing data. It fits gene-wise zero-inflated negative binomial mixed
models with subject-level random intercepts and constructs gene-set
co-expression statistics from quantile-residual correlations.

Two statistical tests are available:

- **Sample-phenotype test:** compares overall gene-set co-expression between
  two sample phenotypes while including cell type as a fixed effect.
- **Cell-type interaction test:** assesses whether co-expression contrasts
  between two cell types differ across sample phenotypes. Use the two-sided
  absolute-delta permutation P value, `Permutation.Pvalue.AbsDelta`, for this test.

The program produces pathway-level statistics and permutation P values.
See the input-processing and statistical-behavior notes below before use.

## Implementation

Source distribution of scGSCA for single-cell gene-set
co-expression analysis. This directory contains the C++ command-line
implementation, not an installable R package.

## Requirements

- C++17 compiler
- CMake 3.16 or later
- Eigen3
- NLopt (including development headers)
- OpenMP (optional, for parallel execution)
- Python 3 (only for the optional integration test)

For macOS with Homebrew:

```sh
brew install cmake eigen nlopt libomp
```

For Debian or Ubuntu:

```sh
sudo apt-get install build-essential cmake libeigen3-dev libnlopt-dev
```

## Build and test

Run from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
./build/scgsca_run --help
```

For dependencies installed outside standard locations, pass their installation
prefixes through `-DCMAKE_PREFIX_PATH` when configuring CMake.

## Input

Supply one tab-separated `.tsv` file per pathway. Rows represent cells. The
gastric case-study inputs inspected for this release contain **relative-count
(RC) normalized expression**, rather than raw UMI counts. The located RC
preprocessing scripts divide each gene's count by the cell's total count
and multiply by 10,000, without log transformation. Compute cell totals from
the full expression matrix, not separately within each pathway.

The executable does not perform RC normalization. It rounds input expression
values to integers before ZINB model fitting and quantile-residual calculation.
This describes the implemented workflow, not a guarantee that rounded RC
values follow a ZINB distribution. Rounding does not recover raw counts;
model adequacy and permutation calibration require separate validation.
Do not substitute log-normalized expression or standardized (z-score) values.

With the default `--gene-start 3`, the first
three columns are metadata and gene columns begin at the fourth column:

```text
Cell_ID    Celltype    Samples    GENE1    GENE2    ...
```

The header above illustrates the column layout; actual files must use tabs.
Use `Cell_ID`, `Celltype`, and `Samples` as metadata column names. The current
C++ CLI also recognizes these aliases automatically:

| Canonical column | Supported aliases |
| --- | --- |
| `Cell_ID` | `cell_id`, `CellID`, `cellid` |
| `Celltype` | `cell_type`, `celltype` |
| `Samples` | `samples`, `sample` |

There are currently no CLI arguments for supplying arbitrary metadata column
names. Rename other column names before running; `--gene-start` only selects
the zero-based index of the first gene column, not metadata column names. Cell-type
names supplied to `--ct-pair` must match the input values exactly.

The final character of each `Samples` value encodes its phenotype group
(for example, `sample01_A` and `sample02_B`). Exactly two suffix levels are
required. `--ref-letter A` selects A as the reference group. The current
implementation uses `Samples` as its subject grouping variable; it does not
read a separate donor ID column. Verify that this encoding and the permutation
design are appropriate for your study, particularly for repeated-donor data.

The historical option name `--kegg-dir` is also used for pathway files from
other collections, including Hallmark. No participant-level data or pathway
database files are bundled.

## Run

Sample-phenotype test only (unpaired example):

```sh
./build/scgsca_run \
  --kegg-dir /path/to/pathway_inputs \
  --out-dir results/sample_phenotype \
  --nperm 500 --seed 416627 --cores 4 \
  --gene-start 3 --ref-letter A --unpaired --skip-ct
```

Cell-type interaction test only (replace cell-type names with your labels):

```sh
./build/scgsca_run \
  --kegg-dir /path/to/pathway_inputs \
  --out-dir results/cell_type_interaction \
  --nperm 500 --seed 416627 --cores 4 \
  --gene-start 3 --ref-letter A --unpaired --skip-dc \
  --ct-pair "Epithelial cell,Immune"
```

Repeat `--ct-pair` to request additional contrasts. Specify the desired pairs
explicitly, using labels present in your own data. Both cell types must be
present in both phenotype groups for a comparison to run. The current executable
still attempts `Astrocyte,Endothelial_cells` and `Astrocyte,Macrophage` if no
pair is supplied; it does not automatically select the first two cell types.
These defaults are not general recommendations.

Use a new output directory for each analysis. `--paired` selects the paired
permutation procedure (implemented by subject-level phenotype-label flips);
it does not infer donor matching from a separate metadata table. `--unpaired`
selects the unpaired subject-label permutation procedure.

## Statistical behavior

The cell-type interaction test used for the manuscript is the two-sided
absolute-delta permutation test. With
`Delta = D(non-reference) - D(reference)`, it counts permutations satisfying
`abs(Delta_perm) >= abs(Delta_observed)` and applies the plus-one correction.
Read `Permutation.Pvalue.AbsDelta` in `*_delta_result.tsv`; this is the
absolute-delta P value referred to as `abs.p.value` in the analysis workflow.

See [Output fields](README_TESTS.md#output-fields) for the result columns and
their meanings, including the two-sided absolute-delta P value and the
additional one-sided signed-delta P value.

The CT procedure permutes subject-level phenotype labels and then cell-type
labels within subjects. Packaging and passing software tests do not establish
exchangeability or statistical calibration for a particular experimental design.

## Software tests

These checks test the implementation; they are distinct from the
sample-phenotype test and cell-type interaction test used for statistical inference.

CTest runs the fixed-right-tail unit test. An additional integration test
requires a suitable input directory with groups A/B and cell-type labels
`Epithelial cell` and `Immune`:

```sh
python3 tests/smoke_ct_tail.py /path/to/test_inputs
```

This integration test checks the CLI, CT output metadata, and reference reversal.
It is not a statistical calibration test.

## Distribution contents

`src/` and `include/` contain the implementation; `tests/` contains the existing
tests. Generated plots, compiled binaries, bundled dynamic libraries, build
caches, local server scripts, and analysis results are excluded.

## License

This project's source code is distributed under the [MIT License](LICENSE).
Third-party dependencies retain their respective licenses.
