# scGSCA statistical tests

## Sample-phenotype test

The sample-phenotype test compares gene-set co-expression between two
phenotypes. Larger values of its nonnegative difference statistic indicate
greater differences. Its permutation P value includes ties and uses the
plus-one correction.

## Cell-type interaction test

Within each phenotype, D measures the co-expression contrast between the two
specified cell types. Group1 is the reference phenotype selected with
`--ref-letter`; Group2 is the other phenotype.

```text
Delta = D(Group2) - D(Group1)
p_abs = (1 + count(abs(Delta_perm) >= abs(Delta_observed))) / (B + 1)
```

B is the number of valid permutations. This absolute-delta test detects
differences in either direction and is the two-sided test used for the
manuscript. It is not calculated by doubling the signed right-tail P value.

Read `Permutation.Pvalue.AbsDelta` from `*_delta_result.tsv` for this result
(referred to as `abs.p.value` in the analysis workflow).

## Output fields

The main statistical columns in `*_delta_result.tsv` are listed below.
Group1 is the reference phenotype and Group2 is the other phenotype. Group
identifiers are stored as numeric character codes for the input phenotype suffixes
(for example, A = 65 and B = 66).

| Field | Meaning |
| --- | --- |
| `Group1` | Reference phenotype identifier |
| `Group2` | Non-reference phenotype identifier |
| `Observed.TS.Group1` | Observed co-expression contrast D between the two cell types within Group1 |
| `Observed.TS.Group2` | Observed co-expression contrast D between the two cell types within Group2 |
| `Observed.Delta` | D(Group2) minus D(Group1) |
| `Observed.AbsDelta` | Absolute value of the observed Delta |
| `Perm.Delta.Mean` | Mean signed Delta across permutations |
| `Perm.AbsDelta.Mean` | Mean absolute Delta across permutations |
| `Permutation.Pvalue.AbsDelta` | Two-sided absolute-delta permutation P value |
| `Permutation.Pvalue.Delta` | One-sided signed-delta permutation P value for D(Group2) > D(Group1), including ties |
| `Nperm.Paired` | Number of matched permutation draws used to calculate Delta from the two groups; this does not indicate a paired study design |
| `Delta.Definition` | Definition of the signed difference: D(Group2)-D(Group1) |

Use `Permutation.Pvalue.AbsDelta` for the two-sided cell-type interaction test.
Both P-value columns are computed in the same run.
Changing the reference reverses the observed Delta but not its magnitude.
Separate runs with reversed references can have different finite-permutation
draws, so their estimated absolute-delta P values need not match exactly.

## Example

```sh
./build/scgsca_run \
  --kegg-dir /path/to/pathway_inputs \
  --out-dir results/cell_type_interaction \
  --nperm 500 --seed 416627 --cores 4 \
  --ref-letter A --gene-start 3 --unpaired --skip-dc \
  --ct-pair "Cancer Epithelial,Myeloid"
```

Use your own cell-type labels and a new output directory. The cell-type
interaction procedure permutes subject-level phenotype labels and then
cell-type labels within subjects. Absolute-value testing does not itself
establish exchangeability or statistical calibration. See the README for
input-processing limitations and software checks.
