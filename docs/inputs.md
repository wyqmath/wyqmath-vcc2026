# Prediction inputs

## Controls

`--controls` takes the extracted VCC 2026 final bundle. Obtain it through the
[VCC data page](https://virtualcellchallenge.org/datasets) or
[official CLI](https://vcc-cli-wiki.virtualcellchallenge.org/).

The bundle must contain `manifest.json`, `gene_names.csv`, `pert_counts.csv`,
and `context_D.h5ad`, `context_E.h5ad`, and `context_F.h5ad`. Keep the control
rows in their original order and use the gene and target orders in the CSVs.

## Source bank

Pass a NumPy `.npz` file to `--source-bank`. Use numeric, Boolean, or Unicode
arrays.

| Array | Shape and type | Meaning |
| --- | --- | --- |
| `raw` | `(8,T,G)`, float32 | Projected source effects |
| `admitted` | `(8,T)`, bool | Source has an eligible row for the target |
| `present` | `(8,G)`, bool | Gene maps to the source axis |
| `targets` | `(T,)`, strings | Unique labels covering the final panel |
| `genes` | `(G,)`, strings | Exact order from `gene_names.csv` |
| `sources` | optional `(8,)`, strings | Channel names in the order below |

Channel order is `gw`, `hct116`, `hek293t`, `cd4_Rest`, `cd4_Stim8hr`,
`cd4_Stim48hr`, `dld1v2`, `h1scmg`. The submission used
[Replogle K562](https://doi.org/10.1016/j.cell.2022.05.013), the HCT116 and
HEK293T Orion screens, three CD4 conditions, a pooled DLD1 screen, and
[H1 SCMG](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295214).

Each target needs at least one admitted source. Admitted rows must be finite
and contain zeros at unmapped genes. An admitted row of zeros counts toward the
average. The code starts from projected effect arrays. The original raw data
processing is outside the package.

## Optional count overrides

`--patches` takes a directory with one subdirectory per context and a
`<target>.npz` file for each of the [15 groups](overrides.md). Supply all 15
blocks in the official gene order and original cell order. Each block must
be canonical CSR with int32 counts and shape `(400, 18533)`. The loader accepts
SciPy CSR archives and arrays named `data`, `indices`, `indptr`, and `shape`.
