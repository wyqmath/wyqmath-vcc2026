# wyqmath-vcc2026

A reference-transfer method for predicting single-cell RNA count distributions
after gene perturbation in the 2026 Virtual Cell Challenge.

| Method | Prediction recipe | Final test result |
| --- | --- | --- |
| Backbone | Same-target effects from eight reference channels, transferred into the official controls | Base generator for Final15 |
| Final15 | Backbone with 15 context-target blocks replaced using matched public perturbation studies | **Overall 0.1631, rank 228 of 806**, our best official result |
| Experimental18 | Final15 plus three independently selected replacements in context E | Lower error on held-out source data for all three targets. No official score. |

[Official final leaderboard](https://virtualcellchallenge.org/leaderboard#final)

## Reference transfer and count generation

The final panel comprises three anonymized backgrounds (D, E, and F), 300
perturbations per background, and 18,533 genes. The input for each background
is its own set of non-targeting control cells. The source cell lines in the
tables identify the public reference experiments.

1. **Combine same-target effects.** Project each eligible reference vector
   onto the final gene axis, with zeros for unmapped genes. Divide each
   nonzero vector by its L2 norm and average across eligible channels. An
   eligible all-zero row enters the denominator. An unavailable target does
   not enter the average. The sources use different effect estimators, so
   vector normalization limits the influence of their original units.
2. **Set the effect magnitude by background.** Multiply the consensus by
   0.85 times a scale of about 67.77 for D or 24.61 for E and F.
   D also uses a dose factor of about 0.838. In D, halve a gene's
   effect if at least three sources provide nonzero votes and at least
   two-thirds of those votes oppose the consensus.
3. **Generate 400 count profiles per target.** Draw 400 distinct control
   cells in each background with PCG64 seed 1, sort by their input row
   index, and shrink each template cell to 45% of its original profile
   plus 55% of the template mean. For gene `g` with mean control count
   `lambda[g]` and transferred effect `d[g]`, round
   `400 * lambda[g] * 2**d[g]` to the nearest integer. Allocate that total
   across cells by multinomial sampling in proportion to their template
   counts for the gene. PCG64 seed 0 starts a continuous sampling stream
   for each background in official target order. A gene-level replacement
   uses a copy of the original draw state, preserving the stream for the
   rest of the prediction.
4. **Calibrate D and apply study-specific replacements.** In D, 16 rounds
   of integer exchanges between cells adjust the mean of the per-cell
   expression proportions. The target starts from the mean per-cell
   proportions in the controls, reweights them by `2**d[g]`, and renormalizes
   over genes with positive simulated totals.
   Each exchange preserves each cell's library size and each gene's total
   count. E and F use the sampled counts directly. Final15 then replaces
   only the blocks listed below. Experimental18 extends these replacements
   to three more E targets.

The predictor writes a sparse H5AD with 360,000 cells and 18,533 genes.
`X` contains raw int32 counts, and `obs` identifies the context and target.
The method uses estimated source effects and fixed transfer settings. It
has no neural-network training step.

## Reference datasets

The eight channels are listed in bank order. Their historical effect
estimators are part of the reference input.

| Channel | Dataset | Effect estimator |
| --- | --- | --- |
| `gw` | Replogle et al., K562 genome-wide CRISPRi. [Paper](https://doi.org/10.1016/j.cell.2022.05.013), [data](https://plus.figshare.com/ndownloader/files/35775507) | Pooled counts rescaled to a common depth, then a log2 ratio to control with pseudocount 1 |
| `hct116` | HCT116, [X-Atlas-Orion](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) | Log2 ratio of per-cell means to pooled non-targeting controls in the same background, pseudocount 1 |
| `hek293t` | HEK293T, [X-Atlas-Orion](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) | Same estimator with the HEK293T control pool |
| `cd4_Rest` | Resting CD4 T cells, [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) | Author-estimated differential effects from `GWCD4i.DE_stats.h5ad`, layer `log_fc` |
| `cd4_Stim8hr` | CD4 T cells after 8 h of stimulation, [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) | Author effects for the 8 h condition |
| `cd4_Stim48hr` | CD4 T cells after 48 h of stimulation, [GSE314342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) | Author effects for the 48 h condition |
| `dld1v2` | Pooled DLD-1 screen, [GSE337988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337988) | Unique highest-UMI guide assignment across 44 expression/guide units, then log2 ratio of per-cell means to pooled controls, pseudocount 1 |
| `h1scmg` | H1 embryonic stem-cell CRISPRi screen, [GSE295214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295214) | Log2 ratio of guide-assigned per-cell means to pooled controls, pseudocount 1 |

The CD4 input takes the first row per target in each condition after
excluding rows marked `low_target_gex`. See the
[author data guide](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/data_sharing_readme.md).
The H1 channel comes from the SCMG study, separate from VCC 2025 training
data.

## Final15 replacements

Each replacement covers one context-target group and contains 400 cells
on the full gene axis.

| Context / target | Source and effect input | Evidence or scope |
| --- | --- | --- |
| D / ANP32B | [GSE242849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242849), MV4-11, relative FPKM | Deposited NC/S3 contrast. Viral-control equivalence is unresolved |
| D / CARM1 | [GSE103528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103528), MOLM13, RSEM expected reads | Day 3 contrast with both hairpins |
| D / PHF6 | [GSE277149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277149), THP-1, featureCounts | Clone/guide assignments for individual repeats are unresolved |
| D / RBM5 | [GSE225633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225633), MOLM13, RSEM estimates | Three replicate-label holdouts |
| D / SMARCD1 | [GSE200201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200201), MOLM13, single-cell gene-expression counts | The final effect uses channels 1 to 10, with channels 11 to 15 as comparison data |
| E / PHF23 | [GSE300095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300095), basal HeLa, featureCounts | KO/WT contrast. RNA clone independence and WT editing equivalence are unresolved |
| E / PRMT6 | [GSE221971](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221971), U87, deposited gene counts | Three library holdouts using deposited counts |
| E / VIM | [GSE286944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286944), U251-MG, featureCounts | Nine correlated holdout views |
| F / CUX1 | [GSE129487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129487), four rheumatoid arthritis synovial fibroblast donors, kallisto estimated counts | Refines supported coordinates of the earlier [GSE154148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154148) author-effect block |
| F / E2F1 | [GSE77260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77260), Wharton's jelly mesenchymal cells, Partek estimates | Six whole-donor holdouts |
| F / GATAD2B | [GSE154148-associated workbook](https://doi.org/10.1038/s41556-021-00725-7), author-corrected differential effect | Consistency was checked against the supplied author effect |
| F / HMGB2 | [GSE269856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269856), proliferating IMR90, gene counts | Library holdouts. WT transfection identity is unresolved |
| F / PHF14 | [GSE131576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131576), proliferating WI38, HTSeq counts | Paired library holdouts |
| F / PRDM15 | [GSE154148-associated workbook](https://doi.org/10.1038/s41556-021-00725-7), author-corrected differential effect | Consistency was checked against the supplied author effect |
| F / PTBP2 | [GSE210131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210131), AG22056 conversion day 3, featureCounts | AMnp/AMp conditional comparison. Viral-burden balance is unresolved |

## Experimental18

The candidate retains all 15 Final15 replacements and adds the following
three E-context blocks. Each comparison fits a source effect without the
held-out libraries. The source effect is a log2 ratio of pooled counts
normalized to 10,000 across the common mapped genes, with a pseudocount
of 0.1. The effect replaces the backbone value only at measured genes
with an absolute source effect of at least 0.5. Other gene columns retain
their original counts. The final local blocks use all listed libraries
after holdout evaluation.

| E target and source | Holdout design | Response MSE | All mapped non-target MSE |
| --- | --- | ---: | ---: |
| NR2C2, [GSE299760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299760). MCF7 and T47D, 72 h siRNA vs. scrambled, 3 + 3 STAR count libraries per background | Fit on either background and hold out the other | 0.280 to 0.140 | 0.043 to 0.038 |
| ADNP, [GSE79395](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79395). HCT116, siRNA vs. control, 3 + 3 HTSeq count libraries | Nine crossed treatment-control holdouts, with 2 + 2 libraries for each fit | 0.378 to 0.112 | 0.028 to 0.025 |
| ZBTB48, [GSE228606](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228606). HEK293, 72 h siRNA vs. non-targeting control, 2 + 2 RSEM expected-count libraries | Four crossed treatment-control holdouts, with 1 + 1 libraries for each fit | 0.574 to 0.183 | 0.077 to 0.049 |

The values are mean squared errors before and after the replacement,
averaged across holdouts. They compare log2 fold changes calculated from
the generated counts with held-out responses. Response MSE uses the
training-defined response genes. Both metrics exclude the perturbed gene.
The ADNP and ZBTB48 holdouts reuse libraries across comparisons. These
experiments support a candidate for future evaluation on independent
perturbation measurements from the target backgrounds.

To build the three input tables, use GSM9045715-17 and GSM9045721-23 as
the NR2C2 treatment libraries and GSM9045712-14 and GSM9045718-20 as
controls. Keep MCF7 and T47D separate for the holdout comparison. Use
GSM2094364-66 vs. GSM2094367-69 for ADNP, and GSM7133917-18 vs. GSM7133919-20
for ZBTB48, with treatment first. Match the sample columns to the deposited
tables, and use only genes with a unique source-to-panel mapping.

## Run

Use Python 3.11 or newer.

```sh
python -m pip install .
wyqmath-vcc smoke --output example.h5ad
```

The smoke test generates and reads back a synthetic H5AD with 48 cells
and six genes.

### Prepare inputs

Obtain the extracted final control bundle through the
[VCC data page](https://virtualcellchallenge.org/datasets) or the
[official CLI](https://vcc-cli-wiki.virtualcellchallenge.org/). The bundle
contains `manifest.json`, `gene_names.csv`, `pert_counts.csv`, and
`context_D.h5ad`, `context_E.h5ad`, and `context_F.h5ad`. Preserve the
original control row order.

To estimate a new effect from processed expression data, provide a
tab-separated count table with a `gene` column and one column per sample,
and a tab-separated sample table with `sample`, `target`, and
`condition` columns. Use `perturbed` or `control` as the condition and
`*` as the target for shared controls. Use unique gene symbols on the
mapped analysis axis before normalization. Fractional expected counts,
such as RSEM estimates, are supported. The default pools samples by
condition before normalization. To average normalized libraries instead,
use `--normalization library-mean-cp10k`.

```sh
wyqmath-vcc estimate-source \
  --counts source-counts.tsv \
  --samples source-samples.tsv \
  --output source-effects.npz
```

For the reference bank, supply one `<channel>.npz` file for each of the
eight channel names above. Each file holds Unicode `genes` and `targets`
vectors and a float32 or float64 `effects` array of shape `(targets, genes)`.
Include only eligible targets and measured genes. Use the historical estimator
listed for each channel when reconstructing Final15. The generic count
estimator is suitable for the pooled-count recipe in Experimental18.

```sh
wyqmath-vcc prepare-bank \
  --controls controls-final \
  --sources reference-effects \
  --output source-bank.npz
```

The bank stores `raw` (projected float32 effects, channel by target by gene),
`admitted` (channel by target), `present` (channel by gene), and the
`targets`, `genes`, and `sources` axes. Each final target needs at least
one eligible reference row.

### Generate predictions

```sh
wyqmath-vcc predict \
  --controls controls-final \
  --source-bank source-bank.npz \
  --variant final15 \
  --patches final-overrides \
  --output final15.h5ad
```

A count override is a file named `<context>/<target>.npz` under the override
directory. Use canonical CSR int32 counts in official gene and sampled cell
order and shape `(400, 18533)`. The loader accepts SciPy sparse NPZ files
or the arrays `data`, `indices`, `indptr`, and `shape`. To generate the
backbone, use `--variant backbone` and omit `--patches`.

```sh
wyqmath-vcc predict \
  --controls controls-final \
  --source-bank source-bank.npz \
  --variant experimental18 \
  --patches final-overrides \
  --experimental-effects experimental-effects \
  --output experimental18.h5ad
```

The experimental effects directory contains `NR2C2.npz`, `ADNP.npz`,
and `ZBTB48.npz`. For each file, use `estimate-source` with the default
pooled normalization and pseudocount on the common, uniquely mapped genes
from that study. Each file contains one target. If all 18 count overrides
are available, point `--patches` to their directory and omit
`--experimental-effects`.

The repository includes count-table effect estimation, bank construction,
and prediction generation. Obtain the controls, study inputs, historical
effect arrays, and Final15 count overrides separately. Rebuilding the
original source-specific inputs from raw sequencing data requires
preprocessing outside this repository. To evaluate predictions against
matched experimental data, use
[Arc's cell-eval2](https://github.com/ArcInstitute/cell-eval2).
