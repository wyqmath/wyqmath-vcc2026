# wyqmath-vcc2026

Reference transfer and single-cell count generation for the 2026 Virtual Cell
Challenge. Entry `cQLMHwkibt1kKx86ssNh` scored 0.163092 and placed 228th of 806.
[Official final leaderboard](https://virtualcellchallenge.org/leaderboard#final)

## Run

Use Python 3.11 or newer.

```sh
python -m pip install .
wyqmath-vcc smoke --output example.h5ad
```

The smoke test generates and reads back a synthetic H5AD with 48 cells and
6 genes.

```sh
wyqmath-vcc predict \
  --controls controls-final \
  --source-bank source_bank.npz \
  --output prediction.h5ad
```

See the [input specification](docs/inputs.md). Add `--patches final-overrides`
to apply the 15 [count overrides](docs/overrides.md) from the submitted recipe.

## Method

For each target, the model averages normalized response vectors from reference
studies. It samples integer counts from the target context's control templates
with a fixed total per gene. In context D, opposing source votes can attenuate
an effect, and 16 rounds of integer exchanges adjust mean cell-normalized
expression while preserving gene totals and cell depths.

A full prediction has 360,000 cells and 18,533 genes. Supply the controls,
prepared source effects, and any override blocks separately. Exact agreement
with the archived final matrix has not been verified. Use
[Arc's cell-eval2](https://github.com/ArcInstitute/cell-eval2) with matched
experimental reference data for evaluation.
