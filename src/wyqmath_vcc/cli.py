import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .core import SOURCES, predict_context


CONTEXTS = {"D": "A", "E": "B", "F": "B"}
GENES = 18_533
TARGETS_PER_CONTEXT = 300
CELLS_PER_TARGET = 400
MAX_NNZ = 4_750_000_000
PATCH_TARGETS = {
    "D": ("ANP32B", "CARM1", "PHF6", "RBM5", "SMARCD1"),
    "E": ("PHF23", "PRMT6", "VIM"),
    "F": ("CUX1", "E2F1", "GATAD2B", "HMGB2", "PHF14", "PRDM15", "PTBP2"),
}


def read_table(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(required).issubset(reader.fieldnames or ()):
            raise ValueError(f"{path.name} requires columns {', '.join(required)}")
        return list(reader)


def unique_index(values, name):
    values = list(values)
    index = {value: i for i, value in enumerate(values)}
    if len(index) != len(values) or not all(values):
        raise ValueError(f"{name} must have unique, nonempty labels")
    return index


def read_panel(bundle):
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    expected = {
        "partition": "test",
        "panel_id": "vcc2026-test-1",
        "contexts": list(CONTEXTS),
        "n_genes": GENES,
        "cells_per_pert": CELLS_PER_TARGET,
    }
    if any(manifest.get(key) != value for key, value in expected.items()):
        raise ValueError("Use the official 2026 final controls bundle (D/E/F)")

    genes = [row["gene_name"] for row in read_table(bundle / "gene_names.csv", ("gene_name",))]
    gene_index = unique_index(genes, "Gene axis")
    if len(genes) != GENES:
        raise ValueError(f"Expected {GENES} genes")
    rows = read_table(bundle / "pert_counts.csv", ("target_gene",))
    has_context = bool(rows and "context" in rows[0])
    if has_context and {row["context"] for row in rows} != set(CONTEXTS):
        raise ValueError("Target table does not cover D/E/F")
    panels = {}
    for context in CONTEXTS:
        targets = [
            row["target_gene"] for row in rows
            if not has_context or row["context"] == context
        ]
        target_index = unique_index(targets, f"{context} targets")
        if len(targets) != TARGETS_PER_CONTEXT or not set(target_index) <= gene_index.keys():
            raise ValueError(f"Invalid target panel for {context}")
        panels[context] = targets
    return genes, panels


def load_controls(bundle, context, genes):
    controls = ad.read_h5ad(bundle / f"context_{context}.h5ad")
    if (controls.n_obs != 18_400 or controls.n_vars != GENES
            or "context" not in controls.obs or set(controls.obs["context"]) != {context}
            or "target_gene" not in controls.obs
            or set(controls.obs["target_gene"]) != {"non-targeting"}):
        raise ValueError(f"{context}: expected labeled non-targeting controls")
    gene_index = unique_index(controls.var_names, f"{context} control genes")
    if gene_index.keys() != set(genes) or controls.X is None:
        raise ValueError(f"{context}: control genes differ from the official gene axis")

    counts = sparse.csr_matrix(controls.X)
    if list(controls.var_names) != genes:
        counts = counts[:, [gene_index[gene] for gene in genes]]
    return counts


def load_bank(path, genes, panels):
    with np.load(path, allow_pickle=False) as bank:
        required = {"raw", "admitted", "present", "targets", "genes"}
        if not required <= set(bank.files):
            raise ValueError(f"Source bank requires keys {', '.join(sorted(required))}")
        if bank["genes"].astype(str).tolist() != genes:
            raise ValueError("Source bank gene order differs from the official axis")
        targets = bank["targets"]
        if targets.ndim != 1:
            raise ValueError("Source bank targets must be one-dimensional")
        index = unique_index(targets.astype(str), "Source bank targets")
        if "sources" in bank.files and bank["sources"].astype(str).tolist() != list(SOURCES):
            raise ValueError("Source bank uses an unexpected channel order")
        raw, admitted, present = bank["raw"], bank["admitted"], bank["present"]
    if (raw.shape != (len(SOURCES), len(targets), len(genes))
            or admitted.shape != raw.shape[:2] or present.shape != raw.shape[::2]):
        raise ValueError("Source bank arrays do not match the declared axes")
    missing = set().union(*panels.values()) - index.keys()
    if missing:
        raise ValueError(f"Source bank is missing {len(missing)} target labels: {', '.join(sorted(missing)[:5])}")
    return raw, admitted, present, index


def checked_block(matrix, cells, genes):
    if not sparse.isspmatrix_csr(matrix) or matrix.shape != (cells, genes):
        raise ValueError(f"Count block must be CSR with shape {(cells, genes)}")
    matrix.check_format(full_check=True)
    if (matrix.dtype != np.int32 or not matrix.has_canonical_format
            or (matrix.data <= 0).any()):
        raise ValueError("Count block must use sorted, unique, positive int32 stored values")
    depth = np.asarray(matrix.sum(axis=1, dtype=np.int64)).ravel()
    if (depth <= 0).any() or (depth >= 1_000_000).any():
        raise ValueError("Each generated row must have positive depth below 1,000,000")
    return matrix


def read_patch(path, cells, genes):

    with np.load(path, allow_pickle=False) as z:
        if "format" in z.files:
            matrix = sparse.load_npz(path)
        else:
            matrix = sparse.csr_matrix(
                (z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])
            )
    return checked_block(matrix, cells, genes)


def load_patches(directory, panels):
    if directory is None:
        return {}
    paths = {
        (context, target): directory / context / f"{target}.npz"
        for context, targets in PATCH_TARGETS.items()
        for target in targets
    }
    for (context, target), path in paths.items():
        if target not in panels[context]:
            raise ValueError(f"Published override {context}/{target} is not on this panel")
        if not path.is_file():
            raise ValueError(f"Expected all 15 published overrides. Missing {path}")
    return {key: read_patch(path, CELLS_PER_TARGET, GENES) for key, path in paths.items()}


def predicted_blocks(bundle, panels, genes, bank, patches):
    raw, admitted, present, index = bank
    for context, branch in CONTEXTS.items():
        indices = [index[target] for target in panels[context]]
        blocks = predict_context(
            load_controls(bundle, context, genes), raw[:, indices], admitted[:, indices],
            present, branch, cells_per_target=CELLS_PER_TARGET,
        )
        for target, block in zip(panels[context], blocks, strict=True):

            yield patches.get((context, target), block)


def prediction_annotations(panels, genes, cells):
    rows = [(context, target) for context, targets in panels.items() for target in targets]
    obs = pd.DataFrame(
        {
            "context": pd.Categorical(np.repeat([context for context, _ in rows], cells),
                                      categories=list(panels)),
            "target_gene": pd.Categorical(np.repeat([target for _, target in rows], cells)),
        },
        index=[
            f"{context}_{i:03d}_{j:03d}"
            for context, targets in panels.items()
            for i in range(len(targets))
            for j in range(cells)
        ],
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene_name"))
    return obs, var


def write_prediction(output, panels, genes, cells, blocks):
    if output.exists():
        raise FileExistsError(output)
    obs, var = prediction_annotations(panels, genes, cells)
    shape = (len(obs), len(var))
    with tempfile.NamedTemporaryFile(dir=output.parent, suffix=".h5ad", delete=False) as temp:
        staged = Path(temp.name)
    try:
        ad.AnnData(sparse.csr_matrix(shape, dtype=np.int32), obs=obs, var=var).write_h5ad(staged)
        with h5py.File(staged, "r+") as handle:
            del handle["X"]
            group = handle.create_group("X")
            group.attrs.update({
                "encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": shape,
            })
            chunk = (min(1_048_576, shape[0] * shape[1]),)
            data = group.create_dataset(
                "data", (0,), maxshape=(None,), dtype="i4", chunks=chunk, compression="lzf"
            )
            indices = group.create_dataset(
                "indices", (0,), maxshape=(None,), dtype="i4", chunks=chunk, compression="lzf"
            )
            pointers = group.create_dataset("indptr", (shape[0] + 1,), dtype="i8")
            pointers[0] = 0
            offset = 0
            for i, (_, block) in enumerate(zip(range(shape[0] // cells), blocks, strict=True)):
                checked_block(block, cells, shape[1])
                end = offset + block.nnz
                if end > MAX_NNZ:
                    raise ValueError("Prediction exceeds the final panel sparse-entry cap")
                data.resize((end,))
                indices.resize((end,))
                data[offset:end] = block.data
                indices[offset:end] = block.indices
                pointers[i * cells + 1:(i + 1) * cells + 1] = block.indptr[1:].astype(np.int64) + offset
                offset = end

        os.link(staged, output)
    finally:
        staged.unlink()


def smoke(output):
    rng = np.random.default_rng(9)
    genes = [f"g{i}" for i in range(6)]
    panels = {"D": genes[:2], "E": genes[:2], "F": genes[:2]}
    cells = 8
    controls = sparse.csr_matrix(rng.poisson(30, size=(16, len(genes))) + 1)
    raw = np.ones((len(SOURCES), 2, len(genes)), dtype=np.float32)
    raw[1::2] *= -1


    raw[:6, 0, 0] = -0.02
    raw[6:, 0, 0] = 0.10
    raw[:6, 1, 1] = 0.02
    raw[6:, 1, 1] = -0.10
    admitted = np.ones(raw.shape[:2], dtype=bool)
    present = np.ones((len(SOURCES), len(genes)), dtype=bool)
    blocks = (
        block
        for branch in CONTEXTS.values()
        for block in predict_context(controls, raw, admitted, present, branch, cells_per_target=cells)
    )
    write_prediction(output, panels, genes, cells, blocks)
    result = ad.read_h5ad(output)
    if (result.shape != (48, 6) or list(result.var_names) != genes
            or not sparse.isspmatrix_csr(result.X) or result.X.dtype != np.int32):
        raise ValueError("Synthetic H5AD readback failed")
    print(f"Synthetic smoke check passed: {result.n_obs} cells, {result.n_vars} genes ({output})")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate single-cell count predictions from reference effects.")
    commands = parser.add_subparsers(dest="command", required=True)
    synthetic = commands.add_parser("smoke", help="run the predictor and H5AD writer on synthetic counts")
    synthetic.add_argument("--output", type=Path, required=True)
    predict = commands.add_parser("predict", help="generate predictions on the official final panel")
    predict.add_argument("--controls", type=Path, required=True)
    predict.add_argument("--source-bank", type=Path, required=True)
    predict.add_argument("--patches", type=Path, help="directory containing all 15 final override blocks")
    predict.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "smoke":
            smoke(args.output)
            return
        if args.output.exists():
            raise FileExistsError(args.output)
        genes, panels = read_panel(args.controls)
        patches = load_patches(args.patches, panels)
        bank = load_bank(args.source_bank, genes, panels)
        blocks = predicted_blocks(args.controls, panels, genes, bank, patches)
        write_prediction(args.output, panels, genes, CELLS_PER_TARGET, blocks)
        print(f"Saved {args.output} with {len(patches)} archived override blocks")
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
