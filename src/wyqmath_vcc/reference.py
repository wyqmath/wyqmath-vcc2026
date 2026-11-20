from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

import numpy as np

from .core import SOURCES

__all__ = ["prepare_bank", "estimate_source"]


def _labels(values, name, *, duplicates=False, empty=False):
    values = np.asarray(values)
    if values.ndim != 1 or values.dtype.kind != "U":
        raise ValueError(f"{name} must be a one-dimensional Unicode array")
    labels = values.tolist()
    if not empty and not labels:
        raise ValueError(f"{name} must contain at least one label")
    seen = set()
    for label in labels:
        if not label.strip():
            raise ValueError(f"{name} contains an empty label")
        if not duplicates and label in seen:
            raise ValueError(f"{name} contains duplicate label {label!r}")
        seen.add(label)
    return labels


def _check_output(output):
    if os.path.lexists(output):
        raise FileExistsError(f"Output already exists: {output.name}")


def _save(output, **arrays):
    _check_output(output)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=output.parent, prefix=f".{output.name}.",
            suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink()


def prepare_bank(
    source_dir: Path, genes: list[str], targets: list[str], output: Path,
) -> None:
    source_dir, output = Path(source_dir), Path(output)
    _check_output(output)
    genes = _labels(genes, "Output genes")
    targets = _labels(targets, "Output targets")
    raw = np.zeros((len(SOURCES), len(targets), len(genes)), dtype=np.float32)
    admitted = np.zeros(raw.shape[:2], dtype=bool)
    present = np.zeros((len(SOURCES), len(genes)), dtype=bool)
    for channel, source in enumerate(SOURCES):
        path = source_dir / f"{source}.npz"
        if not path.is_file():
            raise FileNotFoundError(f"Missing source file: {path.name}")
        with np.load(path, allow_pickle=False) as data:
            if not {"genes", "targets", "effects"} <= set(data.files):
                raise ValueError(f"{source}: requires arrays genes, targets, effects")
            native_genes = _labels(
                data["genes"], f"{source} genes", duplicates=source == "dld1v2",
            )
            native_targets = _labels(data["targets"], f"{source} targets", empty=True)
            effects = data["effects"]
        if (effects.dtype.kind != "f" or effects.dtype.itemsize not in (4, 8)
                or effects.shape != (len(native_targets), len(native_genes))):
            raise ValueError(
                f"{source}: effects must be float32 or float64 with shape (targets, genes)"
            )
        with np.errstate(over="ignore", invalid="ignore"):
            effects = effects.astype(np.float32, copy=False)
        if not np.isfinite(effects).all():
            raise ValueError(f"{source}: effects must be finite after float32 conversion")
        columns = {}
        for column, gene in enumerate(native_genes):
            columns.setdefault(gene, []).append(column)
        if source == "dld1v2":
            with np.errstate(over="ignore", invalid="ignore"):
                for indices in columns.values():
                    if len(indices) > 1:
                        effects[:, indices[0]] = effects[:, indices].mean(
                            axis=1, dtype=np.float32,
                        )
            if not np.isfinite(effects).all():
                raise ValueError("dld1v2: duplicate gene means must be finite")
        rows = {target: row for row, target in enumerate(native_targets)}
        output_rows = [row for row, target in enumerate(targets) if target in rows]
        input_rows = [rows[targets[row]] for row in output_rows]
        output_columns = [column for column, gene in enumerate(genes) if gene in columns]
        input_columns = [columns[genes[column]][0] for column in output_columns]
        admitted[channel, output_rows] = True
        present[channel, output_columns] = True
        raw[channel][np.ix_(output_rows, output_columns)] = effects[
            np.ix_(input_rows, input_columns)
        ]
        del effects
    missing = [target for target, covered in zip(targets, admitted.any(axis=0)) if not covered]
    if missing:
        raise ValueError(f"No admitted source for targets: {', '.join(missing)}")
    _save(
        output, raw=raw, admitted=admitted, present=present,
        genes=np.asarray(genes, dtype=str), targets=np.asarray(targets, dtype=str),
        sources=np.asarray(SOURCES, dtype=str),
    )


def _read_samples(path):
    samples = {}
    targets = []
    seen_targets = set()
    arms = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if (not {"sample", "target", "condition"} <= set(header)
                or len(header) != len(set(header))):
            raise ValueError("Samples TSV requires unique columns sample, target, condition")
        positions = [header.index(name) for name in ("sample", "target", "condition")]
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Samples TSV row {line} has an incorrect number of fields")
            sample, target, condition = [row[position] for position in positions]
            if not sample.strip() or not target.strip():
                raise ValueError(f"Samples TSV row {line} requires nonempty sample and target")
            if sample in samples:
                raise ValueError(f"Samples TSV contains duplicate sample {sample!r}")
            if condition not in ("perturbed", "control"):
                raise ValueError(f"Sample {sample!r}: condition must be perturbed or control")
            if target == "*" and condition != "control":
                raise ValueError(f"Sample {sample!r}: target '*' is reserved for common controls")
            samples[sample] = (target, condition)
            arms.setdefault((target, condition), []).append(sample)
            if condition == "perturbed" and target not in seen_targets:
                targets.append(target)
                seen_targets.add(target)
    if not targets:
        raise ValueError("Samples TSV requires at least one perturbed target")
    for target in targets:
        if (target, "control") not in arms and ("*", "control") not in arms:
            raise ValueError(f"Target {target!r} has no control samples")
    return samples, targets, arms


def _read_counts(path, samples):
    genes = []
    seen_genes = set()
    counts = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 2 or header[0] != "gene":
            raise ValueError("Counts TSV requires gene followed by sample columns")
        sample_names = _labels(header[1:], "Counts samples")
        if set(sample_names) != set(samples):
            raise ValueError("Counts sample columns must match samples TSV assignments exactly")
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Counts TSV row {line} has an incorrect number of fields")
            gene = row[0]
            if not gene.strip():
                raise ValueError(f"Counts TSV row {line} has an empty gene symbol")
            if gene in seen_genes:
                raise ValueError(f"Counts TSV contains duplicate gene symbol {gene!r}")
            try:
                values = np.asarray(row[1:], dtype=np.float64)
            except ValueError as exc:
                raise ValueError(f"Counts TSV row {line} requires numeric counts") from exc
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"Counts TSV row {line} requires finite nonnegative counts")
            seen_genes.add(gene)
            genes.append(gene)
            counts.append(values)
    if not genes:
        raise ValueError("Counts TSV requires at least one gene")
    return genes, sample_names, np.stack(counts)


def _arm_mean(counts, indices, sample_names, target, condition, normalization):
    values = counts[:, indices]
    with np.errstate(over="ignore", invalid="ignore"):
        if normalization == "pooled-cp10k":
            pooled = values.sum(axis=1, dtype=np.float64)
            total = pooled.sum(dtype=np.float64)
            if not np.isfinite(total) or total <= 0:
                raise ValueError(
                    f"Target {target!r} {condition}: pooled library total must be positive and finite"
                )
            scaled = 10000.0 * pooled
            if np.isfinite(scaled).all():
                return scaled / total
            return (pooled / total) * 10000.0
        totals = values.sum(axis=0, dtype=np.float64)
        invalid = np.flatnonzero(~np.isfinite(totals) | (totals <= 0))
        if invalid.size:
            sample = sample_names[indices[int(invalid[0])]]
            raise ValueError(
                f"Sample {sample!r}: library total must be positive and finite for library-mean-cp10k"
            )
        return ((values / totals) * 10000.0).mean(axis=1, dtype=np.float64)


def estimate_source(
    counts_path: Path, samples_path: Path, output: Path, *,
    normalization: str = "pooled-cp10k", pseudocount: float = 0.1,
) -> None:
    output = Path(output)
    _check_output(output)
    if normalization not in ("pooled-cp10k", "library-mean-cp10k"):
        raise ValueError("normalization must be pooled-cp10k or library-mean-cp10k")
    if isinstance(pseudocount, (bool, np.bool_, str, bytes)):
        raise ValueError("pseudocount must be positive and finite")
    try:
        pseudocount = float(pseudocount)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("pseudocount must be positive and finite") from exc
    if not np.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("pseudocount must be positive and finite")
    samples, targets, arms = _read_samples(Path(samples_path))
    genes, sample_names, counts = _read_counts(Path(counts_path), samples)
    sample_index = {sample: column for column, sample in enumerate(sample_names)}
    effects = np.empty((len(targets), len(genes)), dtype=np.float64)
    for row, target in enumerate(targets):
        perturbed = arms[(target, "perturbed")]
        controls = arms.get((target, "control"), arms.get(("*", "control"), []))
        means = [
            _arm_mean(
                counts, [sample_index[sample] for sample in selected], sample_names,
                target, condition, normalization,
            )
            for selected, condition in ((perturbed, "perturbed"), (controls, "control"))
        ]
        case = means[0] + pseudocount
        control = means[1] + pseudocount
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            effect = np.log2(case / control)
        invalid = ~np.isfinite(effect)
        if invalid.any():
            effect[invalid] = np.log2(case[invalid]) - np.log2(control[invalid])
        effects[row] = effect
    if not np.isfinite(effects).all():
        raise ValueError("Estimated effects must be finite float64 values")
    _save(
        output, genes=np.asarray(genes, dtype=str), targets=np.asarray(targets, dtype=str),
        effects=effects,
    )
