from __future__ import annotations

from collections.abc import Iterator
from operator import index

import numpy as np
from scipy import sparse

from ._calibration import calibrate, target_probability
from ._sampling import sample_pair

__all__ = ["SOURCES", "predict_context"]

SOURCES = (
    "gw", "hct116", "hek293t", "cd4_Rest", "cd4_Stim8hr",
    "cd4_Stim48hr", "dld1v2", "h1scmg",
)

_SCALE = {"A": 67.77051441605022, "B": 24.612817569855043}
_A_DOSE = 0.8379070724891897


def _effects(raw, admitted, present, branch):
    nt, ng = raw.shape[1:]
    weights = np.asarray([1, 1, 1, 3, 3, 3, 1, 1] if branch == "A" else [1] * 8,
                         np.float32)
    old = np.zeros((nt, ng), np.float32 if branch == "A" else np.float64)
    nonzero = np.zeros((nt, ng), np.uint8)
    opposed = nonzero.copy()
    for ti in range(nt):
        ids = np.flatnonzero(admitted[:, ti])
        if not len(ids):
            raise ValueError(f"No admitted source for target index {ti}; no control fallback")
        rows = raw[ids, ti]
        if not np.isfinite(rows).all():
            raise ValueError(f"Nonfinite admitted source row for target index {ti}")
        if np.any(rows[~present[ids]] != 0):
            raise ValueError("Unmapped entries in admitted source rows must be zero")
        v = (rows * weights[ids, None]).astype(np.float64)
        if not np.isfinite(v).all():
            raise ValueError("Nonfinite weighted source row")
        norms = np.linalg.norm(v, axis=1)
        active = norms > 1e-12
        u = np.zeros_like(v)
        u[active] = v[active] / norms[active, None]

        cons = u.mean(0) * _SCALE[branch]
        if branch == "A":
            original32 = (0.85 * cons).astype(np.float32)
            original32[~(np.abs(cons) > 0.0)] = 0
            old[ti] = (original32.astype(np.float64) * _A_DOSE).astype(np.float32)
        else:
            keep = np.abs(cons) > 0.0
            old[ti, keep] = 0.85 * cons[keep]
        for j, si in enumerate(ids):
            vote = present[si] & active[j] & (raw[si, ti] != 0)
            nonzero[ti] += vote.astype(np.uint8)
            opposed[ti] += (vote & (np.sign(raw[si, ti]) == -np.sign(old[ti])) &
                            (old[ti] != 0)).astype(np.uint8)
    gate = np.zeros((nt, ng), bool)
    new = old.copy()
    if branch == "A":
        gate = ((nonzero >= 3) &
                (3 * opposed.astype(np.int16) >= 2 * nonzero.astype(np.int16)) & (old != 0))
        new[gate] = (old[gate].astype(np.float64) * 0.5).astype(np.float32)
    return old, new, gate


def _make_query(counts, cells):
    lam = np.maximum(np.asarray(counts.mean(axis=0)).ravel(), 0.0)
    xf = sparse.csr_matrix(counts, dtype=np.float64)
    depths = np.asarray(xf.sum(1)).ravel()
    if not np.isfinite(depths).all() or not (depths > 0).all():
        raise ValueError("Positive finite control cell depths required")
    q0 = np.asarray(xf.multiply(1 / depths[:, None]).sum(0)).ravel() / counts.shape[0]
    selected = np.sort(np.random.Generator(np.random.PCG64(1)).choice(
        counts.shape[0], cells, replace=False))
    template = counts[selected].toarray().astype(np.float64)
    template = 0.45 * template + (1.0 - 0.45) * template.mean(0)[None, :]
    if (not np.isfinite(lam).all() or not np.isfinite(q0).all() or
            abs(float(q0.sum()) - 1) > 1e-12 or not np.isfinite(template).all()):
        raise ValueError("Finite control mean, template, and normalized probability required")
    return lam, q0, template


def predict_context(
    counts: sparse.csr_matrix,
    raw: np.ndarray,
    admitted: np.ndarray,
    present: np.ndarray,
    branch: str,
    *,
    cells_per_target: int = 400,
) -> Iterator[sparse.csr_matrix]:
    if not isinstance(branch, str) or branch not in ("A", "B"):
        raise ValueError("branch must be 'A' or 'B'")
    if isinstance(cells_per_target, (bool, np.bool_)):
        raise ValueError("cells_per_target must be a positive integer")
    try:
        cells = index(cells_per_target)
    except TypeError as exc:
        raise ValueError("cells_per_target must be a positive integer") from exc
    if cells < 1:
        raise ValueError("cells_per_target must be a positive integer")
    if not isinstance(raw, np.ndarray) or raw.dtype != np.float32 or raw.ndim != 3 or raw.shape[0] != 8:
        raise ValueError("raw must be a float32 ndarray of shape (8, targets, genes)")
    nt, ng = raw.shape[1:]
    if ng < 1:
        raise ValueError("At least one gene is required")
    if not isinstance(admitted, np.ndarray) or admitted.dtype != bool or admitted.shape != (8, nt):
        raise ValueError("admitted must be a bool ndarray of shape (8, targets)")
    if not isinstance(present, np.ndarray) or present.dtype != bool or present.shape != (8, ng):
        raise ValueError("present must be a bool ndarray of shape (8, genes)")
    if not sparse.issparse(counts) or counts.format != "csr":
        raise ValueError("counts must be a scipy CSR matrix or array")
    if counts.shape[0] < cells or counts.shape[1] != ng:
        raise ValueError("Control counts need enough cells and the same gene axis as raw")
    if (counts.dtype.kind not in "fiu" or not np.isfinite(counts.data).all() or
            (counts.data < 0).any() or (counts.data != np.floor(counts.data)).any()):
        raise ValueError("Control counts must be finite nonnegative integers")
    counts = sparse.csr_matrix(counts)
    old, new, gate = _effects(raw, admitted, present, branch)
    lam, q0, template = _make_query(counts, cells)
    rng = np.random.Generator(np.random.PCG64(0))
    for ti in range(nt):
        sampled, totals = sample_pair(old[ti], new[ti], gate[ti], lam, template, rng, branch)
        target, _ = target_probability(q0, new[ti], totals)
        depths = sampled.sum(1, dtype=np.int64)
        if not (depths > 0).all():
            raise ValueError("Positive sampled cell depths required")
        result = calibrate(sampled, target)[0] if branch == "A" else sampled
        after_depths = result.sum(1, dtype=np.int64)
        if not np.array_equal(result.sum(0, dtype=np.int64), totals):
            raise AssertionError("Exact integer gene totals lost")
        if not np.array_equal(after_depths, depths):
            raise AssertionError("Exact sampled cell depths lost")
        if not ((after_depths > 0) & (after_depths < 1_000_000)).all():
            raise ValueError("Output cell depths must be positive and below 1000000")
        yield sparse.csr_matrix(result)
