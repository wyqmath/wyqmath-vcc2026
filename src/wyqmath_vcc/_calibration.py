from __future__ import annotations

import numpy as np

ROUNDS = 16


def mean_probability(x, depths=None):
    depths = x.sum(axis=1, dtype=np.int64) if depths is None else depths
    return (x.astype(np.float64) / depths[:, None]).mean(axis=0)


def target_probability(q0, d, totals):
    q = q0 * np.exp2(d.astype(np.float64))
    mass = q.sum()
    if not np.isfinite(q).all() or not np.isfinite(mass) or mass <= 0:
        raise ValueError("Expected target probability must have finite positive mass")
    removed = float(q[totals == 0].sum() / q.sum())
    q[totals == 0] = 0
    if q.sum() <= 0:
        raise ValueError("No target probability remains on positive-total genes")
    q /= q.sum()
    return q, removed


def prefix_take(capacity, total):
    capacity = np.asarray(capacity, np.int64)
    previous = np.cumsum(capacity, dtype=np.int64) - capacity
    return np.minimum(capacity, np.maximum(total - previous, 0))


def calibrate(original, target):
    x = np.asarray(original).copy()
    if x.ndim != 2 or x.shape[0] < 1 or x.dtype.kind not in "iu" or np.any(x < 0):
        raise ValueError("A nonnegative integer cell-by-gene block is required")
    n = x.shape[0]
    depths = x.sum(axis=1, dtype=np.int64)
    totals = x.sum(axis=0, dtype=np.int64)
    target = np.asarray(target)
    if target.shape != totals.shape or np.any(depths <= 0):
        raise ValueError("Positive cell depths and a gene-axis target are required")
    if not np.isfinite(target).all() or np.any(target < 0):
        raise ValueError("Finite nonnegative target required")
    if abs(float(target.sum()) - 1) > 1e-12 or np.any(target[totals == 0] != 0):
        raise ValueError("Normalized supported fixed target required")
    before = mean_probability(x, depths)
    gap = target - before
    ring = np.argsort(depths, kind="stable").tolist()
    if n % 2:
        ring.append(None)
    half = len(ring) // 2
    for _ in range(ROUNDS):
        for i, j in zip(ring[:half], ring[:half - 1:-1]):
            if i is None or j is None:
                continue
            if depths[i] > depths[j]:
                i, j = j, i
            if depths[i] == depths[j]:
                continue
            delta = (1.0 / float(depths[i]) - 1.0 / float(depths[j])) / float(n)
            inc = np.flatnonzero((gap >= delta) & (x[j] > 0))
            dec = np.flatnonzero((-gap >= delta) & (x[i] > 0))
            if not len(inc) or not len(dec):
                continue
            ai = np.floor(np.minimum(x[j, inc], gap[inc] / delta)).astype(np.int64)
            ad = np.floor(np.minimum(x[i, dec], -gap[dec] / delta)).astype(np.int64)
            total = min(int(ai.sum()), int(ad.sum()))
            if total <= 0:
                continue
            ti, td = prefix_take(ai, total), prefix_take(ad, total)
            if int(ti.sum()) != total or int(td.sum()) != total:
                raise AssertionError("Balanced integer transfer failed")
            x[i, inc] += ti
            x[j, inc] -= ti
            x[i, dec] -= td
            x[j, dec] += td
            gap[inc] -= ti * delta
            gap[dec] += td * delta
        gap = target - mean_probability(x, depths)
        ring = [ring[0], ring[-1]] + ring[1:-1]
    after = mean_probability(x, depths)
    if np.any(x < 0) or not np.array_equal(x.sum(0, dtype=np.int64), totals):
        raise AssertionError("Exact integer gene totals lost")
    if not np.array_equal(x.sum(1, dtype=np.int64), depths):
        raise AssertionError("Exact integer cell depths lost")
    if np.max(np.abs(after - target) - np.abs(before - target)) > 2e-13:
        raise AssertionError("A gene moved farther from the fixed CPM target")
    return x, ROUNDS
