from __future__ import annotations

import copy

import numpy as np


def totals_from_mu(mu):
    s = np.asarray(mu, dtype=np.float64).sum(axis=0)
    if not np.isfinite(s).all() or (s < 0).any():
        raise ValueError("Expected gene totals must be finite and nonnegative")
    rounded = np.maximum(np.rint(s), 0.0)
    if np.any(rounded > np.iinfo(np.int32).max):
        raise ValueError("Expected gene total exceeds int32 capacity")
    return rounded.astype(np.int64)


def allocation_probs(template, mu):
    mu = np.asarray(mu, dtype=np.float64)
    p = np.array(template, dtype=np.float64, copy=True)
    s = p.sum(axis=0)
    bad = ~(s > 0)
    if bad.any():
        p[:, bad] = mu[:, bad]
        s = p.sum(axis=0)
    bad = ~(s > 0)
    if bad.any():
        p[:, bad] = 1.0 / p.shape[0]
        s = p.sum(axis=0)
    p /= s
    sm = p[:-1].sum(axis=0)
    over = sm > 1.0
    if over.any():
        p[:-1, over] *= (1.0 - 8.0 * np.finfo(np.float64).eps) / sm[over]
    return p


def sample_pair(old, new, gate, lam, template, rng, branch):
    if branch == "A":
        mu = lam * np.exp2(old.astype(np.float64))
        mn = lam * np.exp2(new.astype(np.float64))
    else:
        mu = lam * (2.0 ** old.astype(np.float64))
        mn = lam * (2.0 ** new.astype(np.float64))
    mu2 = np.ascontiguousarray(np.broadcast_to(mu, template.shape), dtype=np.float64)
    mn2 = np.ascontiguousarray(np.broadcast_to(mn, template.shape), dtype=np.float64)
    if not np.isfinite(mu2).all() or not np.isfinite(mn2).all():
        raise ValueError("Nonfinite expected counts")
    totals, totals_new = totals_from_mu(mu2), totals_from_mu(mn2)
    prob, prob_new = allocation_probs(template, mu2), allocation_probs(template, mn2)
    baseline = np.zeros(template.shape, np.int32)
    sampled = baseline.copy()
    for g in np.flatnonzero((totals > 0) | ((totals_new > 0) & gate)):
        if gate[g]:
            fork = np.random.Generator(np.random.PCG64(0))
            fork.bit_generator.state = copy.deepcopy(rng.bit_generator.state)
        if totals[g] > 0:
            baseline[:, g] = rng.multinomial(int(totals[g]), prob[:, g])
        if gate[g]:
            if totals_new[g] > 0:
                sampled[:, g] = fork.multinomial(int(totals_new[g]), prob_new[:, g])
        else:
            sampled[:, g] = baseline[:, g]
    if not np.array_equal(baseline.sum(0, dtype=np.int64), totals):
        raise AssertionError("Original integer gene totals lost")
    if not np.array_equal(sampled.sum(0, dtype=np.int64), totals_new):
        raise AssertionError("Changed integer gene totals lost")
    return sampled, totals_new
