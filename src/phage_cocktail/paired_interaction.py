"""Convex, explicit bilinear host–phage logistic models without outer-row storage."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import roc_auc_score

from .paired_baseline import weighted_ap

KINDS = ("phage_only", "additive", "bilinear")


def parameter_size(kind, dh, dp):
    if kind not in KINDS:
        raise ValueError("unsupported interaction model")
    return 1 + dp + (dh if kind != "phage_only" else 0) + (dh * dp if kind == "bilinear" else 0)


def unpack(theta, kind, dh, dp):
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (parameter_size(kind, dh, dp),) or not np.isfinite(theta).all():
        raise ValueError("coefficient dimensions or values differ from the model")
    start, a = 1, np.zeros(dh)
    if kind != "phage_only":
        a, start = theta[start:start + dh], start + dh
    b, start = theta[start:start + dp], start + dp
    W = theta[start:].reshape(dh, dp) if kind == "bilinear" else None
    return theta[0], a, b, W


def loss_gradient(theta, H, P, hi, pi, y, weights, kind, l2):
    """Exact weighted logistic objective; sparse observed cells supply gradients."""
    c, a, b, W = unpack(theta, kind, H.shape[1], P.shape[1])
    z = c + (H @ a)[:, None] + (P @ b)[None, :]
    if W is not None:
        z += (H @ W) @ P.T
    values = z[hi, pi]
    total = weights.sum()
    loss = float(weights @ (np.logaddexp(0., values) - y * values) / total
                 + .5 * l2 * np.sum(theta[1:] ** 2))
    residual = weights * (expit(values) - y) / total
    cells = np.zeros((len(H), len(P)))
    np.add.at(cells, (hi, pi), residual)
    blocks = [np.array([residual.sum()])]
    if kind != "phage_only":
        blocks.append(H.T @ cells.sum(axis=1))
    blocks.append(P.T @ cells.sum(axis=0))
    if W is not None:
        blocks.append((H.T @ cells @ P).ravel())
    gradient = np.concatenate(blocks)
    gradient[1:] += l2 * theta[1:]
    return loss, gradient


def fit_interaction(H, P, hi, pi, y, weights, *, kind, l2, optimizer):
    """Fit slopes and interactions; require convergence under the frozen criteria."""
    H, P = np.asarray(H, dtype=float), np.asarray(P, dtype=float)
    hi, pi, y, weights = map(np.asarray, (hi, pi, y, weights))
    if (H.ndim != 2 or P.ndim != 2 or not H.size or not P.size
            or not np.isfinite(H).all() or not np.isfinite(P).all()
            or hi.dtype.kind not in "iu" or pi.dtype.kind not in "iu"
            or y.ndim != 1 or any(x.ndim != 1 or len(x) != len(y) for x in (hi, pi, weights))
            or not len(y) or set(y) != {0, 1} or not np.isfinite(weights).all()
            or (weights <= 0).any() or hi.min() < 0 or hi.max() >= len(H)
            or pi.min() < 0 or pi.max() >= len(P) or not np.isfinite(l2) or l2 <= 0):
        raise ValueError("invalid fit arrays, labels, weights or regularization")
    initial = np.zeros(parameter_size(kind, H.shape[1], P.shape[1]))
    prior = (weights @ y + .5) / (weights.sum() + 1)
    initial[0] = np.log(prior / (1 - prior))
    options = {k: optimizer[k] for k in ("maxiter", "maxls", "maxcor", "ftol", "gtol")}
    result = minimize(loss_gradient, initial, args=(H, P, hi, pi, y, weights, kind, l2),
                      jac=True, method="L-BFGS-B", options=options)
    value, gradient = loss_gradient(result.x, H, P, hi, pi, y, weights, kind, l2)
    maximum = float(np.max(np.abs(gradient)))
    if (not result.success or not np.isfinite(value) or not np.isfinite(gradient).all()
            or maximum > optimizer["maximum_admitted_gradient_infinity_norm"]):
        raise ValueError(f"optimizer not admitted: {kind} L2={l2}, {result.message}, grad={maximum}")
    return dict(kind=kind, l2=l2, theta=result.x, iterations=int(result.nit),
                objective=value, gradient_infinity_norm=maximum, optimizer_message=str(result.message))


def predict_logits(model, H, P):
    """Fixed feature-order sums preserve exact ties across entity/batch ordering."""
    H, P = np.asarray(H, dtype=float), np.asarray(P, dtype=float)
    if H.ndim != 2 or P.ndim != 2 or not np.isfinite(H).all() or not np.isfinite(P).all():
        raise ValueError("finite two-dimensional host/phage features required")
    c, a, b, W = unpack(model["theta"], model["kind"], H.shape[1], P.shape[1])
    base = np.full(len(H), c)
    for j, value in enumerate(a):
        base += H[:, j] * value
    conditional = np.broadcast_to(b, (len(H), len(b))).copy()
    if W is not None:
        for j in range(H.shape[1]):
            conditional += H[:, j, None] * W[j, None, :]
    z = np.broadcast_to(base[:, None], (len(H), len(P))).copy()
    for j in range(P.shape[1]):
        z += conditional[:, j, None] * P[None, :, j]
    if not np.isfinite(z).all():
        raise ValueError("nonfinite prediction")
    return z


def fit_entity_scaler(features, ids, digests):
    """Entity statistics use unique fit assemblies, never held-out rows."""
    unique = {}
    for name in sorted(set(ids)):
        digest = digests[name]
        if digest in unique and not np.array_equal(features.loc[name], features.loc[unique[digest]]):
            raise ValueError("identical assembly has inconsistent features")
        unique.setdefault(digest, name)
    values = features.loc[list(unique.values())].to_numpy(float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("missing or nonfinite fit features")
    scales = values.std(axis=0)
    scales[scales == 0] = 1.
    return dict(columns=features.columns.tolist(), means=values.mean(axis=0).tolist(),
                scales=scales.tolist(), unique_fit_assemblies=len(unique))


def transform_entities(features, scaler):
    if features.columns.tolist() != scaler["columns"]:
        raise ValueError("feature schema/order differs from the fitted scaler")
    values = (features.to_numpy(float) - scaler["means"]) / scaler["scales"]
    if not np.isfinite(values).all():
        raise ValueError("missing or nonfinite transformed features")
    return pd.DataFrame(values, index=features.index, columns=features.columns)


def classification_report(rows, logits, weights):
    y, logits, weights = rows.infection.to_numpy(), np.asarray(logits), np.asarray(weights)
    both = set(y) == {0, 1}
    return dict(rows=len(rows), positives=int(y.sum()), coverage=1.,
                weighted_total=float(weights.sum()), weighted_positive=float(weights @ y),
                average_precision=weighted_ap(y, logits, weights) if both else None,
                auroc=float(roc_auc_score(y, logits, sample_weight=weights)) if both else None,
                uncalibrated_brier=float(np.average((expit(logits) - y) ** 2, weights=weights)),
                log_loss=float(np.average(np.logaddexp(0., logits) - y * logits, weights=weights)))
