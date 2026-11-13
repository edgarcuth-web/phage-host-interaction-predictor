"""Analytic gradients, genuine pairing behavior and fit-only preprocessing."""
import numpy as np
import pandas as pd
import pytest
from scipy.optimize._numdiff import approx_derivative
from scipy.special import expit

from phage_cocktail.paired_interaction import (
    classification_report, fit_entity_scaler, fit_interaction, loss_gradient,
    parameter_size, predict_logits, transform_entities,
)

OPTIMIZER = dict(maxiter=2000, maxls=40, maxcor=20, ftol=1e-12, gtol=1e-6,
                 maximum_admitted_gradient_infinity_norm=1e-4)


@pytest.mark.parametrize("kind", ["phage_only", "additive", "bilinear"])
def test_matrix_free_objective_and_gradient_match_explicit_outer_features(kind):
    rng = np.random.default_rng(7)
    H, P = rng.normal(size=(3, 2)), rng.normal(size=(4, 3))
    hi, pi = np.array([0, 1, 2, 0, 2, 1, 1]), np.array([0, 1, 2, 3, 0, 1, 3])
    y, weights = np.array([1., 0., 1., 0., 0., 1., 1.]), np.array([1., .5, 1., 1., 2., .5, 1.])
    theta = rng.normal(size=parameter_size(kind, 2, 3))
    blocks = [np.ones((len(hi), 1))]
    if kind != "phage_only": blocks.append(H[hi])
    blocks.append(P[pi])
    if kind == "bilinear": blocks.append(np.einsum("ij,ik->ijk", H[hi], P[pi]).reshape(len(hi), -1))
    X = np.column_stack(blocks)
    z, penalty = X @ theta, theta.copy()
    penalty[0] = 0.
    expected = weights @ (np.logaddexp(0., z) - y*z) / weights.sum() + .05 * np.sum(penalty**2)
    gradient = X.T @ (weights*(expit(z) - y)) / weights.sum() + .1 * penalty
    actual, actual_gradient = loss_gradient(theta, H, P, hi, pi, y, weights, kind, .1)
    assert actual == pytest.approx(expected, abs=1e-12)
    np.testing.assert_allclose(actual_gradient, gradient, atol=1e-12, rtol=0)
    finite_difference = approx_derivative(lambda t: loss_gradient(t, H, P, hi, pi, y, weights, kind, .1)[0], theta).ravel()
    np.testing.assert_allclose(actual_gradient, finite_difference, atol=1e-8, rtol=0)


def test_bilinear_learns_pairing_that_balanced_additive_margins_cannot_learn():
    H, P = np.array([[-1.], [1.]]), np.array([[-1.], [1.]])
    hi, pi, y = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1]), np.array([1., 0., 0., 1.])
    fitted = {kind: fit_interaction(H, P, hi, pi, y, np.ones(4), kind=kind, l2=.01, optimizer=OPTIMIZER)
              for kind in ("additive", "bilinear")}
    additive, bilinear = (predict_logits(fitted[k], H, P) for k in ("additive", "bilinear"))
    np.testing.assert_allclose(additive, 0., atol=1e-10)
    assert bilinear[0, 0] > bilinear[0, 1] and bilinear[1, 1] > bilinear[1, 0]
    assert fitted["bilinear"]["objective"] < fitted["additive"]["objective"]


def test_predictions_preserve_exact_feature_ties_under_batching_and_reordering():
    rng = np.random.default_rng(3)
    H, P = rng.normal(size=(4, 5)), rng.normal(size=(3, 6))
    H[3], P[2] = H[0], P[0]
    model = dict(kind="bilinear", theta=rng.normal(size=parameter_size("bilinear", 5, 6)))
    z = predict_logits(model, H, P)
    assert np.array_equal(z[0], z[3]) and np.array_equal(z[:, 0], z[:, 2])
    assert np.array_equal(z[[2, 0]], predict_logits(model, H[[2, 0]], P))


def test_scaler_uses_unique_fit_assemblies_and_ignores_extreme_held_out_values():
    values = pd.DataFrame([[1., 2.], [1., 2.], [3., 2.], [10000., -10000.]],
                          index=["a", "a_alias", "b", "held"], columns=["x", "y"])
    digests = dict(a="same", a_alias="same", b="other", held="held")
    scaler = fit_entity_scaler(values, ["a", "a_alias", "b"], digests)
    assert scaler["unique_fit_assemblies"] == 2
    assert scaler["means"] == [2., 2.] and scaler["scales"] == [1., 1.]
    changed = values.copy(); changed.loc["held"] *= 100
    assert scaler == fit_entity_scaler(changed, ["a", "a_alias", "b"], digests)
    with pytest.raises(ValueError, match="schema"):
        transform_entities(values[["y", "x"]], scaler)
    changed.loc["a_alias", "x"] += 1
    with pytest.raises(ValueError, match="inconsistent"):
        fit_entity_scaler(changed, ["a", "a_alias", "b"], digests)


@pytest.mark.parametrize("failure", ["weight", "label", "index", "nan"])
def test_invalid_fit_inputs_are_not_silently_imputed(failure):
    H, P, hi, pi = np.ones((2, 2)), np.ones((2, 2)), np.array([0, 1]), np.array([0, 1])
    y, weights = np.array([0., 1.]), np.ones(2)
    if failure == "weight": weights[0] = 0
    if failure == "label": y[0] = 2
    if failure == "index": hi[0] = -1
    if failure == "nan": H[0, 0] = np.nan
    with pytest.raises(ValueError, match="invalid"):
        fit_interaction(H, P, hi, pi, y, weights, kind="bilinear", l2=.1, optimizer=OPTIMIZER)


def test_raw_logit_metrics_do_not_merge_saturated_sigmoid_ranks():
    rows = pd.DataFrame(dict(infection=[1, 0]))
    result = classification_report(rows, np.array([1000., 999.]), np.ones(2))
    assert expit(1000.) == expit(999.) == 1.
    assert result["average_precision"] == 1. and np.isfinite(result["log_loss"])
