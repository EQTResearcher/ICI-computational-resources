"""Bayesian calibration tools for ICI constants."""
from __future__ import annotations

from typing import List

import numpy as np


def calibrate_ici_constants(reference_systems: List[dict], n_samples: int = 2000, n_chains: int = 4, target_accept: float = 0.9) -> dict:
    """Calibrate k and alpha using PyMC.

    PyMC and ArviZ are imported lazily so that the main package remains usable
    without the Bayesian stack installed.
    """
    import pymc as pm
    import arviz as az

    DCS_vals = np.array([np.log10(s["D"] * s["C"] * s["S"]) for s in reference_systems])
    FWM_ratio_vals = np.array([s["F"] * s["W"] * s["M"] / 7.52e11 for s in reference_systems])
    ICI_obs = np.array([s["ICI_observed"] for s in reference_systems])
    ICI_sigma = np.array([s.get("ICI_sigma", 0.1) for s in reference_systems])

    with pm.Model() as model:
        k_param = pm.Normal("k", mu=1.26, sigma=0.1)
        log_alpha = pm.Normal("log_alpha", mu=np.log(1.02e5), sigma=0.5)
        alpha_param = pm.Deterministic("alpha", pm.math.exp(log_alpha))
        emergence = pm.math.sqrt(alpha_param * FWM_ratio_vals)
        ICI_pred = k_param * DCS_vals * (1 + emergence)
        pm.Normal("obs", mu=pm.math.log(pm.math.abs(ICI_pred) + 1e-10), sigma=ICI_sigma, observed=np.log(ICI_obs))
        trace = pm.sample(draws=n_samples, chains=n_chains, target_accept=target_accept, return_inferencedata=True)

    k_post = trace.posterior["k"].values.flatten()
    alpha_post = trace.posterior["alpha"].values.flatten()
    rhat = az.rhat(trace)
    ess = az.ess(trace)
    converged = all(rhat.data_vars[v].values.max() < 1.01 for v in ["k", "log_alpha"]) and all(
        ess.data_vars[v].values.min() > 400 for v in ["k", "log_alpha"]
    )
    return {
        "k_mean": float(np.mean(k_post)),
        "k_std": float(np.std(k_post)),
        "k_ci95": (float(np.percentile(k_post, 2.5)), float(np.percentile(k_post, 97.5))),
        "alpha_mean": float(np.mean(alpha_post)),
        "alpha_std": float(np.std(alpha_post)),
        "alpha_ci95": (float(np.percentile(alpha_post, 2.5)), float(np.percentile(alpha_post, 97.5))),
        "converged": bool(converged),
        "r_hat_max": float(max(rhat.data_vars[v].values.max() for v in ["k", "log_alpha"])),
        "trace": trace,
    }
