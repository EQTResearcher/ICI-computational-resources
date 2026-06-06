"""Core ICI calculations.

This module implements the six-parameter ICI calculation described in the
companion appendix. It supports Monte Carlo uncertainty propagation and R(t)
calculation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional
import warnings

import numpy as np

ICI_CONSTANTS = {
    "k": 1.259,
    "alpha": 1.02e5,
    "FWM_h": 7.52e11,
    "R_critical": -10.0,
    "R_resilient": -9.0,
    "R_collapse": -11.0,
}

DOMAIN_UNCERTAINTY = {
    "biological_lab": {"D": 0.02, "C": 0.04, "S": 0.13, "F": 0.09, "W": 0.11, "M": 0.17},
    "biological_field": {"D": 0.06, "C": 0.11, "S": 0.22, "F": 0.26, "W": 0.17, "M": 0.35},
    "ecological": {"D": 0.06, "C": 0.11, "S": 0.22, "F": 0.26, "W": 0.17, "M": 0.35},
    "economic_modern": {"D": 0.04, "C": 0.06, "S": 0.17, "F": 0.22, "W": 0.13, "M": 0.26},
    "historical_empire": {"D": 0.17, "C": 0.22, "S": 0.30, "F": 0.35, "W": 0.30, "M": 0.39},
    "ai_architecture": {"D": 0.04, "C": 0.02, "S": 0.04, "F": 0.22, "W": 0.26, "M": 0.30},
}


@dataclass
class ICIParameters:
    """Six ICI parameters with log-scale uncertainty."""

    D: float
    C: float
    S: float
    F: float
    W: float
    M: float
    D_sigma: float = 0.10
    C_sigma: float = 0.10
    S_sigma: float = 0.15
    F_sigma: float = 0.20
    W_sigma: float = 0.20
    M_sigma: float = 0.25
    system_name: str = "unnamed"
    domain: str = "unknown"
    time_point: Optional[str] = None
    data_sources: Dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ["D", "C", "S", "F", "W", "M"]:
            val = getattr(self, name)
            if val <= 0:
                raise ValueError(f"Parameter {name} must be positive; got {val}.")
            if val < 1e-50:
                warnings.warn(f"Parameter {name} is extremely small ({val:.2e}); check normalization.")

    @classmethod
    def from_domain_defaults(cls, domain: str, **kwargs) -> "ICIParameters":
        uncertainty = DOMAIN_UNCERTAINTY.get(domain, DOMAIN_UNCERTAINTY["biological_field"])
        defaults = {f"{k}_sigma": v for k, v in uncertainty.items()}
        defaults.update(kwargs)
        defaults["domain"] = domain
        return cls(**defaults)


@dataclass
class ICIResult:
    """Result container for an ICI calculation."""

    ICI: float
    ICI_lower: float
    ICI_upper: float
    ICI_lower_95: float
    ICI_upper_95: float
    R_t: float
    R_t_lower: float
    R_t_upper: float
    baseline: float
    emergence: float
    alert_level: str
    params: ICIParameters

    def summary(self) -> str:
        return "\n".join(
            [
                f"System: {self.params.system_name}",
                f"Domain: {self.params.domain}",
                "",
                f"ICI = {self.ICI:.1f} [{self.ICI_lower:.1f}, {self.ICI_upper:.1f}] (68% CI)",
                f"      [{self.ICI_lower_95:.1f}, {self.ICI_upper_95:.1f}] (95% CI)",
                "",
                f"R(t) = {self.R_t:.2f} [{self.R_t_lower:.2f}, {self.R_t_upper:.2f}] (68% CI)",
                "",
                f"Baseline lg(DCS) = {self.baseline:.2f}",
                f"Emergence sqrt(alpha*FWM/FWM_h) = {self.emergence:.4f}",
                "",
                f"Alert level: {self.alert_level.upper()}",
            ]
        )


def determine_alert_level(Rt: float) -> str:
    """Determine alert level from R(t)."""
    if Rt > ICI_CONSTANTS["R_resilient"]:
        return "green"
    if Rt > -9.5:
        return "yellow"
    if Rt > ICI_CONSTANTS["R_collapse"]:
        return "orange"
    return "red"


def compute_ici(
    params: ICIParameters,
    n_samples: int = 10_000,
    reference_Rt: Optional[float] = None,
    seed: int = 42,
) -> ICIResult:
    """Compute ICI, R(t), uncertainty intervals, and alert level."""
    k = ICI_CONSTANTS["k"]
    alpha = ICI_CONSTANTS["alpha"]
    FWM_h = ICI_CONSTANTS["FWM_h"]

    DCS = params.D * params.C * params.S
    FWM = params.F * params.W * params.M
    baseline = np.log10(DCS)
    fwm_ratio = FWM / FWM_h
    emergence = np.sqrt(alpha * fwm_ratio) if fwm_ratio > 0 else 0.0
    ICI_point = k * baseline * (1.0 + emergence)
    R_t_abs = np.log10(FWM) - np.log10(DCS)
    R_t_point = R_t_abs if reference_Rt is None else R_t_abs - reference_Rt

    rng = np.random.default_rng(seed)
    log_D = rng.normal(np.log10(params.D), params.D_sigma, n_samples)
    log_C = rng.normal(np.log10(params.C), params.C_sigma, n_samples)
    log_S = rng.normal(np.log10(params.S), params.S_sigma, n_samples)
    log_F = rng.normal(np.log10(params.F), params.F_sigma, n_samples)
    log_W = rng.normal(np.log10(params.W), params.W_sigma, n_samples)
    log_M = rng.normal(np.log10(params.M), params.M_sigma, n_samples)

    DCS_s = (10 ** log_D) * (10 ** log_C) * (10 ** log_S)
    FWM_s = (10 ** log_F) * (10 ** log_W) * (10 ** log_M)
    baseline_s = np.log10(np.clip(DCS_s, 1e-300, None))
    ratio_s = np.clip(FWM_s / FWM_h, 0, None)
    ICI_s = k * baseline_s * (1.0 + np.sqrt(alpha * ratio_s))

    Rt_s = np.log10(np.clip(FWM_s, 1e-300, None)) - np.log10(np.clip(DCS_s, 1e-300, None))
    if reference_Rt is not None:
        Rt_s = Rt_s - reference_Rt

    pcts = np.percentile(ICI_s, [2.5, 16, 50, 84, 97.5])
    Rt_pcts = np.percentile(Rt_s, [16, 50, 84])

    return ICIResult(
        ICI=ICI_point,
        ICI_lower=pcts[1],
        ICI_upper=pcts[3],
        ICI_lower_95=pcts[0],
        ICI_upper_95=pcts[4],
        R_t=R_t_point,
        R_t_lower=Rt_pcts[0],
        R_t_upper=Rt_pcts[2],
        baseline=baseline,
        emergence=emergence,
        alert_level=determine_alert_level(Rt_pcts[1]),
        params=params,
    )


def compute_ici_timeseries(params_list: list[ICIParameters], reference_index: int = 0) -> list[ICIResult]:
    """Compute an ICI/R(t) time series using one point as the R(t) reference."""
    ref_result = compute_ici(params_list[reference_index])
    return [compute_ici(p, reference_Rt=ref_result.R_t) for p in params_list]
