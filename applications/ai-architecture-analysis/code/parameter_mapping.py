"""
AI-domain parameter mapping for the ICI framework.

This module provides lightweight helper functions for mapping AI architectures and
model systems into the six ICI parameters. It is not a substitute for the shared
ICI core implementation.
"""

from dataclasses import dataclass
from math import log10, sqrt
from typing import Optional, Dict


K_NORM = 1.259
ALPHA = 1.02e5
FWM_H = 7.52e11


@dataclass
class AIICIParameters:
    """Six-parameter container for an AI architecture case."""
    D: float  # independent computational primitive / functional-type diversity
    C: float  # effective model/data capacity after information-bottleneck correction
    S: float  # effective parallel processing throughput
    F: float  # directed feedback-loop or recurrent/attention integration density
    W: float  # normalized state-refresh / online-update frequency
    M: float  # effective memory complexity / parameter-state residence depth
    case_id: str = "unnamed_model"
    architecture_type: Optional[str] = None
    model_family: Optional[str] = None
    notes: str = ""


def validate_positive(params: AIICIParameters) -> None:
    for name in ("D", "C", "S", "F", "W", "M"):
        value = getattr(params, name)
        if value <= 0:
            raise ValueError(f"{name} must be positive; got {value}.")


def compute_ici(params: AIICIParameters,
                k: float = K_NORM,
                alpha: float = ALPHA,
                fwm_h: float = FWM_H) -> Dict[str, float]:
    """Compute ICI and core components using the shared v4 linear-FWM structure."""
    validate_positive(params)
    dcs = params.D * params.C * params.S
    fwm = params.F * params.W * params.M
    baseline = log10(dcs)
    emergence = sqrt(alpha * fwm / fwm_h)
    ici = k * baseline * (1.0 + emergence)
    rt = log10(fwm) - log10(dcs)
    return {
        "ICI": ici,
        "R_t": rt,
        "baseline_log10_DCS": baseline,
        "FWM": fwm,
        "DCS": dcs,
        "emergence": emergence,
    }


def classify_rt_band(rt: float) -> str:
    """Classify the R(t) value into a broad AI-complexity risk band."""
    if rt > -9:
        return "resilient"
    if rt > -10:
        return "watch"
    if rt > -11:
        return "critical"
    return "overload-risk"
