"""Economic-domain parameter mapping for the ICI framework."""
from dataclasses import dataclass
from math import log10, sqrt
from typing import Optional, Dict

K_NORM = 1.259
ALPHA = 1.02e5
FWM_H = 7.52e11

@dataclass
class EconomicICIParameters:
    D: float  # asset / strategy diversity
    C: float  # effective credit capacity or interaction scale
    S: float  # clearing / information throughput rate
    F: float  # closed feedback-loop density
    W: float  # normalized interaction / response frequency
    M: float  # institutional memory / resilience buffer
    case_id: str = "unnamed_case"
    period: Optional[str] = None
    notes: str = ""

def validate_positive(params: EconomicICIParameters) -> None:
    for name in ("D", "C", "S", "F", "W", "M"):
        value = getattr(params, name)
        if value <= 0:
            raise ValueError(f"{name} must be positive; got {value}.")

def compute_ici(params: EconomicICIParameters, k: float = K_NORM, alpha: float = ALPHA, fwm_h: float = FWM_H) -> Dict[str, float]:
    validate_positive(params)
    dcs = params.D * params.C * params.S
    fwm = params.F * params.W * params.M
    baseline = log10(dcs)
    emergence = sqrt(alpha * fwm / fwm_h)
    ici = k * baseline * (1.0 + emergence)
    rt = log10(fwm) - log10(dcs)
    return {"ICI": ici, "R_t": rt, "baseline_log10_DCS": baseline, "FWM": fwm, "DCS": dcs, "emergence": emergence}

def classify_rt_band(rt: float) -> str:
    if rt > -9: return "resilient"
    if rt > -10: return "watch"
    if rt > -11: return "critical"
    return "collapse-risk"
