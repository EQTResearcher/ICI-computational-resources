"""Minimal R(t) helpers for empire-collapse application records."""
from __future__ import annotations
import math
from typing import Iterable

K_NORM = 1.259
ALPHA = 1.02e5
FWM_H = 7.52e11


def compute_ici(D: float, C: float, S: float, F: float, W: float, M: float) -> dict:
    if min(D, C, S, F, W, M) <= 0:
        raise ValueError("All ICI parameters must be positive")
    baseline = math.log10(D * C * S)
    fwm = F * W * M
    emergence = 1.0 + math.sqrt(ALPHA * fwm / FWM_H)
    return {
        "ICI": K_NORM * baseline * emergence,
        "baseline": baseline,
        "FWM": fwm,
        "emergence": emergence,
        "fwm_ratio": fwm / FWM_H,
    }


def compute_rt(D: float, C: float, S: float, F: float, W: float, M: float) -> float:
    if min(D, C, S, F, W, M) <= 0:
        raise ValueError("All ICI parameters must be positive")
    return math.log10(F * W * M) - math.log10(D * C * S)


def rt_band(rt: float) -> str:
    if rt > -9:
        return "resilient"
    if rt > -11:
        return "critical"
    return "collapse-risk"


def compute_series(rows: Iterable[dict]) -> list[dict]:
    out = []
    for row in rows:
        vals = {k: float(row[k]) for k in ["D", "C", "S", "F", "W", "M"]}
        rt = compute_rt(**vals)
        ici = compute_ici(**vals)
        out.append({**row, **ici, "R_t": rt, "R_band": rt_band(rt)})
    return out
