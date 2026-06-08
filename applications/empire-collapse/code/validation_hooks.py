"""Validation hooks for empire-collapse application data."""
from __future__ import annotations

def check_formula_version(version: str) -> None:
    if version != "linear_FWM_v4":
        raise ValueError("Use formula version linear_FWM_v4 for this volume")


def check_rt_plausibility(rt: float) -> str:
    if rt < -20 or rt > 0:
        return "check: R(t) is outside the typical empire-domain range"
    return "ok"
