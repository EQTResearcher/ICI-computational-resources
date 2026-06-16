"""
Validation hooks for the AI architecture ICI application.
"""

from typing import Dict


def check_required_fields(record: Dict) -> None:
    required = ["case_id", "architecture_type", "D", "C", "S", "F", "W", "M", "source_notes"]
    missing = [k for k in required if k not in record or record[k] in ("", None)]
    if missing:
        raise ValueError(f"Missing required fields: {missing}")


def formula_version_guard(formula_version: str) -> None:
    if formula_version != "linear_FWM_v4":
        raise ValueError(
            "Formula mismatch. This application uses linear_FWM_v4: "
            "ICI = k*log10(DCS)*(1+sqrt(alpha*FWM/FWM_h)); "
            "R(t)=log10(FWM)-log10(DCS)."
        )


def ai_falsification_flag(rt_value: float, irreversible_collapse: bool) -> str:
    """Simple domain-level guard for the AI R(t) threshold claim."""
    if rt_value > -9 and irreversible_collapse:
        return "CHECK_THRESHOLD_CLAIM"
    return "OK"
