"""
R(t) trajectory helpers for AI architecture analysis.
"""

from typing import Iterable, List, Dict
from parameter_mapping import AIICIParameters, compute_ici, classify_rt_band


def compute_rt_series(records: Iterable[AIICIParameters]) -> List[Dict]:
    """Compute ICI/R(t) values for a sequence of AI architecture records."""
    out = []
    for record in records:
        result = compute_ici(record)
        out.append({
            "case_id": record.case_id,
            "architecture_type": record.architecture_type,
            "model_family": record.model_family,
            "ICI": result["ICI"],
            "R_t": result["R_t"],
            "band": classify_rt_band(result["R_t"]),
            "baseline_log10_DCS": result["baseline_log10_DCS"],
            "FWM": result["FWM"],
            "DCS": result["DCS"],
        })
    return out
