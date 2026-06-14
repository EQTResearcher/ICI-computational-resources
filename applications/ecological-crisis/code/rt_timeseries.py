"""R(t) trajectory helpers for ecological-crisis analysis."""
from typing import Iterable, List, Dict
from parameter_mapping import EcologicalICIParameters, compute_ici, classify_rt_band

def compute_rt_series(records: Iterable[EcologicalICIParameters]) -> List[Dict]:
    out = []
    for record in records:
        result = compute_ici(record)
        out.append({'case_id': record.case_id, 'ecosystem_type': record.ecosystem_type, 'period': record.period, 'ICI': result['ICI'], 'R_t': result['R_t'], 'band': classify_rt_band(result['R_t']), 'baseline_log10_DCS': result['baseline_log10_DCS'], 'FWM': result['FWM'], 'DCS': result['DCS']})
    return out
