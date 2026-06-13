import sys
from pathlib import Path
CODE_DIR = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE_DIR))
from parameter_mapping import EconomicICIParameters, compute_ici, classify_rt_band

def test_compute_ici_smoke():
    p = EconomicICIParameters(D=10, C=1000, S=100, F=10, W=0.01, M=100)
    r = compute_ici(p)
    assert "ICI" in r
    assert "R_t" in r
    assert classify_rt_band(r["R_t"]) in {"resilient", "watch", "critical", "collapse-risk"}
