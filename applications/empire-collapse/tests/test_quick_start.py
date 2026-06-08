import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1] / "code"))
from rt_timeseries import compute_rt, rt_band


def test_rt_basic():
    rt = compute_rt(120, 5e7, 1.5e6, 3000, 1e-5, 1e4)
    assert isinstance(rt, float)
    assert rt_band(rt) in {"resilient", "critical", "collapse-risk"}
