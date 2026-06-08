"""Quick-start example for applications/empire-collapse."""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1] / "code"))

from rt_timeseries import compute_ici, compute_rt, rt_band

case = {
    "D": 120,
    "C": 5e7,
    "S": 1.5e6,
    "F": 3e3,
    "W": 1e-5,
    "M": 1e4,
}

rt = compute_rt(**case)
ici = compute_ici(**case)
print(f"Example empire case: ICI={ici['ICI']:.3f}, R(t)={rt:.3f}, band={rt_band(rt)}")
