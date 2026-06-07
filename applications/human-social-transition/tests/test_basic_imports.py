from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from rt_timeseries import compute_rt
from parameter_mapping import HumanSocialTransitionRecord, validate_positive_parameters


def test_compute_rt():
    assert isinstance(compute_rt(100, 1_000_000, 10_000, 100, 0.0001, 1000), float)


def test_positive_validation():
    record = HumanSocialTransitionRecord("x", 0, 1, 1, 1, 1, 1, 1, 1)
    validate_positive_parameters(record)
