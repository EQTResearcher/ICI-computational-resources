"""Domain-specific parameter record validation for ICI empire-collapse cases."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

@dataclass
class EmpireICIRecord:
    case_id: str
    period_start: int
    period_end: int
    D: float
    C: float
    S: float
    F: float
    W: float
    M: float
    D_log_std: float = 0.3
    C_log_std: float = 0.3
    S_log_std: float = 0.4
    F_log_std: float = 0.4
    W_log_std: float = 0.4
    M_log_std: float = 0.5
    quality_rating: str = "C"
    notes: Optional[str] = None

    def validate(self) -> None:
        values = {"D": self.D, "C": self.C, "S": self.S, "F": self.F, "W": self.W, "M": self.M}
        bad = [k for k, v in values.items() if v <= 0]
        if bad:
            raise ValueError(f"ICI parameters must be positive: {bad}")
        if self.period_end < self.period_start:
            raise ValueError("period_end must be greater than or equal to period_start")
        if self.quality_rating not in {"A", "B", "C", "D"}:
            raise ValueError("quality_rating must be A, B, C, or D")


def to_parameter_tuple(record: EmpireICIRecord) -> tuple[float, float, float, float, float, float]:
    record.validate()
    return (record.D, record.C, record.S, record.F, record.W, record.M)
