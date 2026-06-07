"""Domain-specific parameter mapping for the human-social-transition ICI volume.

This module is intentionally lightweight. It should call the shared ICI core
implementation from the main repository once the package path is finalized.
"""

from dataclasses import dataclass


@dataclass
class HumanSocialTransitionRecord:
    """Minimal template for a human-social-transition ICI case."""

    case_id: str
    period_start: int
    period_end: int
    D: float
    C: float
    S: float
    F: float
    W: float
    M: float
    quality_rating: str = "draft"
    notes: str = ""


def validate_positive_parameters(record: HumanSocialTransitionRecord) -> None:
    """Validate that all six ICI parameters are strictly positive."""
    for name in ("D", "C", "S", "F", "W", "M"):
        value = getattr(record, name)
        if value <= 0:
            raise ValueError(f"{name} must be positive; got {value!r} for {record.case_id}")
