"""R(t) helpers for human-social-transition case studies."""

import math
from typing import Iterable, Sequence


def compute_rt(D: float, C: float, S: float, F: float, W: float, M: float) -> float:
    """Compute R(t) = lg(F*W*M) - lg(D*C*S)."""
    params = [D, C, S, F, W, M]
    if any(x <= 0 for x in params):
        raise ValueError("All ICI parameters must be positive.")
    return math.log10(F * W * M) - math.log10(D * C * S)


def compute_rt_series(rows: Iterable[Sequence[float]]) -> list[float]:
    """Compute R(t) for an iterable of six-parameter rows."""
    return [compute_rt(*row) for row in rows]
