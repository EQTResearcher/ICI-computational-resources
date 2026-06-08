"""Placeholder interfaces for structural-break diagnostics.

This application can later connect to Bai-Perron or equivalent breakpoint tests.
Keep statistical implementation details here rather than in the printed book.
"""
from __future__ import annotations
from dataclasses import dataclass

@dataclass
class BreakpointDiagnostic:
    case_id: str
    method: str
    break_year: int | None
    p_value: float | None
    notes: str = ""


def placeholder_breakpoint_diagnostic(case_id: str, notes: str = "Not yet estimated") -> BreakpointDiagnostic:
    return BreakpointDiagnostic(case_id=case_id, method="placeholder", break_year=None, p_value=None, notes=notes)
