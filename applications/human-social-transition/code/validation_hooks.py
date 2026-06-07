"""Validation hooks for the human-social-transition application.

These checks are deliberately conservative. They are intended to catch obvious
parameter-entry or formula-version errors before case studies are published.
"""


def check_formula_version(formula_version: str) -> None:
    expected = "linear_FWM"
    if expected not in formula_version:
        raise ValueError(
            f"Unexpected formula version: {formula_version!r}. "
            "Use the linear F*W*M formulation."
        )
