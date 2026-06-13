"""Validation hooks for the economic-crisis ICI application."""
def check_required_fields(record):
    required = ["case_id", "D", "C", "S", "F", "W", "M", "source_notes"]
    missing = [k for k in required if k not in record or record[k] in ("", None)]
    if missing:
        raise ValueError(f"Missing required fields: {missing}")

def formula_version_guard(formula_version: str) -> None:
    if formula_version != "linear_FWM_v4":
        raise ValueError("Formula mismatch. This application uses linear_FWM_v4.")
