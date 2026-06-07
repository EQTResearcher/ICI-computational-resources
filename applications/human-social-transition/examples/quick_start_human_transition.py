"""Minimal quick-start example for the human-social-transition application."""

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from rt_timeseries import compute_rt  # noqa: E402


def main() -> None:
    data_path = ROOT / "data" / "human_social_transition_cases_template.csv"
    with data_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rt = compute_rt(
                float(row["D"]), float(row["C"]), float(row["S"]),
                float(row["F"]), float(row["W"]), float(row["M"]),
            )
            print(f"{row['case_id']}: R(t) = {rt:.3f}")


if __name__ == "__main__":
    main()
