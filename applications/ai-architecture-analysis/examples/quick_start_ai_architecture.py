"""
Quick-start example for the AI architecture ICI application.

Run from this directory:
    python examples/quick_start_ai_architecture.py
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
CODE_DIR = HERE.parents[1] / "code"
sys.path.insert(0, str(CODE_DIR))

from parameter_mapping import AIICIParameters, compute_ici, classify_rt_band


def main():
    # Illustrative placeholder values only. Replace with documented, sourced estimates.
    case = AIICIParameters(
        D=1.0e4,          # computational primitive / functional diversity
        C=1.0e12,         # effective model/data capacity
        S=1.0e10,         # effective parallel processing throughput
        F=1.0e5,          # attention/recurrent/feedback integration density
        W=1.0e-6,         # normalized online-refresh frequency
        M=1.0e12,         # effective memory complexity
        case_id="example_transformer_case",
        architecture_type="transformer",
        model_family="illustrative",
    )
    result = compute_ici(case)
    print(f"{case.case_id}: ICI={result['ICI']:.3f}, R(t)={result['R_t']:.3f}, "
          f"band={classify_rt_band(result['R_t'])}")


if __name__ == "__main__":
    main()
