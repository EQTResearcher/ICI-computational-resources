"""Quick-start example for the economic-crisis ICI application."""
import sys
from pathlib import Path
CODE_DIR = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE_DIR))
from parameter_mapping import EconomicICIParameters, compute_ici, classify_rt_band

def main():
    case = EconomicICIParameters(D=250.0, C=1.0e13, S=1.0e8, F=8.0e3, W=1.0e-5, M=1.0e4, case_id="example_economic_crisis_case", period="illustrative")
    result = compute_ici(case)
    print(f"{case.case_id}: ICI={result['ICI']:.3f}, R(t)={result['R_t']:.3f}, band={classify_rt_band(result['R_t'])}")

if __name__ == "__main__":
    main()
