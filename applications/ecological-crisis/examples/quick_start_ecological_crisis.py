"""Quick-start example for the ecological-crisis ICI application."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1] / 'code'))
from parameter_mapping import EcologicalICIParameters, compute_ici, classify_rt_band

def main():
    case = EcologicalICIParameters(D=500.0, C=1.0e9, S=1.0e5, F=1.0e3, W=1.0e-6, M=1.0e5, case_id='example_coral_reef_case', ecosystem_type='coral_reef', period='illustrative')
    result = compute_ici(case)
    print(f"{case.case_id}: ICI={result['ICI']:.3f}, R(t)={result['R_t']:.3f}, band={classify_rt_band(result['R_t'])}")
if __name__ == '__main__': main()
