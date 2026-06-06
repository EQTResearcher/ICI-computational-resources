# ICI Computational Resources

Companion computational resources for the ICI framework. This repository is designed to replace a long printed programming appendix with an online, maintainable, and citable resource package.

Repository URL after creation:

```text
https://github.com/EQTResearcher/ICI-computational-resources
```

## What is included

```text
src/ici/                         Python package
  core.py                        ICI and R(t) calculation
  csd.py                         Critical slowing down diagnostics
  crossdomain.py                 KSG/MI and Procrustes validation tools
  calibration.py                 Bayesian calibration utilities
  visualization.py               Plotting helpers

examples/                        Reproducible examples
  quick_start.py                 Minimal example
  rome_analysis.py               Historical-domain time-series example

data/                            Sample parameter database files
  cs_ici_db_sample.csv           Small sample from Appendix B

templates/                       Data-entry templates
  parameter_record_template.csv  Standard parameter-record form

docs/                            Documentation converted from the appendix
  appendix_A.md                  Computational stack and SOPs
  appendix_B.md                  Cross-species parameter database
  appendix_C.md                  Historical-domain methodology
  appendix_D.md                  Economic/ecological mapping rules
  appendix_E.md                  Bayesian calibration derivation
  appendix_F.md                  Falsification and failure protocols
  appendix_A_F_full.md           Full converted appendix text
  SOP-01_parameter_extraction.md Parameter extraction protocol
  SOP-02_rt_timeseries_monitoring.md R(t) monitoring protocol
  repository_usage_cn.md         Chinese usage notes

archive/
  source_appendix_A_F.docx       Original Word appendix supplied for reference
```

## Installation

Clone the repository and install the local package:

```bash
git clone https://github.com/EQTResearcher/ICI-computational-resources.git
cd ICI-computational-resources
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
pip install -e .
```

Alternatively, with Conda:

```bash
conda env create -f environment.yml
conda activate ici-computational-resources
pip install -e .
```

## Quick start

```python
from ici import ICIParameters, compute_ici

params = ICIParameters.from_domain_defaults(
    domain="biological_lab",
    D=2.0e3,
    C=1.0e8,
    S=1.0e5,
    F=20,
    W=1.0e-3,
    M=1.0e2,
    system_name="E. coli K-12",
)

result = compute_ici(params)
print(result.summary())
```

Or run:

```bash
python examples/quick_start.py
python examples/rome_analysis.py
```

## Core formula

The package implements:

```text
DCS = D × C × S
FWM = F × W × M
ICI = k × lg(DCS) × [1 + sqrt(alpha × FWM / FWM_h)]
R(t) = lg(FWM) - lg(DCS)
```

Default constants:

```text
k = 1.259
alpha = 1.02 × 10^5
FWM_h = 7.52 × 10^11
```

## Documentation status

The Markdown documentation in `docs/` is converted from the supplied Word appendix and should be reviewed before formal publication. The Python modules in `src/ici/` have been cleaned into executable files for repository use.

## Citation

Please cite this repository using `CITATION.cff`. After creating a release, consider archiving it through Zenodo to obtain a DOI.

## License

Code is released under the MIT License. Book text and conceptual exposition may be subject to separate copyright terms; adjust this section if you want stricter control over non-code materials.
