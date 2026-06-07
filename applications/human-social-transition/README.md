# ICI Domain Application: Human Social Transition

This directory contains the computational resources for the ICI domain-application volume on human social transition.

The purpose of this application directory is to separate book-specific resources from the shared ICI core framework. The shared formula implementation, SOPs, validation logic, and common templates should remain in the top-level `core/` and `shared/` directories of the main repository. This directory should only contain resources specific to the human-social-transition volume.

## Repository path

```text
applications/human-social-transition/
```

## Recommended contents

```text
human-social-transition/
├── README.md
├── code/
│   ├── README.md
│   ├── parameter_mapping.py
│   ├── rt_timeseries.py
│   └── validation_hooks.py
├── data/
│   ├── README.md
│   └── human_social_transition_cases_template.csv
├── examples/
│   ├── README.md
│   └── quick_start_human_transition.py
├── docs/
│   ├── README.md
│   ├── domain_mapping_notes.md
│   ├── reproducibility_notes.md
│   └── book_appendix_notice.md
├── appendix/
│   ├── README.md
│   ├── source_appendix_system.docx
│   ├── source_appendix_system.md
│   └── short_appendix_repository_notice.docx
├── figures/
│   └── README.md
├── templates/
│   ├── README.md
│   └── parameter_record_template.csv
└── tests/
    ├── README.md
    └── test_basic_imports.py
```

## Scope

This volume applies the ICI framework to large-scale transitions in human society. Typical resources in this directory may include:

- human-society parameter mapping rules;
- historical or socio-institutional case templates;
- R(t) time-series analysis scripts;
- crisis-threshold and transition-warning examples;
- book-specific documentation and reproducibility notes.

## Dependency on the shared ICI core

This directory should not duplicate the core ICI formula. All executable examples should import the standard implementation from the shared core package once it is available in the repository.

The standard formula used by the ICI framework is:

```text
ICI = k * lg(D*C*S) * (1 + sqrt(alpha * F*W*M / FWM_h))
```

The older formula based on `sqrt(D*C*S*F^W*M)` should not be used.

## Suggested citation

For now, cite the main repository:

```text
Li, K., & Li, L. ICI Computational Resources. GitHub repository:
https://github.com/EQTResearcher/ICI-computational-resources
```

When a Zenodo DOI is created, replace this temporary citation with the DOI-based citation.
