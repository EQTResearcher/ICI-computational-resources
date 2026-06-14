# ICI Application: Ecological Crisis Analysis

This directory contains the application-specific computational resources for the ICI volume on ecological and environmental crisis analysis.

The book applies the Integrated Complexity Index (ICI) and the balance factor R(t) to ecosystems, biodiversity loss, ecological regime shifts, feedback collapse, critical slowing down, and resilience reconstruction. Long code blocks, reproducibility scripts, data templates, and validation notes should be maintained here rather than printed in full in the book.

## Repository location

```text
ICI-computational-resources/applications/ecological-crisis/
```

## Scope

This application directory is intended to support:

- ecological-domain mapping of the six ICI parameters: D, C, S, F, W, and M;
- R(t) trajectory construction for ecosystem degradation and regime-shift cases;
- early-warning interpretation using critical slowing down signals;
- templates for ecological parameter records;
- reproducibility notes for figures, tables, and statistical checks;
- migration of long source-program sections from the manuscript into version-controlled resources.

## Relationship to the shared ICI core

This application should call the shared ICI implementation maintained in the repository core. The local scripts here are lightweight domain adapters, not an independent theory or formula fork.

The standard formula used across the ICI application volumes is:

```text
ICI(t) = k * log10(D(t) * C(t) * S(t)) *
         (1 + sqrt(alpha * F(t) * W(t) * M(t) / FWM_h))

R(t) = log10(F(t) * W(t) * M(t)) - log10(D(t) * C(t) * S(t))
```

## Suggested workflow

1. Read `docs/domain_mapping_notes.md`.
2. Copy `templates/ecological_parameter_record_template.csv`.
3. Fill ecosystem-specific values for D, C, S, F, W, and M.
4. Run `examples/quick_start_ecological_crisis.py`.
5. Record sources, assumptions, uncertainty bands, and data-processing choices in `docs/reproducibility_notes.md`.
6. Use `appendix/source_ecological_crisis_v3.docx` only as the source manuscript archive.

## Status

Initial application-directory scaffold for the ecological crisis volume.
