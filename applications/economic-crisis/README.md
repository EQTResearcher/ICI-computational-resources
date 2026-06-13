# ICI Application: Economic Crisis Analysis

This directory contains application-specific computational resources for the ICI volume on economic crisis analysis.

The book applies the Integrated Complexity Index (ICI) and the balance factor R(t) to major economic and financial crises. Long program listings, reproducibility scripts, data templates, and validation notes should be maintained here rather than printed in full in the book.

## Repository location

```text
ICI-computational-resources/applications/economic-crisis/
```

## Scope

- economic-domain mapping of D, C, S, F, W, and M;
- R(t) trajectory construction for economic and financial crises;
- crisis-band classification and early-warning interpretation;
- templates for crisis-specific parameter records;
- reproducibility notes for figures, tables, and statistical checks;
- migration of long source-program sections from the manuscript into version-controlled resources.

## Relationship to the shared ICI core

This application should call the shared ICI implementation maintained in the repository core. Local scripts here are lightweight domain adapters, not an independent theory or formula fork.

```text
ICI(t) = k * log10(D(t) * C(t) * S(t)) *
         (1 + sqrt(alpha * F(t) * W(t) * M(t) / FWM_h))

R(t) = log10(F(t) * W(t) * M(t)) - log10(D(t) * C(t) * S(t))
```

## Suggested workflow

1. Read `docs/domain_mapping_notes.md`.
2. Copy `templates/economic_parameter_record_template.csv`.
3. Fill crisis-specific values for D, C, S, F, W, and M.
4. Run `examples/quick_start_economic_crisis.py`.
5. Record sources, assumptions, and uncertainty bands in `docs/reproducibility_notes.md`.
6. Use `appendix/source_economic_crisis_v3.docx` only as the source manuscript archive.
