# ICI Application: Empire Collapse

This directory contains the application-specific computational resources for the ICI volume on empire collapse.

The book applies the Integrated Complexity Index (ICI) and the balance factor R(t) to historical imperial systems. Long code blocks, data templates, validation scripts, and reproducibility notes should be maintained here rather than printed in full in the book.

## Repository location

```text
applications/empire-collapse/
```

## Core formula

```text
ICI(t) = k * lg(D(t) * C(t) * S(t)) * (1 + sqrt(alpha * F(t) * W(t) * M(t) / FWM_h))
R(t) = lg(F(t) * W(t) * M(t)) - lg(D(t) * C(t) * S(t))
```

This application uses the shared ICI core implementation in the repository root or `core/` directory. Application files here should not redefine the canonical formula unless explicitly marked as a test, wrapper, or domain adapter.

## Directory map

```text
code/       Domain-specific parameter mapping and R(t) utilities
examples/   Minimal runnable examples for empire-collapse cases
data/       Case templates and reproducible input tables
docs/       Domain mapping notes, reproducibility notes, migration plan
appendix/   Source manuscript appendix/material and extracted markdown
figures/    Generated figures and exported charts
templates/  Parameter record templates
tests/      Basic tests for application scripts
```

## Recommended workflow

1. Use `data/empire_collapse_cases_template.csv` to enter historical case parameters.
2. Use `code/parameter_mapping.py` to validate six-parameter records.
3. Use `code/rt_timeseries.py` to compute R(t) and warning bands.
4. Use `code/breakpoint_notes.py` as the placeholder interface for Bai-Perron / structural-break validation.
5. Keep book-facing explanations in `docs/book_appendix_notice.md`; keep executable content in `code/` and `examples/`.

## Citation

Please cite the book volume and this repository directory together when using the empire-collapse application resources.
