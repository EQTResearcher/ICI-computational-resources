# ICI Application: AI Architecture and Complexity Analysis

This directory contains the application-specific computational resources for the ICI volume on AI complexity and AI architecture analysis.

The book applies the Integrated Complexity Index (ICI) and the balance factor R(t) to AI systems, including large language models, Transformer architectures, recurrent networks, diffusion models, neuromorphic systems, quantum-classical hybrids, and future adaptive architectures. Long code blocks, reproducibility scripts, data templates, and validation notes should be maintained here rather than printed in full in the book.

## Repository location

```text
ICI-computational-resources/applications/ai-architecture-analysis/
```

## Scope

This application directory is intended to support:

- AI-domain mapping of the six ICI parameters: D, C, S, F, W, and M;
- ICI/R(t) computation for AI architectures and model families;
- architecture-complexity comparison across model generations;
- static-versus-dynamic architecture diagnostics;
- templates for AI-system parameter records;
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
2. Copy `templates/ai_parameter_record_template.csv`.
3. Fill architecture-specific values for D, C, S, F, W, and M.
4. Run `examples/quick_start_ai_architecture.py`.
5. Record sources, assumptions, uncertainty bands, and extraction choices in `docs/reproducibility_notes.md`.
6. Use `appendix/source_ai_complexity_v3.docx` only as the source manuscript archive.

## Status

Initial application-directory scaffold for the AI complexity volume.
