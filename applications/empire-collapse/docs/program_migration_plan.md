# Program Migration Plan for the Empire Collapse Volume

## Principle

The printed book should not carry long executable code blocks. Code should be migrated into this directory and referenced by path. The book should retain only the conceptual formula, short usage notes, and a repository notice.

## What to move out of the book

- Long Python or Julia implementations
- Monte Carlo, Bootstrap, and uncertainty-propagation code
- Bai-Perron / breakpoint-test implementation details
- Kalman filter, particle filter, or online calibration code
- CI/CD, Docker, and environment-locking content
- Long parameter tables and machine-readable case records
- Visualization scripts and figure-generation code

## What should remain in the book

- The ICI and R(t) formulas
- The historical meaning of D, C, S, F, W, and M
- A short description of R(t) warning bands
- One small illustrative table, if necessary
- A link to this directory
- Citation and version information

## Recommended replacement line in the manuscript

```text
The full code, data templates, validation scripts, and reproducibility notes for this volume are maintained in the online repository under applications/empire-collapse/.
```
