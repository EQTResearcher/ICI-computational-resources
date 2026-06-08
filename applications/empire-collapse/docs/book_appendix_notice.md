# Short Appendix Notice for the Book

## Computational Resources and Online Repository

The full computational resources for this volume are maintained in the ICI computational repository:

```text
https://github.com/EQTResearcher/ICI-computational-resources
```

The resources for this book are located in:

```text
applications/empire-collapse/
```

This directory contains the domain-specific code, data templates, historical parameter records, reproducibility notes, validation hooks, and source appendix materials for the empire-collapse application of the ICI framework.

To keep the printed and electronic book readable, long code listings and extensive machine-readable tables are not reproduced in full in the book. The book retains the mathematical definitions, historical interpretation, and methodological explanation; the repository provides executable materials for replication and extension.

Canonical formulas used by this volume:

```text
ICI(t) = k * lg(D(t) * C(t) * S(t)) * (1 + sqrt(alpha * F(t) * W(t) * M(t) / FWM_h))
R(t) = lg(F(t) * W(t) * M(t)) - lg(D(t) * C(t) * S(t))
```

Version note: readers should use the release tag or commit hash associated with the edition of the book they are reading.
