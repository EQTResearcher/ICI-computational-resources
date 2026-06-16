# Reproducibility Notes: AI Architecture Analysis

For every AI architecture case, record:

1. model or architecture name;
2. architecture type and model family;
3. parameter estimates for D, C, S, F, W, and M;
4. source notes for each parameter;
5. uncertainty ranges or log-standard deviations;
6. formula version;
7. ICI/R(t) computation method;
8. code version or commit hash;
9. any deviations from the shared ICI core implementation.

## Minimum metadata

```text
case_id
architecture_type
model_family
D, C, S, F, W, M
D_log_std, C_log_std, S_log_std, F_log_std, W_log_std, M_log_std
source_notes
formula_version = linear_FWM_v4
```
