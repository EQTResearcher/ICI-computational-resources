# Domain Mapping Notes: AI Architecture Analysis

This file records how the six ICI parameters are mapped into the AI architecture and AI-complexity domain.

## D: computational primitive and functional-type diversity

D should measure independent functional diversity rather than raw parameter count. Candidate proxies include computational graph primitives, distinct operator families, modality-specific modules, attention/recurrence/state-space components, and independent functional clusters.

## C: effective capacity

C should measure effective model/data capacity, not raw dataset size or raw parameter count alone. Candidate proxies include active parameter count, information-bottleneck-corrected training corpus capacity, context capacity, latent-state capacity, and effective representational volume.

## S: effective processing throughput

S should measure realized processing throughput, not theoretical peak FLOPs alone. Candidate proxies include roofline-corrected throughput, token/s, inference/training throughput, memory-bandwidth-adjusted compute, and effective parallel update rate.

## F: feedback and integration-loop density

F should measure functional closed loops or effective integration pathways. Candidate proxies include directed cycles in computation graphs, recurrent loops, attention-mediated integration paths, tool-use feedback, self-reflection loops, and controller-feedback loops.

## W: normalized state-refresh and online adaptation frequency

W measures whether feedback and memory can update in time. Static pretrained systems typically have very low W during deployment, while online-learning, event-driven, recurrent, or neuromorphic systems may have higher W.

## M: memory complexity

M measures effective memory complexity, not raw stored bytes alone. Candidate proxies include effective parameter memory, persistent state, context memory, retrieval memory, weight-update capacity, and Hessian- or redundancy-corrected representational memory.

## Formula discipline

Use only the linear-FWM formula version:

```text
R(t) = log10(F * W * M) - log10(D * C * S)
```

Do not use the obsolete `F^W * M` structure.
