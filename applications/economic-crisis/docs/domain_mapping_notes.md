# Domain Mapping Notes: Economic Crisis

This file records how the six ICI parameters are mapped into the economic and financial-crisis domain.

## D: asset-category diversity and strategy-space width
D should measure economically meaningful diversity, not nominal product counts. Use independent asset classes, market-participant types, strategy clusters, or derivative-structure families only when they have distinct pricing, liquidity, or regulatory behavior.

## C: effective credit capacity and interaction scale
C should measure the effective resource and interaction capacity of the system. Candidate proxies include leverage-adjusted credit stock, liquid balance-sheet capacity, open interest, central-bank liquidity, core counterparty exposures, and payment-network scale.

## S: clearing throughput and information-processing rate
S should measure the rate at which the system processes economically meaningful transactions or information updates. Candidate proxies include payment-settlement throughput, trading-message flow, order-book refresh rate, clearing volume, or policy-operation frequency.

## F: feedback-loop density and regulatory/market topology
F should measure closed and functional feedback loops, not merely institution count. Candidate proxies include leverage-price loops, margin-call loops, central-bank reaction loops, dealer-inventory feedback, cross-asset arbitrage loops, and regulatory countercyclical loops.

## W: normalized interaction frequency and response speed
W is the normalized response or refresh frequency of feedback loops. It captures whether feedback can close in time. A high-S system with low W can become unstable because market throughput outruns regulatory or institutional response.

## M: institutional memory and expectation anchoring
M measures stored, accessible, and operational crisis memory. Candidate proxies include regulatory frameworks, stress-test history, crisis-resolution playbooks, historical risk-model depth, and the persistence of policy credibility.

## Formula discipline
Use only the linear-FWM formula version:

```text
R(t) = log10(F * W * M) - log10(D * C * S)
```
