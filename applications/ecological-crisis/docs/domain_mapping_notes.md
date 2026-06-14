# Domain Mapping Notes: Ecological Crisis

## D: effective biodiversity and functional diversity
D should measure effective diversity, not raw species count alone. Candidate proxies include Shannon effective species number, phylogenetic diversity, functional-trait space volume, and ecological niche differentiation.

## C: effective biomass and ecological carrying capacity
C should measure active ecological capacity. Candidate proxies include biomass, effective interacting biomass, resource stock, habitat capacity, and corrected remote-sensing estimates. Inactive, dead, or non-participating biomass should be discounted.

## S: metabolic flux and energy-throughput rate
S should measure the rate of effective energy conversion and material cycling. Candidate proxies include NPP/GPP converted into consistent throughput units, eddy-covariance fluxes, respiration/metabolism rates, and nutrient-flow rates.

## F: feedback-loop density and food-web regulation
F should measure closed and functional feedback loops rather than static connectivity. Candidate proxies include strongly connected components in food webs, cycle-basis counts, predator-prey feedbacks, mutualistic feedbacks, microbial recycling loops, and negative-feedback regulation pathways.

## W: normalized interaction frequency
W measures the normalized refresh rate of ecological interactions. Candidate proxies include species-interaction frequency, pollination-network timing, trophic exchange cadence, seasonal interaction windows, and response delays.

## M: ecological memory
M measures genetic, seed-bank, hydrological, soil, microbial, and adaptive memory. It should capture information residence time and recoverability, not archive size alone.

## Formula discipline

Use only the linear-FWM formula version:

```text
R(t) = log10(F * W * M) - log10(D * C * S)
```

Do not use the obsolete `F^W * M` structure.
