# Statistical notes

## Biological replication
The animal is the biological replicate. Cells, ROIs, sections and assay panels are nested/repeated measurements and do not increase biological n.

## Endpoint semantics before models
Each endpoint is explicitly typed. Continuous HALO cluster-derived estimates are not genuine counts and are never rounded or floored to satisfy a count-model requirement.

## Two-part count analysis
For genuine integer counts only, the framework separates:
1. occurrence above a declared threshold;
2. the positive-count distribution among positive observations.

The positive-count component uses a zero-truncated negative-binomial family. A regular negative-binomial model with `ziformula` is zero-inflated and should not be described as a hurdle model.

## Small n
With few animals, individual animal points, effect size, consistency and sensitivity analyses are primary. Very large cell counts do not create a larger biological n.

## Normalization
Target-mRNA values can be expressed as percent remaining relative to the matched control mean. The reference is computed within the configured panel × anatomical × cell-type stratum. Stored values are never clipped.

## Global summaries
`equal_region` and pooled-cell summaries are different estimands. The integration pipeline labels equal-region summaries explicitly.

## Multi-panel studies
Consecutive sections from the same animal are repeated measurements, not independent replicates. Main figures preserve panel identity. Panel averaging, if ever used, must be an explicitly labeled sensitivity/consensus analysis.
