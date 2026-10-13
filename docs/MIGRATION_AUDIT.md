# Audit-driven migration

This branch implements the high-priority findings from the repository audit.

## Corrected in the modular framework
1. Two-part count analysis uses a zero-truncated NB positive component rather than calling a conventional zero-inflated NB model a hurdle model.
2. Continuous estimated spot indices are never floored/rounded to make count models run.
3. Schema-v1 configuration is operational and separates single-panel analysis from multipanel integration.
4. Missing-value handling is explicit per endpoint.
5. Multipanel integration is first-class and does not increase biological n.
6. Region/subregion/layer anatomy is configurable.
7. Spatial RNA plots use magma and percentile-based display limits.
8. QC flags are machine-readable and do not automatically exclude REVIEW samples.
9. Provenance is saved alongside input hashes and exact config.
10. Legacy root scripts remain temporarily for traceability.

## Pre-merge requirement
Run the modular pipeline on the real SRS4119 raw exports and compare all validated numerical endpoints against the accepted v23 results. Any expected difference caused by a corrected estimand/model must be documented before merge.
