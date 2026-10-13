# QC reference

`QC_STATUS.csv` uses:
- **PASS** — no configured concern detected.
- **REVIEW** — inspect the source assay/image before strong interpretation.
- **FAIL** — a precondition for valid analysis is not satisfied.
- **NOT_RUN** — not applicable or required inputs were unavailable.

A REVIEW flag never automatically excludes an animal.

## Required QC concepts
- animal/file roster completeness;
- explicit missing-value handling;
- invalid/negative endpoint values;
- phenotype ambiguity/unclassified cells;
- cell-proportion QC by animal and panel;
- marker intensity separation when intensities are supplied;
- coordinate completeness for spatial output;
- positivity-threshold sensitivity;
- panel concordance for shared measurements.

The cell-proportion stacked bar is an ROI/segmentation QC and must not be described as a true whole-brain cell-composition estimate.
