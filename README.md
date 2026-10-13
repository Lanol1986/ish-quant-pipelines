# Quantitative ISH / RNAscope analysis framework

Reusable R pipelines for quantitative in-situ hybridisation and RNAscope analyses in oligonucleotide-therapeutic studies.

## Core rules

- **Animal is the biological replicate.** Cells, ROIs, sections and assay panels do not increase biological n.
- **Single-panel analysis and multi-panel integration are separate.** Pipelines 01-04 validate one physical assay panel. Pipeline 05 integrates already-validated standardized outputs.
- **Measurement semantics are explicit.** Genuine integer puncta counts, continuous HALO estimated-spot indices, area fractions and intensities are different endpoint types. Continuous estimates are never rounded/floored to make a count model run.
- **Panel identity is preserved.** Consecutive sections from the same animal remain separate measurements in integrated figures.
- **QC flags trigger review, not automatic exclusion.**
- **Stored normalized values are never clipped.**
- **Spatial maps use cell centroids and magma by default.** Consecutive sections are never treated as registered cell-to-cell maps.

## Pipelines

| Pipeline | Purpose |
|---|---|
| `pipelines/01_biodistribution_chromogenic.R` | Single-panel compound/oligonucleotide distribution |
| `pipelines/02_target_kd_chromogenic.R` | Single-panel target-mRNA pharmacodynamics |
| `pipelines/03_multiplex_celltype.R` | Single-panel compound/siRNA + target + cell-type analysis |
| `pipelines/04_compound_distribution_ihc.R` | Single-panel compound + IHC phenotype analysis |
| `pipelines/05_integrated_multipanel_analysis.R` | Integrated study-level analysis across validated panels |

## SRS4119

SRS4119 is the reference multi-panel study.

- Panel 1: Hprt1 siRNA passenger-strand signal + Hprt1 mRNA + neurons + microglia.
- Panel 2: Hprt1 siRNA passenger-strand signal + Hprt1 mRNA + neurons + oligodendrocytes + endothelial cells.
- Panels are consecutive 5-µm sections from the same animals.
- Panel 1 and Panel 2 are **not independent replicates**.
- Shared endpoints are shown visually separated in the same plot.
- Panel-specific cell types retain panel provenance.
- The integrated pipeline never matches cells across consecutive sections.

Run each panel first, then integrate:

```bash
Rscript pipelines/03_multiplex_celltype.R config/SRS4119_panel1.yml
Rscript pipelines/03_multiplex_celltype.R config/SRS4119_panel2.yml
Rscript pipelines/05_integrated_multipanel_analysis.R config/SRS4119_integration.yml
```

## Endpoint types

`endpoints[].kind` controls what analysis is allowed:

- `integer_count` — genuine discrete puncta counts; a two-part count analysis may be fitted.
- `estimated_spot_index` — cluster/intensity-derived estimate; treat as continuous and never floor/round.
- `area_fraction`
- `integrated_intensity`
- `continuous`
- `proportion`

For genuine integer counts, the framework defines a two-part analysis as:
1. occurrence above the positivity threshold;
2. a **zero-truncated negative-binomial** model among positive observations.

A regular negative-binomial model with `ziformula` is zero-inflated, not a hurdle model.

## Standardized result table

Each panel exports an animal-level tidy table with fields such as:

`study_id, animal_id, group, panel, section_id, organ, region, subregion, layer, cell_type, endpoint, metric, raw_value, normalized_value, unit, n_cells`.

This is the interface to pipeline 05 and to future scRNA-seq, snRNA-seq, qPCR or LC-MS concordance analyses.

## Presentation versus scientific figures

Important endpoints can generate:

- `FULL_SCIENTIFIC`: unrestricted values, full annotation.
- `PRESENTATION_COMPACT`: slide-ready, minimal internal title, individual animal points, and a 0-100% display for mRNA remaining where desired.

The compact mRNA renderer does not change stored values. Out-of-range points are flagged.

## Spatial analysis

Spatial plots use actual cell centroids, `magma`, and configurable percentile display limits. Display limits are recorded in metadata.

## Migration status

The original root scripts are retained temporarily as **legacy reference code**. New development occurs under `R/`, `pipelines/`, `config/`, `docs/`, and `tests/`.

The legacy scripts should only be removed after real-study regression testing confirms that the modular framework reproduces validated outputs where the scientific estimand is unchanged.

## Validation

Run:

```bash
Rscript tests/run_tests.R
```

SRS4119 raw exports are not stored in this public repository. End-to-end numerical regression on the real study therefore remains an explicit pre-merge validation step in the study environment.
