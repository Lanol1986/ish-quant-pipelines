# Integrated multipanel analysis

Pipeline 05 is intentionally separate from the single-panel assay pipelines.

Each physical panel has its own staining, classifier, section, and QC. These are validated first. Pipeline 05 then consumes standardized animal-level outputs.

## Rules
- Panel identity is never discarded in the main visualization.
- Consecutive sections do not increase biological n.
- Shared endpoints can be connected within animal to show repeated-panel concordance.
- Panel-specific cell types remain panel-specific; missing measurements are never invented.
- Consecutive sections are not cell-to-cell registered unless a separate explicit registration method is performed.

## SRS4119
Panel 1 contains neurons and microglia; Panel 2 contains neurons, oligodendrocytes and endothelial cells. Hprt1 siRNA passenger-strand signal and Hprt1 mRNA are shared endpoints.

The integration figure therefore shows Panel 1 and Panel 2 separately in the same plot, matching the intended presentation style.
