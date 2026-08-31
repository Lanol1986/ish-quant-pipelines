# QC reference

What each QC output means and what to do about it. Read `qc/` before `tables/`,
every time.

---

## Present in every pipeline

### `qc/file_manifest.csv`
Every input file, its animal id, size and modification time. Check the animal
count matches the design before anything else.

The roster assertion runs at the same point and **stops the pipeline** if an
animal in the config has no file, or a file has no config entry. To exclude an
animal, put it in `design$expected_absent` with a reason. Silent absence of an
animal is one of the easiest ways to change a result.

### `qc/read_qc.csv`
Per file: rows read, cells kept, and the number of NA values coerced to zero.
`n_na_est` should be small. Anything above half the cells in a file means the
export settings were wrong for that slide, and the pipeline warns.

Parsing problems are not tolerated at all. A misparsed column becomes character,
then NA, then zero, and the slide reads as having no signal.

### `qc/expression_regime.csv`
Per stratum: percent positive, the fraction of positive cells above the
countable range, and a `regime` label.

| regime | Meaning | Action |
|---|---|---|
| `countable` | Per-cell counts are defensible | proceed |
| `mixed` | A minority of positive cells exceed the countable range | % positive as primary; counts reported as an index |
| `dense_or_saturated` | Values are cluster-derived, not counts | switch primary to % labelled area or integrated intensity |

`is_integer_valued` tells you whether the column is a genuine count. If it is
FALSE, the column is a model output built on an assumed single-dot value, and
that divisor belongs in `readout$single_dot_calibration`.

### `qc/threshold.yml`
The threshold derived from control material, with the rule used and the number
of non-zero background values it rests on. Compare against the value in the
config; if they differ substantially, decide which to use and record why.

### `qc/threshold_sweep.csv` and `plots/*_threshold_sweep.png`
Percent positive per animal at thresholds 1 through 6. Flat group separation
across the sweep means the threshold is not driving the result. Collapsing
separation means the result depends on a number that was chosen rather than
measured.

### `qc/cluster_rate.csv`
Percent of cells containing at least one cluster, per animal and region, with
`fold_vs_group_median` and a `CHECK STAINING` flag at tenfold. The cluster rate
is a sensitive indicator of signal density, so two animals in one group that
differ by an order of magnitude are not staining replicates of each other.

### `qc/staining.csv` and `plots/*_staining_vs_readout.png`
Per-slide median and 90th-percentile optical density or intensity. Plot the
primary readout against these before believing between-animal spread is
biological. If the two are correlated, the spread is at least partly chemistry.

If this table is empty, the export contained no optical-density or intensity
columns. Adding them costs nothing and is the fastest way to separate biology
from staining.

### `qc/denominators.csv`
Cells and cells/mm2 per stratum. A per-area readout is confounded by
cellularity, so cells/mm2 travels beside it. A per-cell readout depends entirely
on segmentation, so the cell count has to be visible.

### `qc/degenerate_check.csv`
Distinct values taken by each summary statistic across animals. Anything flagged
`DEGENERATE` takes two or fewer distinct values and must not be used as an
endpoint. This is what catches a median-derived quantity that has collapsed to a
coin flip.

### `qc/design_power.yml`
Animals per group and the smallest achievable two-sided Wilcoxon p. If
`inference_possible` is FALSE, the run is descriptive: no boxplots, no group
contrasts.

### `qc/unresolved_regions.csv`
Annotations that matched none of the declared regions. Those cells are excluded,
and this is where you see how many. Add entries to `design$region_aliases` if
annotation names vary.

### `tables/null_effect_summary.csv`
Not strictly QC, but read it with the QC. It is the apparent effect the pipeline
produces from the control animals alone, using a leave-one-out reference. Any
treated effect smaller than this is not established.

### `run_manifest/run_manifest.yml`
`unrecorded_parameters` lists every config field left blank. Fill them in and
rerun.

---

## Pipeline 03 only: dual-positive characterisation

A cell scored positive for two supposedly exclusive markers arises three ways,
and they have separable signatures.

| Check | Output | Reads as a segmentation artefact when |
|---|---|---|
| Rate | `qc/dual_rate.csv` | exceeds `celltype$dual_rate_alarm_pct` |
| Intensity signature | `plots/QC_dual_intensity_signature.png` | duals form a low shoulder between the negatives and the true positives |
| Compartment logic | `qc/dual_compartments.csv` | marker B is cytoplasm-positive and nucleus-negative on a marker-A nucleus-positive cell |
| Density dependence | `models/dual_density_model.csv` | odds ratio per SD of local density well above 1 |
| Densest region | `qc/dual_by_region.csv` | the most densely packed region has by far the highest rate |

A large density coefficient means the fix is upstream — a smaller expansion
radius, or protein-defined boundaries — not a filter applied afterwards. The
definitive test is a dilation sweep: re-run the classification at two or three
radii and plot dual rate against radius. Linear in radius means geometry.

The alternative explanations are genuine biological association (satellite
cells, phagocytosed material), which is weakly density-dependent, and spectral
bleed-through or autofluorescence, which correlates across all channels rather
than two.

`celltype$dual_policy` declares in advance what happens to them: `report_only`,
`own_class`, or `exclude` (recording how many per group per region). There is no
silent option, because dropping duals into an unused class and renormalising the
denominator makes every composition percentage depend on how many duals the
slide produced.

### `qc/composition.csv` and `qc/marker_intensity_drift.csv`
Per-animal marker-positive fractions and per-slide intensity distributions. If a
marker positivity threshold moves between slides, the single-positive rate, the
dual rate and the composition denominator all move with it. Wide spread in the
control group is the signal to look at the staining before anything else.

---

## Pipeline 04 only: phenotype partition and marker failure

### `qc/phenotype_partition.csv`
Do the four mutually exclusive phenotype counts sum to the total cell count? If
not, either a phenotype is missing from the export or the phenotypes overlap,
and every percentage computed from them is wrong.

### `qc/marker_failure.csv` and `plots/05_marker_detection_qc.png`
A stratum where an antibody yielded zero positive cells, while every other
animal in that region yielded thousands, is a staining or classification
failure. It is flagged and every marker-derived quantity there is masked to NA.

This matters because the failure is otherwise invisible. A percentage divides by
the marker count, gets zero, becomes NA and vanishes from its panel with a
"Removed n rows" warning. A density divides by area, gets a clean zero, and is
**plotted on the axis**, where it reads as a real observation that the animal
had none of that cell type. Masking makes the treatment consistent: a failed
slide is absent from every panel.

Marker-**independent** endpoints are unaffected, because with a complete
phenotype partition the total compound-positive count is the sum of the two
compound-positive phenotypes however the marker split them. Separate the two in
any write-up.

### `qc/effective_coverage.csv`
What survives exclusions and marker failures, per region and group, with the
sexes represented. In a design with one animal per group per sex, losing one
animal from a region turns a level into a single observation of a single sex.

### `qc/background.csv`
Control animals carry no compound, so whatever the classifier calls positive in
them is the empirical false-positive rate. That number, not zero, is the floor
against which treated values are read, and it is printed on the relevant figure
captions.

### `qc/annotation_recodes.csv`
Every annotation rename applied, with its reason and scope. Renaming an
anatomical annotation asserts that the tissue is not what the annotation said it
was; two different tissues can be entirely different exposure compartments. A
recode without a recorded reason is rejected at load time, and a recode is
always scoped to specific animals so that a future sample with a correct label
is not silently caught.
