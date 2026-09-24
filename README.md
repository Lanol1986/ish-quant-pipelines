# RNAscope / ISH quantification toolkit

Reusable R pipelines for quantifying RNAscope and related in situ hybridisation
data for oligonucleotide therapeutic programmes.

Four pipelines, one self-contained R script each. Everything study-specific
lives in a YAML config, so a new study needs a new config and no change to any
script.

---

## The four pipelines

| # | Script | Assay | Endpoint |
|---|---|---|---|
| 01 | `01_biodistribution_chromogenic.R` | Probe against the compound, chromogenic, brightfield | Compound **distribution** |
| 02 | `02_target_kd_chromogenic.R` | Probe against the target transcript, chromogenic | **Knockdown** vs control, or dose **trend** |
| 03 | `03_multiplex_celltype.R` | Multiplex compound + target + antibodies, object level | Cell-type-resolved knockdown and uptake |
| 04 | `04_compound_distribution_ihc.R` | Compound probe + one antibody, four-way phenotype | Distribution by compartment |
| 05 | `05_multiplex_two_panel_integration.R` | Two or more independently analysed HALO AI multiplex panels | Panel-aware integration and concordance |

Pipeline 02 has two modes, set by `design$type`:

- `two_group` / `multi_group` — group contrasts with the animal as the unit.
- `dose_response` — where there is one animal per dose, no group contrast is
  possible and none is attempted. A **trend across dose** is estimable and is
  fitted from every animal.

Pipeline 04 has two input modes, set by `detection$export_level`:

- `summary` — aggregated exports, one row per image and analysis region. Prints
  at startup what it cannot do rather than approximating it.
- `object` — one row per cell. Enables the hurdle model and the subcellular
  pattern statistic.

---

## Install

R 4.1 or later.

```r
install.packages(c(
  "tidyverse", "fs", "yaml", "scales", "digest",     # required
  "glmmTMB", "emmeans", "performance",               # modelling
  "FNN", "ragg"                                      # QC and graphics
))
```

Modelling packages are optional. Without them the descriptive and QC parts of
every pipeline still run and report what was skipped.

`ragg` is worth having on a cluster. Without a UTF-8 locale, R silently replaces
en-dashes and `>=` in plot labels with dots and only warns through `mbcsToSbcs`.

---

## Run

```bash
Rscript 02_target_kd_chromogenic.R config/my_study.yml
```

From a Jupyter notebook with the R kernel, see `notebooks/`:

```r
CONFIG_PATH <- "config/my_study.yml"
source("02_target_kd_chromogenic.R")
```

The notebooks are thin drivers. Keeping the logic in `.R` files means git diffs
are readable and the same code runs unattended.

---

## Output

```
results/
├── tables/          animal-level summaries, relative-to-control, null distribution
├── plots/           figures, one point per animal
├── qc/              every QC table; read these first
├── models/          fitted models, coefficient tables, ICC
├── run_manifest/    config copy + hash, input file hashes, sessionInfo, git commit
├── METHODS_STUB.md  methods section populated from the config
└── EVIDENCE_BOUNDARIES.txt   (compound-detection pipelines)
```

The run manifest lists every config field left blank under
`unrecorded_parameters`. That list is itself a finding: the commonest
reproducibility failure in this field is that the decisive parameter tuning
happened inside a GUI and was never written down.

---

## What the pipelines enforce

**One row of the analysed data frame is one animal.** Cells are nested in
sections within animals. Cells from one animal are not independent observations
of the treatment; the animal is. Treating them as independent inflates type I
error severely, and using animal as a batch covariate does not fix it.
Cell-level tables exist here only as QC.

**No boxplot over fewer than four animals.** A box across two points invents
quartiles. `plot_animal_points()` falls back to points and a mean bar and says
so in the caption.

**Nothing is clipped.** A one-sided clip such as `pmax(100 - remaining, 0)`
censors every animal above the control mean to exactly zero while animals below
it contribute fully. The mean of that is biased upward whenever there is
between-animal spread and can never be negative, so it reports an effect even
when the groups are identical. A negative knockdown is information.

**Every ratio comes with its null.** `control_null_distribution()` applies the
same normalisation to the control animals with a leave-one-out reference, and
`null_effect_summary()` reports the apparent effect that produces. Run it before
quoting any effect size: if the treated effect sits inside the
control-versus-control spread, there is no effect to report.

**Counts are counts, or they are called an index.** An "estimated spots" or
"copies" column derived from cluster area or intensity is a model output built
on an assumed single-dot value. `qc_expression_regime()` reports which regime
each stratum is in, and the modelling functions refuse non-integer input rather
than fitting a negative binomial to a continuous estimate.

**Degenerate statistics are flagged.** When most cells are at zero, a per-animal
median is pinned to a small integer and a ratio built on it can take only a
handful of values. `check_degenerate_summaries()` catches that before it becomes
a figure.

**Thresholds come from the negative control**, as mean plus one SD of non-zero
puncta. `qc_threshold_sweep()` then shows whether the result survives a
different choice.

**Absences are declared.** `assert_roster()` stops if an animal in the config
has no file, or a file has no config entry. An exclusion goes in
`design$expected_absent` with an image-based reason recorded before the numbers
were seen.

**Missing measurements are not zeros.** Pipeline 04 flags a marker count of zero
against its stratum median and masks every marker-derived quantity there.
Otherwise a percentage divides by zero and vanishes from its panel with a
"Removed n rows" warning, while a density divides by area, comes out as a clean
zero, and is plotted on the axis where it reads as a real observation.

**Never divide by a housekeeping probe.** Degradation hits high-expressing
housekeepers hardest and low-to-moderate expressors least, so the control's
response to a pre-analytical insult is not proportional to the target's and the
ratio adds noise instead of cancelling artefact. Use one as a **sample-level
inclusion gate** (`housekeeping_inclusion_gate()`): if a sample fails, drop it;
do not rescue it by dividing.

---

## Statistical approach

Cells are nested in sections within animals; the animal remains the experimental
unit. The count workflow uses a **true two-part hurdle analysis** when discrete
integer puncta are scientifically defensible:

- **Part 1:** binomial mixed model for detectable signal.
- **Part 2:** zero-truncated Poisson/NB mixed model among positive cells.

This is intentionally not implemented with `glmmTMB(..., ziformula=...)`.
A zero-inflated mixture and a hurdle model answer different questions. The
pipeline exports the two components separately so distribution breadth and
conditional burden cannot be conflated.

Continuous or cluster-derived burden indices are not silently rounded into
counts for model fitting. Use an appropriate continuous positive-burden model
or a genuine integer single-spot column.

Multiplicity is controlled over a predeclared grid, and cell-level inference
never substitutes for animal replication.

### Platform split

- **QuPath / brightfield:** chromogenic RED RNAscope and DIG-DAB.
- **HALO AI / fluorescence:** duplex/multiplex smRNA/mRNA RNAscope + antibodies.

The platform-specific import/QC layers stay separate. Harmonisation starts at
validated animal/region/cell-type summaries.

### Nested brain anatomy

Pipeline 03 supports an optional `subregion` level. For cerebellum this can
retain the parent `Cerebellum` result while additionally reporting molecular,
Purkinje-cell, granular and white-matter layers. Layers are repeated anatomical
measurements within an animal, never additional biological replicates.

### Two-panel integration

Run each fluorescent panel independently through Pipeline 03 first. Pipeline 05
then performs panel-aware integration at the animal level, checks shared-cell
types for concordance, and normalises target mRNA to matched controls within
panel before optional combination. Raw fluorescence intensity remains
panel-specific unless cross-panel calibration is explicitly documented.

---

## Interpretive limits for oligonucleotide work

Written to `EVIDENCE_BOUNDARIES.txt` by pipelines 01, 03 and 04. In short:

- **Detected compound is not active compound.** Two uptake pathways exist, one
  ending in the endolysosome where the oligonucleotide is inert. Knockdown does
  not correlate with bulk intracellular accumulation; measured endosomal escape
  is 1–2%, while fewer than 2,000 cytosolic copies suffice to silence. Report
  presence and pattern. Perinuclear versus peripheral distribution has
  distinguished productive from non-productive uptake where total amount did not.
- **No limit of detection has been published** for oligonucleotide ISH in any
  units. A negative region means "below an unknown threshold", not "no compound".
- **No correlation against LC-MS or hybridisation-ELISA** tissue concentration
  has been published. Do not phrase the readout as a concentration.
- **Full-length compound versus n-1 or chain-shortened metabolites** is untested
  by these assays. Only LC-MS/MS answers that.
- For siRNA, **state which strand the probe targets**. Passenger-strand
  persistence is a weaker activity proxy than guide, since the passenger is
  discarded on RISC loading.

---

## Layout

```
pipelines/    one self-contained script per assay type
config/       study_config_TEMPLATE.yml + one worked example per pipeline
notebooks/    thin Jupyter drivers
docs/         methods checklist, QC reference, statistical notes, first-run guide
```

Each pipeline script has the same eleven sections: environment, utilities,
configuration, input, quality control, aggregation, relative to control, models,
figures, provenance, main. The shared sections are identical across the four
files by design: each script is portable on its own, and a change to shared
logic is a deliberate, visible edit in each file rather than an invisible
side-effect of editing a library.

---

## Status

**Not yet executed against real data.** The code is written against the column
names and value ranges these exports normally carry, but it has not been run end
to end. Do a dry run on one study before relying on it, and expect to adjust the
column mappings in the config first. `docs/FIRST_RUN_CHECKLIST.md` lists what to
check and in what order.
