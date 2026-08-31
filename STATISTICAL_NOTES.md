# Statistical notes

Background for the model sections. Written so a reviewer can see why the
pipelines do what they do without reading the code.

---

## 1. What one row is

The dominant error in this literature is not the choice of test. It is that
analyses do not say what one row of the analysed data frame represents.

Cells are nested in sections, sections in animals. **Cells from one animal are
not independent observations of the treatment; the animal is.** Nested designs
are common across the field, and ignoring the nesting inflates type I error far
beyond the nominal level. Using individual as a batch-effect covariate does not
fix it.

Note the tension with the classic protocol advice that a small number of animals
with a few sections each suffices. That is a statement about precision *within*
an animal, and it is exactly the reasoning that produces pseudoreplication.

**Power comes from animals, not cells.**

`design_power_note()` prints the smallest achievable two-sided Wilcoxon p:

| Group sizes | Minimum p |
|---|---|
| 3 v 3 | 0.100 |
| 4 v 4 | 0.029 |
| 4 v 5 | 0.016 |
| 6 v 6 | 0.002 |

At 3 v 3, no amount of separation reaches 0.05. Worth knowing before anyone asks
for a p-value.

---

## 2. Distribution

Transcription is bursty, so steady-state counts are **negative binomial, not
Poisson**. Four further sources of overdispersion stack on top: cell
heterogeneity, section-plane truncation, and detection-efficiency variation
across the slide.

Assume NB and test it. `compare_count_families()` fits Poisson, NB1 and NB2 and
compares by AIC.

**Never `log(x+1)` count data.** Use a log *link* on untransformed counts.

---

## 3. Two kinds of zero

- **Structural** — the cell does not express the gene.
- **Sampling** — it does, but this plane missed it.

A single NB cannot separate them. A hurdle can, and its two parts are exactly
the two quantities usually reported by hand:

| Model part | Reported by hand as |
|---|---|
| Logistic on P(count > 0) | "% positive cells" |
| Zero-truncated NB on the rest | "dots per positive cell" |

Fitting both says **which component moved**. A treatment that recruits more
cells into expression is different biology from one that raises output per
expressing cell, and a t-test on mean dots per cell cannot tell them apart.

```r
library(glmmTMB)

m <- glmmTMB(
  dots ~ group + (1 | animal/section),
  family    = nbinom2,        # compare vs nbinom1 and poisson by AIC
  ziformula = ~ 1,            # or ~ group for a group-dependent hurdle
  data      = cells)

# unequal ROI areas: use an offset, never a pre-computed ratio
m_area <- glmmTMB(dots ~ group + offset(log(area_mm2)) + (1 | animal),
                  family = nbinom2, data = rois)

performance::icc(m)   # report this
```

Python equivalents: `statsmodels` GEE with `groups="animal"` and a
NegativeBinomial family for a population-average effect; `bambi` or PyMC for a
Bayesian NB GLMM with nesting.

---

## 4. Rules the pipelines enforce

1. Random effect for animal, and section within animal. Always.
2. `offset(log(area))` or `offset(log(n_cells))`, never a pre-computed ratio. A
   ratio forces the exposure coefficient to exactly 1 and discards the precision
   information.
3. Report the ICC. At 0.3, 500 cells per animal carry roughly the information of
   three independent observations. This is the number that makes the point land.
4. Never `log(x+1)` count data.
5. Sanity-check against pseudobulk: aggregate to one value per animal and
   re-test. Disagreement in **direction** means the model is misspecified;
   merely more power is expected.

---

## 5. Multiplicity

No systematic practice exists in this literature. Apply BH-FDR across the full
probe by region grid, **declared before analysis**.

Regions within an animal are correlated. If region structure is the question,
model it — region as a fixed effect with animal as a random effect, then
contrasts — rather than fitting separate models per region.

Correcting within a figure panel but not across panels is not a correction.

---

## 6. Ratios and their nulls

Expressing a treated value as a fraction of a control mean treats that mean as
known exactly. It is not, and in this assay it is often the noisiest quantity in
the analysis: control coefficients of variation of 30–60% between animals are
routine.

Two consequences.

Get the ratio from a model with a log link and exponentiate the coefficient, so
the interval includes control variance.

Show what the same normalisation does to the **control animals themselves**.
`control_null_distribution()` uses a leave-one-out reference so each control
animal is compared against the others rather than partly against itself.

**Never clip.** A one-sided clip censors animals on one side of the reference and
not the other, so the mean is biased and can never be negative. The failure mode
is specific: apply a clipped estimator to control animals compared against their
own group mean and it will return a substantial apparent effect, because every
animal above the mean is forced to zero while every animal below it contributes
its full value. `null_effect_summary()` quantifies exactly that, per study.

---

## 7. Degenerate summary statistics

When most cells are at zero the per-animal median is pinned to a small integer,
and a knockdown computed from it then takes only a handful of values. In the
worst case it alternates between 0 and 100 per animal: a coin flip presented as
a measurement.

`check_degenerate_summaries()` counts distinct values per group and flags
anything taking two or fewer across animals. Use `median_among_pos`, the median
among positive cells, which is the hurdle's second component and does not
collapse.

---

## 8. Where the literature is thin

Worth knowing, and worth being the group that fills it.

- **Negative binomial, zero-inflated and hurdle models applied to puncta
  counts** appear to be essentially absent from the published RNAscope
  literature.
- **No dose–response (Emax, 4PL) or PK/PD model** has been fitted to a spatial
  ISH readout. That is directly relevant to tissue-exposure questions in
  oligonucleotide programmes, and the `dose_response` mode in pipeline 02 is a
  first step toward it.
- **No published head-to-head comparison** of the major image-analysis platforms
  for this assay. Choose on workflow and validation burden, not on accuracy
  claims.
- **No reporting standard** for quantitative ISH. The nearest usable substitutes
  are the Schmied 2024 image-analysis checklist, QUAREP-LiMi and ARRIVE 2.0. See
  `METHODS_CHECKLIST.md`.
