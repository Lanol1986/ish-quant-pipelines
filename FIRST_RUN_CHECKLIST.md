# First-run checklist

The code has not been executed against real data. Work through this once per
pipeline before relying on the output.

## Before running

1. **Column names.** Open one export and compare its headers against the config.
   Probe names appear verbatim inside the column headers of a subcellular
   export, and object-export header text varies between software versions.
   Adjust `assay$probes[[1]]$export_name` or `export_columns` until they match
   exactly, including any double spaces. The readers fail with the full list of
   available columns, so the first failure tells you what to write.

2. **`paths$id_regex`.** Check the first capture group extracts the animal id
   from your filenames. `build_file_manifest()` stops with the offending
   filenames if not.

3. **Roster.** `design$animals` must list every animal, and every animal must
   either have a file or be listed in `design$expected_absent` with a reason.

4. **Regions.** Run once and read `qc/unresolved_regions.csv` before trusting
   anything. Add entries to `design$region_aliases` where annotation names do
   not match the declared region names.

5. **Delimiter.** `paths$delimiter` if your export is not the pipeline default
   (tab for subcellular, comma for object and summary).

## After the first run

6. **`qc/read_qc.csv`** — `n_na_est` should be small. More than half the values
   in a file being NA and set to zero means the export settings were wrong.

7. **`qc/expression_regime.csv`** — is a per-cell count defensible? If `regime`
   is `dense_or_saturated`, switch the primary readout for that stratum. If
   `is_integer_valued` is FALSE, record the divisor in
   `readout$single_dot_calibration`.

8. **`qc/threshold.yml`** — compare the threshold derived from control material
   against the config value. If they differ substantially, decide and record.

9. **`qc/cluster_rate.csv`** — anything flagged `CHECK STAINING` developed
   differently and is not a staining replicate.

10. **`qc/degenerate_check.csv`** — anything flagged `DEGENERATE` must not be an
    endpoint.

11. **`tables/null_effect_summary.csv`** — the apparent effect the pipeline
    produces from control animals alone. Any treated effect smaller than this is
    not established.

12. **`run_manifest/run_manifest.yml`** — read `unrecorded_parameters`, fill them
    in the config, rerun.

## Things that will probably need adjusting

- `read_subcellular_export()` builds optional column names from the probe name.
  If your channels are named differently, pass them through `extra_cols`.
- `qc_staining_drift()` guesses metric columns by pattern (`_od_mean$`,
  `_intensity$`, `hematoxylin`). Pass `metric_cols` explicitly if the guess
  misses.
- Pipeline 03 derives the intensity column for the second marker by name
  matching. Check `qc/dual_intensity.csv` is non-empty; if not, rename the key
  in `export_columns$marker_intensities` so it contains the marker name.
- Pipeline 04 summary mode expects the four phenotype columns in the declared
  order: marker+/compound+, marker+/compound-, marker-/compound+,
  marker-/compound-. Getting the order wrong swaps the compartments silently, so
  check `qc/phenotype_partition.csv` and the marker totals against the export.
- `glmmTMB` convergence on very large cell-level tables can be slow. Consider
  subsampling cells per animal for model fitting — the animal-level summaries
  are unaffected — and record that in the methods.

## Sanity tests worth running once

- Fit the same comparison two ways and check they agree in direction: the
  hurdle GLMM in `models/`, and the animal-level Wilcoxon in
  `qc/pseudobulk_check.csv`. Disagreement in direction means the model is
  misspecified. The model having more power is expected.
- In pipeline 03, read `qc/dual_rate.csv`. If `n_dual` is zero everywhere, the
  classifier was set to mutually exclusive upstream and resolved the duals by a
  rule you did not choose. That is worth knowing either way.
- In pipeline 04, confirm `qc/phenotype_partition.csv` shows
  `partition_ok = TRUE` for every row before reading any percentage.
