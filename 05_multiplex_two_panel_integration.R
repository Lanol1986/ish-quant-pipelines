## =============================================================================
## PIPELINE 05 -- TWO-PANEL MULTIPLEX INTEGRATION
##
## Integrates animal-level outputs from two or more independently analysed
## HALO AI fluorescent RNAscope/IF panels. Panels may come from adjacent sections
## and may carry different antibody sets.
##
## IMPORTANT:
##   - panels are repeated assays from the same animals, not additional n;
##   - each panel must be QC'd and summarised by 03_multiplex_celltype.R first;
##   - target mRNA is normalised to its matched control WITHIN PANEL before any
##     shared-cell-type integration;
##   - raw fluorescence intensity is never averaged across panels unless the
##     configuration explicitly states that acquisition/amplification/unmixing
##     were calibrated across panels.
##
## Usage:
##   Rscript 05_multiplex_two_panel_integration.R config.yml
## =============================================================================

PKGS_REQUIRED <- c("tidyverse", "yaml", "fs", "scales")

check_packages <- function() {
  miss <- PKGS_REQUIRED[
    !vapply(PKGS_REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing required packages: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(tidyverse); library(yaml); library(fs); library(scales)
  })
}

`%||%` <- function(a, b) {
  if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

read_config <- function(path) {
  if (!file.exists(path)) stop("Config not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  need <- c("study", "panels", "integration", "output")
  miss <- setdiff(need, names(cfg))
  if (length(miss)) {
    stop("Config missing sections: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  if (length(cfg$panels) < 2) stop("At least two panels are required.", call. = FALSE)
  cfg
}

read_panel_summary <- function(panel) {
  path <- path.expand(panel$summary_csv)
  if (!file.exists(path)) stop("Panel summary not found: ", path, call. = FALSE)

  d <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("group", "animal_id", "region", "celltype")
  miss <- setdiff(required, names(d))
  if (length(miss)) {
    stop(panel$name, " is missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  d |>
    dplyr::mutate(
      panel = panel$name,
      animal_id = as.character(animal_id),
      group = as.character(group),
      region = as.character(region),
      celltype = as.character(celltype))
}

panel_roster_qc <- function(all_panels) {
  all_panels |>
    dplyr::distinct(panel, animal_id, group) |>
    tidyr::complete(panel, animal_id, fill = list(group = NA_character_)) |>
    dplyr::mutate(status = ifelse(is.na(group), "MISSING", "PRESENT"))
}

normalise_within_panel <- function(d, value_col, control_group) {
  if (!value_col %in% names(d)) {
    stop("Target metric '", value_col, "' is absent.", call. = FALSE)
  }

  refs <- d |>
    dplyr::filter(group == control_group) |>
    dplyr::group_by(panel, region, celltype) |>
    dplyr::summarise(
      control_mean = safe_mean(.data[[value_col]]),
      n_control = sum(is.finite(.data[[value_col]])),
      .groups = "drop")

  d |>
    dplyr::left_join(refs, by = c("panel", "region", "celltype")) |>
    dplyr::mutate(
      ratio_to_control = .data[[value_col]] / control_mean,
      pct_remaining = 100 * ratio_to_control,
      pct_reduction = 100 * (1 - ratio_to_control),
      normalised_metric = value_col)
}

shared_anchor_table <- function(norm, shared_celltypes) {
  norm |>
    dplyr::filter(celltype %in% shared_celltypes) |>
    dplyr::select(panel, group, animal_id, region, celltype,
                  pct_remaining, pct_reduction,
                  dplyr::any_of(c("pct_compound_pos",
                                  "median_compound_among_pos",
                                  "mean_compound")))
}

pair_shared_panels <- function(anchor) {
  panels <- unique(anchor$panel)
  if (length(panels) != 2) {
    return(tibble::tibble())
  }

  a <- anchor |>
    dplyr::filter(panel == panels[1]) |>
    dplyr::select(-panel) |>
    dplyr::rename(
      pct_remaining_panel1 = pct_remaining,
      pct_reduction_panel1 = pct_reduction)

  b <- anchor |>
    dplyr::filter(panel == panels[2]) |>
    dplyr::select(-panel) |>
    dplyr::rename(
      pct_remaining_panel2 = pct_remaining,
      pct_reduction_panel2 = pct_reduction)

  by <- intersect(c("group", "animal_id", "region", "celltype"),
                  intersect(names(a), names(b)))

  dplyr::inner_join(a, b, by = by) |>
    dplyr::mutate(
      panel_difference = pct_remaining_panel1 - pct_remaining_panel2,
      panel_mean = (pct_remaining_panel1 + pct_remaining_panel2) / 2)
}

panel_concordance_summary <- function(paired) {
  if (!nrow(paired)) return(tibble::tibble())

  paired |>
    dplyr::group_by(region, celltype) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      pearson_r = if (dplyr::n() >= 3)
        stats::cor(pct_remaining_panel1, pct_remaining_panel2,
                   use = "complete.obs", method = "pearson") else NA_real_,
      spearman_rho = if (dplyr::n() >= 3)
        stats::cor(pct_remaining_panel1, pct_remaining_panel2,
                   use = "complete.obs", method = "spearman") else NA_real_,
      mean_difference = safe_mean(panel_difference),
      sd_difference = if (sum(is.finite(panel_difference)) >= 2)
        stats::sd(panel_difference, na.rm = TRUE) else NA_real_,
      loa_low = mean_difference - 1.96 * sd_difference,
      loa_high = mean_difference + 1.96 * sd_difference,
      .groups = "drop")
}

combine_shared_target <- function(anchor, method = "mean") {
  if (!method %in% c("mean", "median", "keep_panel")) {
    stop("integration$shared_target_method must be mean, median or keep_panel.",
         call. = FALSE)
  }
  if (method == "keep_panel") return(anchor)

  fun <- if (method == "mean") safe_mean else function(x)
    stats::median(x[is.finite(x)], na.rm = TRUE)

  anchor |>
    dplyr::group_by(group, animal_id, region, celltype) |>
    dplyr::summarise(
      n_panels = dplyr::n_distinct(panel),
      pct_remaining_integrated = fun(pct_remaining),
      pct_reduction_integrated = fun(pct_reduction),
      .groups = "drop")
}

plot_panel_concordance <- function(paired) {
  ggplot2::ggplot(
    paired,
    ggplot2::aes(x = pct_remaining_panel1,
                 y = pct_remaining_panel2,
                 shape = group)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::facet_grid(celltype ~ region, scales = "free") +
    ggplot2::labs(
      x = "Panel 1 target remaining (%)",
      y = "Panel 2 target remaining (%)",
      title = "Shared-cell-type agreement across panels",
      caption = "Same animals measured in adjacent/independent panels; points are paired, not additional biological replicates.") +
    ggplot2::theme_classic(base_size = 12)
}

plot_bland_altman <- function(paired) {
  ggplot2::ggplot(
    paired,
    ggplot2::aes(x = panel_mean, y = panel_difference, shape = group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::facet_grid(celltype ~ region, scales = "free") +
    ggplot2::labs(
      x = "Mean target remaining across panels (%)",
      y = "Panel 1 - Panel 2 (%)",
      title = "Bland-Altman view of panel agreement",
      caption = "Use this as a comparability check before combining shared-cell-type target readouts.") +
    ggplot2::theme_classic(base_size = 12)
}

main <- function(config_path) {
  check_packages()
  cfg <- read_config(config_path)

  out <- path.expand(cfg$output$dir)
  dirs <- list(root = out,
               tables = file.path(out, "tables"),
               qc = file.path(out, "qc"),
               plots = file.path(out, "plots"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  all <- purrr::map_dfr(cfg$panels, read_panel_summary)

  qc_roster <- panel_roster_qc(all)
  readr::write_csv(qc_roster, file.path(dirs$qc, "panel_roster.csv"))

  target_metric <- cfg$integration$target_metric %||% "mean_all_cells"
  ctrl <- cfg$integration$control_group
  norm <- normalise_within_panel(all, target_metric, ctrl)
  readr::write_csv(norm, file.path(dirs$tables, "all_panels_normalised.csv"))

  shared <- unlist(cfg$integration$shared_celltypes)
  anchor <- shared_anchor_table(norm, shared)
  readr::write_csv(anchor, file.path(dirs$tables, "shared_celltype_panel_values.csv"))

  paired <- pair_shared_panels(anchor)
  concordance <- panel_concordance_summary(paired)
  readr::write_csv(paired, file.path(dirs$qc, "shared_celltype_paired_values.csv"))
  readr::write_csv(concordance, file.path(dirs$qc, "panel_concordance_summary.csv"))

  if (nrow(paired)) {
    ggplot2::ggsave(
      file.path(dirs$plots, "QC_panel_concordance.png"),
      plot_panel_concordance(paired),
      width = 11, height = 7, dpi = 300, bg = "white")
    ggplot2::ggsave(
      file.path(dirs$plots, "QC_panel_bland_altman.png"),
      plot_bland_altman(paired),
      width = 11, height = 7, dpi = 300, bg = "white")
  }

  integrated <- combine_shared_target(
    anchor, cfg$integration$shared_target_method %||% "keep_panel")
  readr::write_csv(
    integrated, file.path(dirs$tables, "shared_celltype_integrated_target.csv"))

  # Compound distribution remains panel-specific by default. A panel may carry
  # different amplification, exposure, unmixing and antibody/segmentation
  # context. Only target metrics normalised to matched controls are integrated
  # here unless the study explicitly validates cross-panel comparability.
  if (isTRUE(cfg$integration$combine_raw_intensity)) {
    if (!isTRUE(cfg$integration$intensity_calibrated_across_panels)) {
      stop("combine_raw_intensity=TRUE but intensity_calibrated_across_panels is FALSE.",
           call. = FALSE)
    }
    warning("Raw-intensity integration was explicitly enabled. Document the ",
            "cross-panel calibration in the methods.", call. = FALSE)
  }

  message("Two-panel integration complete: ", dirs$root)
  invisible(list(all_panels = all, normalised = norm, paired = paired,
                 concordance = concordance, integrated = integrated))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    stop("Usage: Rscript 05_multiplex_two_panel_integration.R <config.yml>",
         call. = FALSE)
  }
  RESULT <- main(args[1])
}
