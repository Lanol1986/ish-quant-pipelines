## =============================================================================
## PIPELINE 04 -- COMPOUND DISTRIBUTION WITH ONE ANTIBODY
##
## Probe against the compound plus a single antibody, scored as a four-way
## phenotype: marker+/compound+, marker+/compound-, marker-/compound+,
## marker-/compound-.
##
## Endpoint: DISTRIBUTION and PATTERN. Not concentration, not activity.
##
## Two input modes, set by detection$export_level:
##   summary  aggregated, one row per image and analysis region
##   object   one row per cell
##
## SELF-CONTAINED. This file has no dependencies on any other file in the
## repository. Everything study-specific lives in a YAML config.
##
## Usage
##   Rscript pipelines/04_compound_distribution_ihc.R config/my_study.yml
##
## or from a notebook with the R kernel:
##   CONFIG_PATH <- "config/my_study.yml"
##   source("pipelines/04_compound_distribution_ihc.R")
##
## Contents
##   1  Environment            7  Relative to control
##   2  Utilities              8  Models
##   3  Configuration          9  Figures
##   4  Input                 10  Provenance
##   5  Quality control       11  Main
##   6  Aggregation
## =============================================================================

## =============================================================================
## SECTION 1 -- ENVIRONMENT
## =============================================================================

PKGS_REQUIRED <- c("tidyverse", "fs", "yaml", "scales", "digest")
PKGS_OPTIONAL <- c("glmmTMB", "emmeans", "performance", "ragg")

check_packages <- function(install = FALSE) {
  missing <- PKGS_REQUIRED[
    !vapply(PKGS_REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    if (install) utils::install.packages(missing)
    else stop("Missing required packages: ", paste(missing, collapse = ", "),
              call. = FALSE)
  }
  absent <- PKGS_OPTIONAL[
    !vapply(PKGS_OPTIONAL, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    message("Optional packages absent; steps needing them will be skipped: ",
            paste(absent, collapse = ", "))
  }
  invisible(absent)
}

## A non-UTF-8 LC_CTYPE makes R silently replace characters such as en-dashes
## and the >= sign in plot labels with a dot, and only warns via mbcsToSbcs.
set_utf8_locale <- function() {
  for (loc in c("en_US.UTF-8", "C.UTF-8", "en_GB.UTF-8")) {
    res <- suppressWarnings(try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE))
    if (!inherits(res, "try-error") && nzchar(res)) return(invisible(TRUE))
  }
  warning("No UTF-8 locale available. Keep plot labels ASCII-only, or install ",
          "ragg, which this pipeline uses automatically when present.",
          call. = FALSE)
  invisible(FALSE)
}

setup_session <- function(install = FALSE) {
  check_packages(install)
  set_utf8_locale()
  suppressPackageStartupMessages({
    library(tidyverse); library(fs); library(yaml); library(scales)
  })
  options(dplyr.summarise.inform = FALSE, stringsAsFactors = FALSE, warn = 1)
  invisible(TRUE)
}

## =============================================================================
## SECTION 2 -- UTILITIES
## =============================================================================

`%||%` <- function(a, b) {
  if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
}

safe_mean   <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_median <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else stats::median(x) }
safe_sd     <- function(x) { x <- x[is.finite(x)]; if (length(x) <= 1) NA_real_ else stats::sd(x) }
safe_pct    <- function(x) { x <- x[!is.na(x)];    if (!length(x)) NA_real_ else mean(x) * 100 }

#' Column headers to safe snake_case, handling the micro sign sensibly.
make_safe_name <- function(x) {
  x |>
    gsub("\u00b5m\u00b2", "um2", x = _, fixed = TRUE) |>
    gsub("\u03bcm\u00b2", "um2", x = _, fixed = TRUE) |>
    gsub("\u00b5m", "um", x = _, fixed = TRUE) |>
    gsub("\u03bcm", "um", x = _, fixed = TRUE) |>
    gsub("[^A-Za-z0-9]+", "_", x = _) |>
    gsub("_+", "_", x = _) |>
    gsub("^_|_$", "", x = _) |>
    tolower()
}

#' Map a free-text annotation name onto one of the declared regions.
#'
#' Annotation names pick up sample suffixes and inconsistent capitalisation.
#' Anything that matches no declared region returns NA so it can be COUNTED
#' rather than quietly dropped.
#'
#' @param x raw annotation names
#' @param regions declared region names from the config
#' @param aliases optional named list, e.g. list(Cortex = c("ctx", "cortical"))
resolve_region <- function(x, regions, aliases = NULL) {
  if (is.null(regions) || !length(regions)) return(as.character(x))

  norm <- function(s) {
    s <- tolower(as.character(s))
    s <- gsub("_[0-9]+\\.[0-9]+$", "", s)   # trailing sample id
    s <- gsub("_[0-9]+$", "", s)            # trailing index
    s <- gsub("[^a-z0-9]+", "_", s)
    s <- gsub("_+", "_", s)
    gsub("^_|_$", "", s)
  }

  xn  <- norm(x)
  out <- rep(NA_character_, length(x))
  for (r in regions) {
    pats <- c(norm(r), if (!is.null(aliases[[r]])) norm(aliases[[r]]))
    for (p in pats) {
      hit <- is.na(out) & grepl(p, xn, fixed = TRUE)
      out[hit] <- r
    }
  }
  out
}

#' Interpret a positivity column of any exported type.
#'
#' Platforms write 0/1 in most exports and "Positive"/"Negative" in some.
#' Comparing with `%in% 1` works for the first and silently returns FALSE for
#' the second, which would turn every cell into "unclassified".
as_positive_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x > 0)
  if (is.character(x) || is.factor(x)) {
    xl <- tolower(trimws(as.character(x)))
    yes <- c("1", "true", "positive", "pos", "yes", "y")
    no  <- c("0", "false", "negative", "neg", "no", "n", "")
    unknown <- setdiff(unique(xl), c(yes, no))
    if (length(unknown)) {
      stop("Unrecognised values in a positivity column: ",
           paste(utils::head(unknown, 10), collapse = ", "), call. = FALSE)
    }
    return(xl %in% yes)
  }
  stop("Cannot interpret a positivity column of class ", class(x)[1], call. = FALSE)
}

## =============================================================================
## SECTION 3 -- CONFIGURATION
##
## Everything study-specific lives in a YAML file. A new study needs a new
## config and no change to this script.
## =============================================================================

CONFIG_REQUIRED_TOP <- c("study", "paths", "design", "assay", "readout", "output")

load_study_config <- function(path) {
  if (!file.exists(path)) stop("Config not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)

  missing_top <- setdiff(CONFIG_REQUIRED_TOP, names(cfg))
  if (length(missing_top)) {
    stop("Config is missing required sections: ",
         paste(missing_top, collapse = ", "),
         "\nCompare against config/study_config_TEMPLATE.yml", call. = FALSE)
  }

  cfg <- validate_design(cfg)
  cfg <- validate_readout(cfg)
  cfg <- validate_paths(cfg)

  cfg$.config_path <- normalizePath(path)
  cfg$.config_hash <- digest::digest(file = path, algo = "sha256")
  cfg
}

validate_design <- function(cfg) {
  d <- cfg$design
  if (is.null(d$animals) || !length(d$animals)) {
    stop("design$animals is required: one entry per animal.", call. = FALSE)
  }

  animals <- dplyr::bind_rows(lapply(d$animals, tibble::as_tibble))
  miss <- setdiff(c("animal_id", "group"), names(animals))
  if (length(miss)) {
    stop("design$animals entries need animal_id and group. Missing: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(animals$animal_id)) {
    stop("Duplicate animal_id in design$animals: ",
         paste(unique(animals$animal_id[duplicated(animals$animal_id)]),
               collapse = ", "), call. = FALSE)
  }

  if (is.null(d$group_levels)) {
    d$group_levels <- unique(animals$group)
    message("design$group_levels inferred as: ",
            paste(d$group_levels, collapse = ", "))
  }
  bad <- setdiff(unique(animals$group), d$group_levels)
  if (length(bad)) {
    stop("Groups present in design$animals but absent from design$group_levels: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }

  if (is.null(d$control_group)) {
    stop("design$control_group is required: the group other groups are ",
         "expressed relative to.", call. = FALSE)
  }
  if (!d$control_group %in% d$group_levels) {
    stop("design$control_group '", d$control_group,
         "' is not one of design$group_levels.", call. = FALSE)
  }

  if (is.null(d$type)) d$type <- "two_group"
  if (!d$type %in% c("two_group", "multi_group", "dose_response")) {
    stop("design$type must be two_group, multi_group or dose_response.",
         call. = FALSE)
  }
  if (identical(d$type, "dose_response") && !"dose" %in% names(animals)) {
    stop("design$type is dose_response but design$animals has no `dose` field.",
         call. = FALSE)
  }

  d$animals_tbl <- animals |>
    dplyr::mutate(animal_id = as.character(animal_id),
                  group     = factor(group, levels = d$group_levels))

  n_per_group <- dplyr::count(d$animals_tbl, group, name = "n")
  if (any(n_per_group$n < 2)) {
    warning("At least one group has a single animal. Group-level inference is ",
            "not possible; the pipeline runs descriptively and will not draw ",
            "boxplots or fit group contrasts. Groups affected: ",
            paste(n_per_group$group[n_per_group$n < 2], collapse = ", "),
            call. = FALSE)
  }

  cfg$design <- d
  cfg
}

validate_readout <- function(cfg) {
  r <- cfg$readout
  allowed <- c("pct_positive", "counts_per_cell", "burden_index",
               "area_fraction", "integrated_intensity")

  if (is.null(r$primary)) stop("readout$primary is required.", call. = FALSE)
  if (!r$primary %in% allowed) {
    stop("readout$primary must be one of: ", paste(allowed, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(r$positivity_threshold)) {
    stop("readout$positivity_threshold is required. Derive it from the ",
         "negative control rather than choosing it by eye.", call. = FALSE)
  }
  if (is.null(r$threshold_source)) {
    stop("readout$threshold_source is required: negative_control_probe, ",
         "control_animal, or fixed_declared.", call. = FALSE)
  }

  ## Any count derived from a cluster is an index built on an assumed
  ## single-dot size or intensity. The divisor has to travel with the results,
  ## or the numbers are not comparable between slides.
  if (is.null(r$single_dot_calibration)) {
    warning("readout$single_dot_calibration is not set. Cluster-derived counts ",
            "are indices, not counts; record the divisor used.", call. = FALSE)
  }
  if (is.null(r$countable_regime_max)) r$countable_regime_max <- 15

  cfg$readout <- r
  cfg
}

validate_paths <- function(cfg) {
  p <- cfg$paths
  if (is.null(p$input_dir))  stop("paths$input_dir is required.",  call. = FALSE)
  if (is.null(p$output_dir)) stop("paths$output_dir is required.", call. = FALSE)
  p$input_dir  <- path.expand(p$input_dir)
  p$output_dir <- path.expand(p$output_dir)
  if (!dir.exists(p$input_dir)) {
    stop("paths$input_dir does not exist: ", p$input_dir, call. = FALSE)
  }
  cfg$paths <- p
  cfg
}

setup_output_dirs <- function(cfg) {
  out  <- cfg$paths$output_dir
  dirs <- list(root   = out,
               tables = file.path(out, "tables"),
               plots  = file.path(out, "plots"),
               qc     = file.path(out, "qc"),
               models = file.path(out, "models"),
               run    = file.path(out, "run_manifest"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

## =============================================================================
## SECTION 4 -- INPUT
##
## Design rule: fail loudly.
##
## The failure this guards against is a detection export that is missing a
## column, or whose column is misparsed. Filling the gap with NA and then
## replacing NA with zero turns "this file was exported with the wrong settings"
## into "every cell in this animal has no signal", and that animal then looks
## like an immaculate control. Nothing here guesses.
## =============================================================================

#' Build a manifest of input files before reading any of them.
#'
#' @param input_dir directory holding the exports
#' @param pattern regex selecting the files
#' @param id_regex regex whose first capture group is the animal id
build_file_manifest <- function(input_dir, pattern, id_regex) {

  files <- fs::dir_ls(input_dir, regexp = pattern, recurse = FALSE, type = "file")
  ## Resource forks travel inside zip archives and look like real exports.
  files <- files[!stringr::str_detect(basename(files), "^\\._")]

  if (!length(files)) {
    stop("No input files matched '", pattern, "' in ", input_dir, call. = FALSE)
  }

  info <- file.info(as.character(files))
  manifest <- tibble::tibble(
    source_file = basename(files),
    path        = as.character(files),
    animal_id   = stringr::str_match(basename(files), id_regex)[, 2],
    size_mb     = round(info$size / 1024^2, 2),
    mtime       = info$mtime
  )

  if (any(is.na(manifest$animal_id))) {
    stop("Could not extract an animal id from: ",
         paste(manifest$source_file[is.na(manifest$animal_id)], collapse = ", "),
         "\nCheck paths$id_regex against the actual filenames.", call. = FALSE)
  }

  dplyr::arrange(manifest, animal_id, source_file)
}

#' Assert that the files found match the animals declared in the config.
#'
#' Silent absence of an animal is one of the easiest ways to change a result. An
#' animal that is legitimately absent belongs in design$expected_absent with a
#' reason, so the exclusion is recorded rather than implied.
assert_roster <- function(manifest, cfg) {
  expected <- cfg$design$animals_tbl$animal_id
  found    <- unique(manifest$animal_id)

  declared <- cfg$design$expected_absent %||% character(0)
  if (is.list(declared)) declared <- names(declared)

  missing    <- setdiff(expected, c(found, declared))
  unexpected <- setdiff(found, expected)

  if (length(missing)) {
    stop("Animals declared in the config have no input file: ",
         paste(missing, collapse = ", "),
         "\nEither add the files, or record them under design$expected_absent ",
         "with the reason they were excluded.", call. = FALSE)
  }
  if (length(unexpected)) {
    stop("Input files found for animals not in the config: ",
         paste(unexpected, collapse = ", "),
         "\nAdd them to design$animals or move them out of input_dir.",
         call. = FALSE)
  }
  if (length(declared)) {
    message("Declared exclusions honoured: ", paste(declared, collapse = ", "))
  }
  invisible(TRUE)
}

#' Read one delimited export with no silent repair.
#'
#' @param required_cols exact column names that must be present, case sensitive
#' @param delim tab for most whole-slide exports, comma for most object exports
read_export_strict <- function(path, required_cols, delim = "\t") {

  df <- readr::read_delim(path, delim = delim, show_col_types = FALSE,
                          progress = FALSE, name_repair = "minimal",
                          guess_max = 100000)

  ## A parsing problem makes a numeric column character; as.numeric() then makes
  ## it NA, and a downstream replace_na() makes it zero. Stop here instead.
  probs <- readr::problems(df)
  if (nrow(probs)) {
    print(utils::head(probs, 20))
    stop(nrow(probs), " parsing problem(s) in ", basename(path),
         ". Resolve them before proceeding.", call. = FALSE)
  }

  missing <- setdiff(required_cols, names(df))
  if (length(missing)) {
    stop("Missing required column(s) in ", basename(path), ": ",
         paste(missing, collapse = ", "),
         "\nAvailable columns:\n  ", paste(names(df), collapse = "\n  "),
         call. = FALSE)
  }
  df
}

#' Coerce to numeric and report, rather than hide, the NAs.
#'
#' In detection exports NA usually does mean zero. The count still belongs in
#' the QC table, so a file that is entirely NA is visible.
numeric_with_na_report <- function(df, col) {
  raw <- df[[col]]
  val <- suppressWarnings(as.numeric(raw))
  list(values = val,
       n_na = sum(is.na(val)),
       n_na_coerced = sum(is.na(val) & !is.na(raw)))
}

#' Read every file in a manifest and bind, carrying per-file QC.
#'
#' @param reader function(path) returning list(data, qc)
read_all <- function(manifest, reader) {
  results <- purrr::pmap(
    list(manifest$path, manifest$source_file, manifest$animal_id),
    function(path, source_file, animal_id) {
      res <- reader(path)
      res$data <- dplyr::mutate(res$data, animal_id = animal_id,
                                source_file = source_file)
      res$qc   <- dplyr::mutate(res$qc, animal_id = animal_id)
      res
    })
  list(data = dplyr::bind_rows(lapply(results, `[[`, "data")),
       qc   = dplyr::bind_rows(lapply(results, `[[`, "qc")))
}

#' Attach design metadata and make the nesting explicit.
attach_design <- function(tbl, cfg) {
  out <- dplyr::left_join(tbl, cfg$design$animals_tbl, by = "animal_id")

  if (any(is.na(out$group))) {
    stop("No design entry for animal(s): ",
         paste(unique(out$animal_id[is.na(out$group)]), collapse = ", "),
         call. = FALSE)
  }

  ## section_id is created even when there is one section per animal, so that
  ## the (1 | animal/section) term can be written without restructuring later.
  if (!"section_id" %in% names(out)) {
    out$section_id <- paste(out$animal_id, out$source_file, sep = ":")
  }
  out
}

#' Read an aggregated summary export (one row per image and analysis region).
#'
#' Summary exports carry no per-cell distribution, no coordinates and no way to
#' fit a hurdle model or check a classification against local density. Use them
#' when that is all there is, and read the matching object-level export when
#' possible.
#'
#' @param columns named list mapping short names to exact header text
#' @param phenotypes named list of mutually exclusive phenotype count columns
read_summary_export <- function(path, columns, phenotypes = NULL, delim = ",") {

  df <- read_export_strict(path,
    required_cols = unname(unlist(c(columns, phenotypes))), delim = delim)

  out <- tibble::tibble(.rows = nrow(df))
  for (nm in names(columns))    out[[nm]] <- df[[columns[[nm]]]]
  if (!is.null(phenotypes)) {
    for (nm in names(phenotypes)) out[[nm]] <- as.numeric(df[[phenotypes[[nm]]]])
  }

  list(data = out,
       qc = tibble::tibble(
         source_file = basename(path), n_rows = nrow(df),
         annotations = paste(unique(as.character(df[[columns$region]])),
                             collapse = " | ")))
}

#' Read an object-level export carrying marker positivity and probe counts.
#'
#' @param marker_cols named list, e.g. list(markera = "MarkerA Positive")
#' @param probe_cols named list, e.g. list(target = "Target_mRNA Copies")
read_object_export <- function(path, marker_cols, probe_cols,
                               extra_cols = character(0), delim = ",") {

  df <- read_export_strict(path,
    required_cols = unname(c(unlist(marker_cols), unlist(probe_cols))),
    delim = delim)

  out <- tibble::tibble(.rows = nrow(df))
  for (nm in names(marker_cols)) {
    out[[paste0(nm, "_pos")]] <- as_positive_flag(df[[marker_cols[[nm]]]])
  }
  for (nm in names(probe_cols)) {
    v <- numeric_with_na_report(df, probe_cols[[nm]])
    out[[nm]] <- tidyr::replace_na(v$values, 0)
  }

  optional <- intersect(unique(c(
    "Image Location", "Image Tag", "Analysis Region", "Object Id",
    "XMin", "XMax", "YMin", "YMax", extra_cols)), names(df))
  for (cc in optional) out[[make_safe_name(cc)]] <- df[[cc]]

  list(data = out,
       qc = tibble::tibble(source_file = basename(path), n_cells = nrow(df),
                           optional_cols = paste(optional, collapse = "; ")))
}

#' Apply annotation renames declared in the config, and log every one.
#'
#' Renaming an anatomical annotation asserts that the tissue is not what the
#' annotation said it was. It is a decision, not a cleaning step: it needs a
#' recorded reason and image verification, and it must be scoped to the animals
#' it applies to. A global rule would silently catch a future sample where the
#' original label was correct.
apply_annotation_recodes <- function(x, animal_id, recodes) {
  out <- as.character(x)
  log <- list()

  if (is.null(recodes) || !length(recodes)) {
    return(list(values = out, log = tibble::tibble()))
  }

  for (r in recodes) {
    if (is.null(r$from) || is.null(r$to) || is.null(r$reason)) {
      stop("Every entry in design$annotation_recodes needs `from`, `to` and ",
           "`reason`. A recode without a recorded reason is a silent ",
           "relabelling of one tissue as another.", call. = FALSE)
    }
    hit <- out == r$from
    if (!is.null(r$animal_id)) hit <- hit & (animal_id %in% r$animal_id)

    if (any(hit)) {
      log[[length(log) + 1]] <- tibble::tibble(
        from = r$from, to = r$to,
        animal_id = paste(unique(animal_id[hit]), collapse = ", "),
        n_rows = sum(hit), reason = r$reason)
      out[hit] <- r$to
    }
  }
  list(values = out,
       log = if (length(log)) dplyr::bind_rows(log) else tibble::tibble())
}

## =============================================================================
## SECTION 5 -- QUALITY CONTROL
##
## Every check here can change what the results mean. They run before any
## result is computed, and their output goes to qc/ so it is read first.
## =============================================================================

#' Classify the expression regime, which determines what readout is legitimate.
#'
#' Sparse, up to roughly 10-15 puncta per cell: a per-cell count is a genuine
#'   count.
#' Dense, puncta merging: cluster area divided by an assumed spot size, which is
#'   an INDEX and must be reported as one, with the divisor.
#' Saturated: integrated intensity or per cent labelled area.
qc_expression_regime <- function(cell_data, count_col,
                                 by = c("group", "animal_id", "region"),
                                 countable_max = 15) {
  by <- intersect(by, names(cell_data))

  cell_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_cells      = dplyr::n(),
      n_positive   = sum(.data[[count_col]] >= 1, na.rm = TRUE),
      pct_positive = 100 * mean(.data[[count_col]] >= 1, na.rm = TRUE),
      median_pos   = stats::median(
                       .data[[count_col]][.data[[count_col]] >= 1], na.rm = TRUE),
      p95          = stats::quantile(.data[[count_col]], 0.95, na.rm = TRUE),
      max_value    = suppressWarnings(max(.data[[count_col]], na.rm = TRUE)),
      n_above_max  = sum(.data[[count_col]] > countable_max, na.rm = TRUE),
      pct_pos_above_max = 100 * sum(.data[[count_col]] > countable_max, na.rm = TRUE) /
                          pmax(sum(.data[[count_col]] >= 1, na.rm = TRUE), 1),
      is_integer_valued = all(
        abs(.data[[count_col]] - round(.data[[count_col]])) < 1e-8, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(
      regime = dplyr::case_when(
        pct_pos_above_max > 25 ~ "dense_or_saturated",
        pct_pos_above_max > 5  ~ "mixed",
        TRUE                   ~ "countable"),
      readout_note = dplyr::case_when(
        regime == "countable" ~ "Counts per cell are defensible.",
        regime == "mixed" ~ paste0(
          "A minority of positive cells exceed the countable range. Use ",
          "% positive as primary and report counts as an index."),
        TRUE ~ paste0(
          "Values here are cluster-derived, not counts. Use % labelled area or ",
          "integrated intensity as primary.")))
}

#' Derive a positivity threshold from negative-control material.
#'
#' The rule is mean plus one SD of puncta on negative-control tissue, counting
#' only non-zero values. A threshold chosen because it is a round number is a
#' preference, not a threshold.
qc_threshold_from_control <- function(control_cells, count_col, n_sd = 1) {
  x  <- control_cells[[count_col]]
  x  <- x[is.finite(x)]
  nz <- x[x > 0]

  if (!length(nz)) {
    message("No non-zero values on the negative control. Background is at or ",
            "below the detection floor; a threshold of >= 1 is defensible.")
    return(list(threshold = 1, n_cells = length(x), n_nonzero = 0,
                mean_nonzero = 0, sd_nonzero = NA_real_,
                rule = "no non-zero background; threshold set to 1"))
  }

  m <- mean(nz)
  s <- stats::sd(nz)
  list(threshold    = m + n_sd * ifelse(is.na(s), 0, s),
       n_cells      = length(x),
       n_nonzero    = length(nz),
       pct_nonzero  = 100 * length(nz) / length(x),
       mean_nonzero = m,
       sd_nonzero   = s,
       rule = sprintf("mean + %g SD of non-zero puncta on negative control", n_sd))
}

#' Sweep the positivity threshold across animals.
#'
#' If group separation is flat across the sweep, the threshold is not driving
#' the result. If it collapses, the result depends on a number that was chosen
#' rather than measured. Either way the reader should see the curve.
qc_threshold_sweep <- function(cell_data, count_col, thresholds = 1:6,
                               by = c("group", "animal_id", "region")) {
  by <- intersect(by, names(cell_data))
  purrr::map_dfr(thresholds, function(k) {
    cell_data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
      dplyr::summarise(threshold = k,
                       n_cells = dplyr::n(),
                       pct_positive = 100 * mean(.data[[count_col]] >= k,
                                                 na.rm = TRUE),
                       .groups = "drop")
  })
}

#' Per-slide staining and intensity metrics.
#'
#' A slide that developed or stained differently moves every downstream number.
#' Plotting the primary readout against these separates biology from chemistry,
#' and the columns needed are usually already in the export.
qc_staining_drift <- function(cell_data, metric_cols = NULL,
                              by = c("group", "animal_id", "section_id")) {
  by <- intersect(by, names(cell_data))

  if (is.null(metric_cols)) {
    metric_cols <- grep("_od_mean$|_intensity$|hematoxylin|haematoxylin",
                        names(cell_data), value = TRUE)
  }
  metric_cols <- intersect(metric_cols, names(cell_data))

  if (!length(metric_cols)) {
    message("No optical-density or intensity columns found; staining drift QC ",
            "skipped. Exporting them costs nothing and is the fastest way to ",
            "separate biology from staining.")
    return(tibble::tibble())
  }

  cell_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(metric_cols),
                    list(median = ~stats::median(.x, na.rm = TRUE),
                         p90    = ~stats::quantile(.x, 0.90, na.rm = TRUE)),
                    .names = "{.col}__{.fn}"),
      .groups = "drop")
}

#' Per-slide cluster rate, flagged against the group median.
#'
#' The fraction of cells containing at least one cluster is a sensitive
#' indicator of signal density. Two animals in the same group differing by an
#' order of magnitude are not staining replicates of each other.
qc_cluster_rate <- function(cell_data, cluster_col = "n_clusters",
                            by = c("group", "animal_id", "region")) {
  by <- intersect(by, names(cell_data))

  out <- cell_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     pct_cluster_pos = 100 * mean(.data[[cluster_col]] >= 1,
                                                  na.rm = TRUE),
                     .groups = "drop")

  grp <- intersect(c("group", "region"), names(out))
  if (length(grp)) {
    out <- out |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
      dplyr::mutate(
        group_median = stats::median(pct_cluster_pos, na.rm = TRUE),
        fold_vs_group_median = pct_cluster_pos / pmax(group_median, 1e-6),
        flag = dplyr::if_else(
          fold_vs_group_median > 10 | fold_vs_group_median < 0.1,
          "CHECK STAINING", "")) |>
      dplyr::ungroup()
  }
  out
}

#' Cells and cells per mm2 per stratum.
#'
#' Anything per unit area is confounded by cellularity, so cells/mm2 has to
#' travel beside it. Anything per cell depends entirely on segmentation, so the
#' cell count has to be visible too.
qc_denominators <- function(cell_data, area_col = NULL,
                            by = c("group", "animal_id", "region")) {
  by <- intersect(by, names(cell_data))

  out <- cell_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  if (!is.null(area_col) && area_col %in% names(cell_data)) {
    areas <- cell_data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
      dplyr::summarise(area_mm2 = dplyr::first(.data[[area_col]]), .groups = "drop")
    out <- out |>
      dplyr::left_join(areas, by = by) |>
      dplyr::mutate(cells_per_mm2 = n_cells / area_mm2)
  } else {
    out$area_mm2      <- NA_real_
    out$cells_per_mm2 <- NA_real_
    message("No annotation area supplied. Without it, differences in how much ",
            "tissue was analysed are invisible.")
  }
  out
}

write_qc_tables <- function(qc_list, qc_dir) {
  for (nm in names(qc_list)) {
    x <- qc_list[[nm]]
    if (is.data.frame(x)) {
      if (nrow(x)) readr::write_csv(x, file.path(qc_dir, paste0(nm, ".csv")))
    } else if (is.list(x)) {
      yaml::write_yaml(x, file.path(qc_dir, paste0(nm, ".yml")))
    }
  }
  invisible(names(qc_list))
}

## -----------------------------------------------------------------------------
## Phenotype partition, marker failure, and effective coverage
## -----------------------------------------------------------------------------

#' Do the mutually exclusive phenotype counts add up to the total?
#'
#' A phenotype set that does not partition the detected cells means either a
#' phenotype is missing from the export, or the phenotypes overlap, or some
#' cells fall into none of them. All three change every percentage computed from
#' them.
qc_phenotype_partition <- function(tbl, phenotype_cols, total_col = "total_cells",
                                   id_cols = c("animal_id", "region")) {
  id_cols        <- intersect(id_cols, names(tbl))
  phenotype_cols <- intersect(phenotype_cols, names(tbl))

  tbl |>
    dplyr::mutate(
      phenotype_total = rowSums(dplyr::across(dplyr::all_of(phenotype_cols)),
                                na.rm = TRUE),
      partition_difference = .data[[total_col]] - phenotype_total,
      partition_ok = partition_difference == 0) |>
    dplyr::select(dplyr::all_of(c(id_cols, total_col)),
                  phenotype_total, partition_difference, partition_ok)
}

#' Distinguish a marker count of zero from a missing marker measurement.
#'
#' A stratum where an antibody yielded zero positive cells while every other
#' animal in that region yielded thousands is not a biological zero. It is a
#' staining or classification failure, and it must become NA rather than
#' entering a figure as a point at the bottom of the axis, where it reads as
#' "this animal had none of that cell type".
#'
#' @param min_fraction_of_median a count below this fraction of the stratum
#'   median is treated as a failure rather than a low value
qc_marker_failure <- function(tbl, marker_col, by = "region",
                              min_fraction_of_median = 0.05) {
  by <- intersect(by, names(tbl))

  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::mutate(
      stratum_median = stats::median(
        .data[[marker_col]][.data[[marker_col]] > 0], na.rm = TRUE),
      fraction_of_median = .data[[marker_col]] / stratum_median,
      marker_failed = .data[[marker_col]] == 0 |
        (is.finite(fraction_of_median) &
           fraction_of_median < min_fraction_of_median),
      failure_note = dplyr::case_when(
        .data[[marker_col]] == 0 ~ paste0(
          "ZERO marker-positive cells: staining or classification failure, ",
          "not a biological zero"),
        marker_failed ~ "Marker count far below the stratum median; check the slide",
        TRUE ~ "")) |>
    dplyr::ungroup()
}

#' Set every marker-derived quantity to NA where the marker failed.
#'
#' Percentages already come out NA through division by zero, which is why they
#' vanish from a figure with a "Removed n rows" warning. Densities and raw
#' fractions do not: they come out as a clean zero and are plotted. This makes
#' the treatment consistent, so a failed slide is absent everywhere rather than
#' absent from some panels and present as zero in others.
mask_failed_marker <- function(tbl, affected_cols, failure_col = "marker_failed") {
  affected_cols <- intersect(affected_cols, names(tbl))
  n <- sum(tbl[[failure_col]], na.rm = TRUE)
  if (n) {
    message("Masking ", length(affected_cols), " marker-derived column(s) in ",
            n, " stratum/strata where the marker failed.")
    for (cc in affected_cols) tbl[[cc]][tbl[[failure_col]]] <- NA
  }
  tbl
}

#' Report the design coverage that survives exclusions and marker failures.
#'
#' In a design with one animal per group per sex, losing one animal from a
#' region turns a level into a single observation of a single sex.
qc_effective_coverage <- function(tbl, by = "region", group_col = "group",
                                  sex_col = "sex", invalid_col = NULL) {
  by <- intersect(by, names(tbl))
  d  <- tbl
  if (!is.null(invalid_col) && invalid_col %in% names(d)) {
    d <- d[!d[[invalid_col]] %in% TRUE, , drop = FALSE]
  }

  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, group_col)))) |>
    dplyr::summarise(
      n_animals = dplyr::n_distinct(animal_id),
      sexes = if (sex_col %in% names(d))
        paste(sort(unique(as.character(.data[[sex_col]]))), collapse = "/")
      else NA_character_,
      .groups = "drop") |>
    dplyr::mutate(flag = dplyr::case_when(
      n_animals == 0 ~ "LEVEL LOST in this region",
      n_animals == 1 ~ "single animal; one sex only",
      TRUE ~ ""))
}

## =============================================================================
## SECTION 6 -- AGGREGATION TO THE ANALYSIS UNIT
##
## Cells are nested in sections, sections in animals. Cells from one animal are
## not independent observations of the treatment; the animal is. Everything
## below is explicit about what one row represents.
## =============================================================================

#' Assign puncta-count bins for an H-score.
#'
#' Two properties matter.
#'
#' Exhaustiveness. Branch conditions such as `x == 0` and `x >= 1 & x <= 3`
#' applied to a continuous, cluster-derived estimate leave values between 0 and
#' 1 matching no branch. Those cells become NA and are then dropped from the
#' denominator, and since a treated group holds more of them than a control, the
#' bias is directional. `cut()` with explicit breaks cannot do that.
#'
#' Provenance. Scoring scales for this assay are not standardised, and published
#' vendor documentation is internally inconsistent about them. Whatever scale
#' you use, name it here and define it in the methods every time.
#'
#' @param x per-cell values; floored first if not integer-valued
#' @param breaks lower bounds of each bin above zero
#' @param scale_name recorded in the run manifest
assign_bins <- function(x, breaks = c(1, 4, 10),
                        scale_name = "custom_0_3_declared_in_methods") {
  if (any(abs(x - round(x)) > 1e-8, na.rm = TRUE)) x <- floor(x)

  bin <- cut(x, breaks = c(-Inf, 0, breaks, Inf), right = TRUE,
             labels = FALSE) - 1L
  bin[is.na(x)] <- NA_integer_

  attr(bin, "scale_name") <- scale_name
  attr(bin, "breaks")     <- breaks
  attr(bin, "max_bin")    <- length(breaks)
  bin
}

#' H-score: sum over bins of bin multiplied by per cent of cells in that bin.
#'
#' Range is 0 to 100 times the top bin. Bin boundaries are roughly geometric, so
#' the score behaves like a log summary: differences are not proportional to
#' abundance, and bin 0 contributes nothing, so two tissues with very different
#' fractions of expressing cells can score identically. Where per-cell counts
#' exist, model the counts and compute an H-score afterwards only for
#' comparability with the pathology literature.
#'
#' @param na_policy "error" refuses to renormalise over a silently shrunken
#'   denominator; "drop" proceeds and reports how many cells were lost.
calc_hscore <- function(bin_vec, na_policy = c("error", "drop")) {
  na_policy <- match.arg(na_policy)
  n_na <- sum(is.na(bin_vec))

  if (n_na > 0) {
    if (na_policy == "error") {
      stop(n_na, " cells have an undefined bin; the denominator would shrink ",
           "silently. Fix the bin coverage, or pass na_policy = 'drop' and ",
           "report the dropped fraction.", call. = FALSE)
    }
    warning(n_na, " cells dropped from the H-score denominator.", call. = FALSE)
  }

  b <- bin_vec[!is.na(bin_vec)]
  if (!length(b)) return(NA_real_)
  n <- length(b)
  sum(vapply(sort(unique(b)), function(k) k * 100 * sum(b == k) / n, numeric(1)))
}

#' Summarise cells to one row per animal (and region), the analysis unit.
#'
#' The two hurdle components are reported separately because they answer
#' different biological questions:
#'   pct_positive     -- did the treatment change how many cells express?
#'   median_among_pos -- did it change output per expressing cell?
#' A test on the overall mean cannot tell them apart.
summarise_to_animal <- function(cell_data, count_col, threshold = 1,
                                by = c("group", "animal_id", "region"),
                                bins = c(1, 4, 10),
                                bin_scale_name = "custom_0_3_declared_in_methods") {
  by <- intersect(by, names(cell_data))
  if (!"animal_id" %in% by) {
    stop("`by` must include animal_id: the animal is the analysis unit.",
         call. = FALSE)
  }

  d <- dplyr::mutate(cell_data,
                     .count = .data[[count_col]],
                     .pos   = .data[[count_col]] >= threshold)
  if (!is.null(bins)) {
    d$.bin <- assign_bins(d$.count, breaks = bins, scale_name = bin_scale_name)
  }

  out <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_cells          = dplyr::n(),
      n_positive       = sum(.pos, na.rm = TRUE),
      pct_positive     = safe_pct(.pos),                    # hurdle part 1
      mean_among_pos   = safe_mean(.count[.pos]),           # hurdle part 2
      median_among_pos = safe_median(.count[.pos]),
      mean_all_cells   = safe_mean(.count),
      sd_all_cells     = safe_sd(.count),
      p95_all_cells    = stats::quantile(.count, 0.95, na.rm = TRUE),
      max_all_cells    = suppressWarnings(max(.count, na.rm = TRUE)),
      hscore = if (!is.null(bins)) calc_hscore(.bin, na_policy = "drop")
               else NA_real_,
      .groups = "drop")

  attr(out, "analysis_unit")        <- paste("one row =", paste(by, collapse = " x "))
  attr(out, "count_col")            <- count_col
  attr(out, "positivity_threshold") <- threshold
  out
}

#' Flag summary statistics that are degenerate across animals.
#'
#' An overall median is pinned to a small integer whenever more than half the
#' cells sit at zero. A ratio built on it can then take only a handful of
#' values, sometimes only two, which makes it a coin flip rather than a
#' measurement. This check makes that visible before it becomes a figure.
check_degenerate_summaries <- function(animal_tbl,
                                       cols = c("median_among_pos", "hscore",
                                                "mean_all_cells"),
                                       by = c("group", "region")) {
  by   <- intersect(by, names(animal_tbl))
  cols <- intersect(cols, names(animal_tbl))
  if (!length(cols)) return(tibble::tibble())

  animal_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(cols),
                    list(n_distinct = ~dplyr::n_distinct(.x, na.rm = TRUE),
                         n          = ~sum(!is.na(.x))),
                    .names = "{.col}__{.fn}"),
      .groups = "drop") |>
    tidyr::pivot_longer(-dplyr::all_of(by),
                        names_to = c("statistic", ".value"), names_sep = "__") |>
    dplyr::mutate(flag = dplyr::case_when(
      n >= 3 & n_distinct <= 2 ~ "DEGENERATE: two or fewer distinct values",
      n >= 4 & n_distinct <= 3 ~ "quantised: few distinct values",
      TRUE ~ ""))
}

#' Aggregate to one value per animal across regions.
#'
#' Used as a sanity check against a mixed model. Disagreement in DIRECTION means
#' the model is misspecified; the model having more power is expected.
pseudobulk <- function(cell_data, count_col, by = c("group", "animal_id")) {
  by <- intersect(by, names(cell_data))
  cell_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(n_cells      = dplyr::n(),
                     mean_count   = safe_mean(.data[[count_col]]),
                     pct_positive = 100 * mean(.data[[count_col]] >= 1,
                                               na.rm = TRUE),
                     .groups = "drop")
}

## =============================================================================
## SECTION 7 -- EXPRESSING RESULTS RELATIVE TO CONTROL
##
## Nothing here clips. A one-sided clip such as pmax(100 - remaining, 0)
## censors every animal above the control mean to exactly zero while animals
## below it contribute their full value. The mean of that is biased upward
## whenever there is between-animal spread, and it can never be negative, so it
## reports an effect even when the groups are identical. A negative value is
## information.
## =============================================================================

#' Express values relative to a control reference, without censoring.
#'
#' @param within columns defining the stratum each reference applies to
relative_to_control <- function(animal_tbl, value_col, cfg, within = c("region")) {

  ctrl   <- cfg$design$control_group
  within <- intersect(within, names(animal_tbl))
  if (!"group" %in% names(animal_tbl)) {
    stop("relative_to_control() needs a `group` column.", call. = FALSE)
  }

  ref <- animal_tbl |>
    dplyr::filter(as.character(group) == ctrl) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(within))) |>
    dplyr::summarise(control_mean = safe_mean(.data[[value_col]]),
                     control_sd   = safe_sd(.data[[value_col]]),
                     n_control    = sum(is.finite(.data[[value_col]])),
                     .groups = "drop") |>
    dplyr::mutate(control_cv_pct = 100 * control_sd / control_mean)

  out <- if (length(within)) {
    dplyr::left_join(animal_tbl, ref, by = within)
  } else {
    dplyr::mutate(animal_tbl,
                  control_mean   = ref$control_mean[1],
                  control_sd     = ref$control_sd[1],
                  n_control      = ref$n_control[1],
                  control_cv_pct = ref$control_cv_pct[1])
  }

  out <- dplyr::mutate(out,
    ratio_to_control = .data[[value_col]] / control_mean,
    pct_remaining    = 100 * ratio_to_control,
    pct_change       = 100 * (ratio_to_control - 1),
    pct_reduction    = 100 * (1 - ratio_to_control))   # may be negative

  ## A ratio against a near-zero denominator is not a usable quantity: its value
  ## is set by however many control cells happened to trip the threshold.
  tiny <- is.finite(out$control_mean) &
    (abs(out$control_mean) < .Machine$double.eps^0.5 |
       abs(out$control_mean) < 0.01 * safe_mean(abs(out[[value_col]])))

  if (any(tiny, na.rm = TRUE)) {
    for (cc in c("ratio_to_control", "pct_remaining", "pct_change", "pct_reduction")) {
      out[[cc]][tiny] <- NA_real_
    }
    warning("Control reference is at or near zero in ", sum(tiny, na.rm = TRUE),
            " stratum/strata; ratios set to NA. Report the DIFFERENCE on the ",
            "original scale with a confidence interval instead of a fold change.",
            call. = FALSE)
  }

  if (any(out$n_control < 2, na.rm = TRUE)) {
    warning("One or more strata have a single control animal, so the reference ",
            "carries no error estimate and any oddity in that animal ",
            "propagates to every treated value. Present raw values as primary.",
            call. = FALSE)
  }

  attr(out, "normalised_col") <- value_col
  attr(out, "control_group")  <- ctrl
  out
}

#' What does "no effect" look like in these data?
#'
#' Applies the same normalisation to the control animals with a leave-one-out
#' reference, so each control animal is compared against the others rather than
#' partly against itself. Plot this beside every "% remaining" figure: if the
#' treated points sit inside the control spread there is no effect to report,
#' whatever the point estimate says.
control_null_distribution <- function(animal_tbl, value_col, cfg,
                                      within = c("region"),
                                      leave_one_out = TRUE) {
  ctrl   <- cfg$design$control_group
  within <- intersect(within, names(animal_tbl))

  ctrl_tbl <- dplyr::filter(animal_tbl, as.character(group) == ctrl)
  if (!nrow(ctrl_tbl)) {
    warning("No control animals found for group '", ctrl, "'.", call. = FALSE)
    return(tibble::tibble())
  }

  grouped <- if (length(within)) {
    dplyr::group_by(ctrl_tbl, dplyr::across(dplyr::all_of(within)))
  } else dplyr::group_by(ctrl_tbl)

  grouped |>
    dplyr::mutate(
      .n_ref = sum(is.finite(.data[[value_col]])),
      .ref = if (leave_one_out && .n_ref[1] > 1) {
        (sum(.data[[value_col]], na.rm = TRUE) - .data[[value_col]]) / (.n_ref - 1)
      } else safe_mean(.data[[value_col]])) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      ratio_to_control = .data[[value_col]] / .ref,
      pct_remaining    = 100 * ratio_to_control,
      pct_reduction    = 100 * (1 - ratio_to_control),
      reference_type   = if (leave_one_out) "leave_one_out_control"
                         else "control_mean") |>
    dplyr::select(-.n_ref, -.ref)
}

#' Summarise the apparent effect the pipeline produces from noise alone.
null_effect_summary <- function(null_tbl, within = c("region")) {
  within <- intersect(within, names(null_tbl))
  null_tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(within))) |>
    dplyr::summarise(
      n_control_animals       = dplyr::n(),
      null_mean_pct_reduction = safe_mean(pct_reduction),
      null_sd_pct_reduction   = safe_sd(pct_reduction),
      null_min                = suppressWarnings(min(pct_reduction, na.rm = TRUE)),
      null_max                = suppressWarnings(max(pct_reduction, na.rm = TRUE)),
      .groups = "drop")
}

#' Housekeeping probes are an inclusion gate, not a divisor.
#'
#' Degradation hits high-expressing housekeepers hardest and low-to-moderate
#' expressors least, so the control's response to a pre-analytical insult is not
#' proportional to the target's, and the ratio adds noise instead of cancelling
#' artefact. If a sample fails, drop it; do not rescue it by dividing.
housekeeping_inclusion_gate <- function(animal_tbl, hk_col, min_value,
                                        id_cols = c("animal_id", "region")) {
  id_cols <- intersect(id_cols, names(animal_tbl))
  gate <- animal_tbl |>
    dplyr::select(dplyr::all_of(c(id_cols, hk_col))) |>
    dplyr::mutate(passes_hk_gate = .data[[hk_col]] >= min_value,
                  gate_rule = sprintf("%s >= %g", hk_col, min_value))

  n_fail <- sum(!gate$passes_hk_gate, na.rm = TRUE)
  if (n_fail) message(n_fail, " sample(s) fail the housekeeping gate and should ",
                      "be EXCLUDED, not normalised.")
  gate
}

## =============================================================================
## SECTION 8 -- MODELS
##
## Four properties of these data drive everything here.
##
## Nesting. Cells sit in sections, sections in animals. Treating cells as
##   independent observations of the treatment inflates type I error severely,
##   and using animal as a batch covariate does not fix it. Power comes from
##   animals, not cells.
##
## Distribution. Transcription is bursty, so steady-state counts are negative
##   binomial rather than Poisson. Cell heterogeneity, section-plane truncation
##   and detection-efficiency variation across a slide add further
##   overdispersion.
##
## Two kinds of zero. Structural (the cell does not express) and sampling (it
##   does, but this plane missed it). A single negative binomial cannot separate
##   them; a hurdle can, and its two parts are exactly the two quantities
##   usually reported by hand.
##
## Exposure. Unequal areas or cell counts go in as offset(log(exposure)), never
##   as a pre-computed ratio, which forces the exposure coefficient to exactly 1
##   and discards the precision information.
## =============================================================================

.need_glmmTMB <- function() {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB is required for this step. install.packages('glmmTMB')",
         call. = FALSE)
  }
}

#' Can this design support group inference at all?
#'
#' With n animals per group the smallest achievable two-sided Wilcoxon p is
#' 2 / choose(n1 + n2, n1). At 3 vs 3 that is 0.10, so no amount of separation
#' can reach 0.05. With one animal per group nothing is possible.
design_power_note <- function(cfg) {
  n <- dplyr::count(cfg$design$animals_tbl, group, name = "n")
  min_p <- if (nrow(n) == 2) 2 / choose(sum(n$n), n$n[1]) else NA_real_

  list(
    animals_per_group = as.list(stats::setNames(n$n, as.character(n$group))),
    min_two_sided_wilcoxon_p = min_p,
    inference_possible = all(n$n >= 2),
    note = if (all(n$n >= 2)) {
      if (is.na(min_p)) "More than two groups; use model contrasts."
      else sprintf(paste0("Smallest achievable two-sided Wilcoxon p with these ",
                          "group sizes is %.4f."), min_p)
    } else {
      paste0("At least one group has a single animal. No group-level inference ",
             "is possible. Report descriptively; do not draw boxplots over one ",
             "or two points and do not fit group contrasts.")
    })
}

#' Hurdle negative-binomial GLMM on per-cell counts.
#'
#' Fitting both parts says WHICH COMPONENT moved. A treatment that recruits more
#' cells into expression is different biology from one that raises output per
#' expressing cell, and a test on mean counts per cell cannot distinguish them.
#'
#' @param count_col integer counts. A cluster-derived estimate is an index, not
#'   a count; floor it first or use a genuine single-spot count.
fit_hurdle_nb <- function(cell_data, count_col, fixed = "group",
                          random = "(1 | animal_id/section_id)",
                          zi = "~ group", offset_col = NULL,
                          family = c("nbinom2", "nbinom1")) {
  .need_glmmTMB()
  family <- match.arg(family)

  y <- cell_data[[count_col]]
  if (any(abs(y - round(y)) > 1e-8, na.rm = TRUE)) {
    stop(count_col, " is not integer-valued. Use a genuine count column, or ",
         "floor() the estimate and record that in the methods.", call. = FALSE)
  }

  rhs <- paste(fixed, random, sep = " + ")
  if (!is.null(offset_col)) rhs <- paste0(rhs, " + offset(log(", offset_col, "))")

  glmmTMB::glmmTMB(
    stats::as.formula(paste(count_col, "~", rhs)),
    family    = if (family == "nbinom2") glmmTMB::nbinom2() else glmmTMB::nbinom1(),
    ziformula = stats::as.formula(zi),
    data      = cell_data)
}

#' Compare Poisson, NB1 and NB2 by AIC. Assume NB, but test it.
compare_count_families <- function(cell_data, count_col, fixed = "group",
                                   random = "(1 | animal_id/section_id)",
                                   zi = "~ 1") {
  .need_glmmTMB()
  form <- stats::as.formula(paste(count_col, "~", paste(fixed, random, sep = " + ")))
  zif  <- stats::as.formula(zi)

  fams <- list(poisson = stats::poisson(),
               nbinom1 = glmmTMB::nbinom1(),
               nbinom2 = glmmTMB::nbinom2())

  purrr::imap_dfr(fams, function(f, nm) {
    fit <- try(glmmTMB::glmmTMB(form, family = f, ziformula = zif,
                                data = cell_data), silent = TRUE)
    ok <- !inherits(fit, "try-error")
    tibble::tibble(family = nm, converged = ok,
                   AIC    = if (ok) stats::AIC(fit) else NA_real_,
                   logLik = if (ok) as.numeric(stats::logLik(fit)) else NA_real_)
  }) |>
    dplyr::arrange(AIC)
}

#' Binomial GLMM on positivity, with cells as trials.
#'
#' Watch for complete separation. When a control group sits essentially at zero
#' the model will not converge sensibly, and it should not: report the
#' descriptive difference with an exact interval and reserve the model for
#' contrasts within the treated animals.
fit_positivity_glmm <- function(animal_tbl, n_pos_col = "n_positive",
                                n_tot_col = "n_cells", fixed = "group",
                                random = "(1 | animal_id)") {
  .need_glmmTMB()
  d <- dplyr::mutate(animal_tbl, .n_neg = .data[[n_tot_col]] - .data[[n_pos_col]])

  fit <- glmmTMB::glmmTMB(
    stats::as.formula(paste0("cbind(", n_pos_col, ", .n_neg) ~ ",
                             fixed, " + ", random)),
    family = stats::binomial(), data = d)

  rng <- d |>
    dplyr::group_by(group) |>
    dplyr::summarise(p = mean(.data[[n_pos_col]] / .data[[n_tot_col]]),
                     .groups = "drop")
  if (any(rng$p < 1e-4, na.rm = TRUE) || any(rng$p > 1 - 1e-4, na.rm = TRUE)) {
    warning("At least one group is at or near 0% or 100% positive. The model is ",
            "close to complete separation; treat its intervals with suspicion ",
            "and report descriptives.", call. = FALSE)
  }
  fit
}

#' Intraclass correlation.
#'
#' The number that makes the design point land: at an ICC of 0.3, 500 cells per
#' animal carry roughly the information of three independent observations.
report_icc <- function(fit) {
  if (!requireNamespace("performance", quietly = TRUE)) {
    message("performance not installed; ICC skipped."); return(NULL)
  }
  out <- try(performance::icc(fit), silent = TRUE)
  if (inherits(out, "try-error")) {
    message("ICC could not be computed for this model."); return(NULL)
  }
  out
}

#' Dose-response trend with dose as a continuous predictor.
#'
#' Where a design has one animal per dose there is no within-group replication
#' and no group contrast is possible. A TREND is still estimable and uses every
#' animal, because the slope is informed by all of them. This is the analysis
#' such a design supports.
fit_dose_trend <- function(animal_tbl, value_col, dose_col = "dose",
                           covariates = character(0),
                           random = "(1 | animal_id)", log_dose = TRUE) {
  .need_glmmTMB()
  d <- animal_tbl

  dose_term <- if (log_dose) { d$.dose_term <- log1p(d[[dose_col]]); ".dose_term" }
               else dose_col

  glmmTMB::glmmTMB(
    stats::as.formula(paste(value_col, "~",
                            paste(c(dose_term, covariates, random), collapse = " + "))),
    family = stats::gaussian(), data = d)
}

#' Tidy a fitted model, exponentiating log-link coefficients into ratios.
tidy_model <- function(fit, exponentiate = TRUE, conf_level = 0.95) {
  ct <- summary(fit)$coefficients$cond
  if (is.null(ct)) return(tibble::tibble())
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  out <- tibble::tibble(term = rownames(ct), estimate = ct[, 1],
                        std_error = ct[, 2], statistic = ct[, 3],
                        p_value = ct[, 4]) |>
    dplyr::mutate(conf_low  = estimate - z * std_error,
                  conf_high = estimate + z * std_error)

  if (exponentiate) {
    dplyr::mutate(out,
                  estimate  = exp(estimate),
                  conf_low  = exp(conf_low),
                  conf_high = exp(conf_high),
                  scale = "ratio (exponentiated log-link coefficient)")
  } else {
    dplyr::mutate(out, scale = "link scale")
  }
}

#' Benjamini-Hochberg across the declared grid.
#'
#' Declare the grid before analysis. Regions within an animal are correlated, so
#' where region structure is the question, model it (region fixed effect plus
#' animal random effect, then contrasts) rather than fitting separate models per
#' region. Correcting within a figure panel but not across panels is not a
#' correction.
adjust_multiplicity <- function(results_tbl, p_col = "p_value",
                                grid_cols = c("region"), method = "BH") {
  grid_cols <- intersect(grid_cols, names(results_tbl))
  dplyr::mutate(results_tbl,
    p_adjusted = stats::p.adjust(.data[[p_col]], method = method),
    adjustment = sprintf("%s across %d tests over: %s", method, dplyr::n(),
                         paste(grid_cols, collapse = " x ")))
}

#' Compare a mixed model against an animal-level test on aggregated values.
#'
#' Disagreement in direction means the model is misspecified. The model having
#' more power is expected and is not itself a problem.
pseudobulk_check <- function(pseudobulk_tbl, value_col = "mean_count",
                             group_col = "group") {
  g  <- as.character(pseudobulk_tbl[[group_col]])
  lv <- unique(g)
  if (length(lv) != 2) {
    message("Pseudobulk check needs exactly two groups; skipped.")
    return(NULL)
  }
  a <- pseudobulk_tbl[[value_col]][g == lv[1]]
  b <- pseudobulk_tbl[[value_col]][g == lv[2]]
  wt <- suppressWarnings(stats::wilcox.test(b, a))

  list(groups = lv, n = c(length(a), length(b)),
       median = c(stats::median(a, na.rm = TRUE), stats::median(b, na.rm = TRUE)),
       ratio  = mean(b, na.rm = TRUE) / mean(a, na.rm = TRUE),
       W = unname(wt$statistic), p_value = wt$p.value,
       min_possible_p = 2 / choose(length(a) + length(b), length(a)))
}

## =============================================================================
## SECTION 9 -- FIGURES
##
## Two rules.
##   One point is one animal. Cell-level distributions are QC, never results.
##   No box over fewer than four points, because a box across two animals
##   invents quartiles. The plotting function falls back to points and a mean
##   bar and says so in the caption.
## =============================================================================

theme_rnascope <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key        = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position   = "bottom",
      strip.background  = ggplot2::element_rect(fill = "white", colour = "black",
                                                linewidth = 0.6),
      strip.text        = ggplot2::element_text(face = "bold", colour = "black"),
      axis.text         = ggplot2::element_text(colour = "black"),
      axis.title        = ggplot2::element_text(colour = "black"),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = 0.5),
      plot.title        = ggplot2::element_text(face = "bold", colour = "black"),
      plot.subtitle     = ggplot2::element_text(colour = "grey30"),
      plot.caption      = ggplot2::element_text(colour = "grey40", hjust = 0))
}

theme_rnascope_slide <- function(base_size = 14) {
  theme_rnascope(base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title    = ggplot2::element_text(size = 15),
      axis.text     = ggplot2::element_text(size = 13),
      strip.text    = ggplot2::element_text(size = 14, face = "bold"))
}

#' Save a figure, using ragg where available so UTF-8 labels survive.
save_plot <- function(plot, filename, plot_dir, width = 9, height = 6, dpi = 300) {
  args <- list(filename = file.path(plot_dir, filename), plot = plot,
               width = width, height = height, dpi = dpi, bg = "white")
  if (requireNamespace("ragg", quietly = TRUE)) args$device <- ragg::agg_png
  do.call(ggplot2::ggsave, args)
  message("Saved: ", filename)
  invisible(file.path(plot_dir, filename))
}

#' Animal-level comparison that adapts to how many animals there are.
plot_animal_points <- function(animal_tbl, y_col, x_col = "group", facet = NULL,
                               shape_by = NULL, y_label = y_col, title = NULL,
                               subtitle = NULL, caption = NULL,
                               min_n_for_box = 4, colours = NULL) {

  per_cell <- animal_tbl |>
    dplyr::filter(is.finite(.data[[y_col]])) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(x_col, facet)))) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  use_box <- nrow(per_cell) > 0 && all(per_cell$n >= min_n_for_box)

  p <- ggplot2::ggplot(animal_tbl,
                       ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))

  if (use_box) {
    p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[x_col]]),
                                   outlier.shape = NA, alpha = 0.30, width = 0.6,
                                   colour = "grey30", linewidth = 0.4, fatten = 0)
  }

  ## The mean bar is always drawn: it is the summary a reader can check against
  ## the visible points.
  p <- p + ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                                 linewidth = 0.4, colour = "black")

  point_aes <- if (!is.null(shape_by)) ggplot2::aes(shape = .data[[shape_by]])
               else ggplot2::aes()
  p <- p + ggplot2::geom_point(point_aes,
                               position = ggplot2::position_jitter(width = 0.08,
                                                                   height = 0),
                               size = 2.6, alpha = 0.95, colour = "black")

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  }
  if (!is.null(colours)) {
    p <- p + ggplot2::scale_fill_manual(values = colours, guide = "none")
  }

  n_note <- paste0("One point = one animal (n = ",
                   paste(per_cell$n, collapse = ", "), " per group).")
  if (!use_box) {
    n_note <- paste(n_note, "Boxes omitted: fewer than", min_n_for_box,
                    "animals per group, so quartiles are not estimable.")
  }

  p + ggplot2::labs(x = NULL, y = y_label, title = title, subtitle = subtitle,
                    caption = paste(c(n_note, caption), collapse = "\n")) +
    theme_rnascope_slide()
}

#' Percent-remaining plot with the control null distribution behind it.
#'
#' The grey band is what no effect looks like given the noise in these data, so
#' the reader can judge the treated points against it directly.
plot_with_control_null <- function(relative_tbl, null_tbl, y_col = "pct_remaining",
                                   x_col = "group", facet = NULL,
                                   y_label = "% remaining vs control",
                                   title = NULL, subtitle = NULL) {

  band <- if (nrow(null_tbl)) {
    null_tbl |>
      dplyr::group_by(dplyr::across(dplyr::all_of(intersect(facet, names(null_tbl))))) |>
      dplyr::summarise(lo = suppressWarnings(min(.data[[y_col]], na.rm = TRUE)),
                       hi = suppressWarnings(max(.data[[y_col]], na.rm = TRUE)),
                       .groups = "drop")
  } else tibble::tibble()

  p <- ggplot2::ggplot(relative_tbl,
                       ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))

  if (nrow(band)) {
    p <- p + ggplot2::geom_rect(
      data = band,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
      inherit.aes = FALSE, fill = "grey70", alpha = 0.30)
  }

  p +
    ggplot2::geom_hline(yintercept = 100, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                          linewidth = 0.4, colour = "black") +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0),
                        size = 2.8, colour = "black") +
    (if (!is.null(facet)) ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))) +
    ggplot2::labs(x = NULL, y = y_label, title = title, subtitle = subtitle,
      caption = paste0(
        "Grey band: the same normalisation applied to the control animals ",
        "(leave-one-out reference).\nIt is what no effect looks like in these ",
        "data. Values are not clipped; above 100% and\nbelow 0% are both ",
        "possible and meaningful.")) +
    theme_rnascope_slide()
}

plot_threshold_sweep <- function(sweep_tbl, facet = "region", colour_by = "group") {
  ggplot2::ggplot(sweep_tbl,
                  ggplot2::aes(x = threshold, y = pct_positive,
                               group = animal_id, colour = .data[[colour_by]])) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    (if (!is.null(facet) && facet %in% names(sweep_tbl))
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_y")) +
    ggplot2::labs(x = "Positivity threshold (puncta per cell)",
                  y = "% positive cells",
                  title = "Sensitivity of positivity to the threshold",
      caption = paste0(
        "One line per animal. Flat group separation across the sweep means the ",
        "threshold is not driving\nthe result. Collapsing separation means it ",
        "depends on a number that was chosen, not measured.")) +
    theme_rnascope()
}

plot_staining_vs_readout <- function(animal_tbl, stain_col, readout_col,
                                     colour_by = "group", label_col = "animal_id") {
  ggplot2::ggplot(animal_tbl,
                  ggplot2::aes(x = .data[[stain_col]], y = .data[[readout_col]],
                               colour = .data[[colour_by]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = .data[[label_col]]),
                       vjust = -0.9, size = 3, show.legend = FALSE) +
    ggplot2::labs(x = stain_col, y = readout_col,
                  title = "Readout against per-slide staining metric",
      caption = paste0("A trend here means slide-to-slide staining contributes ",
                       "to the between-animal spread in the readout.")) +
    theme_rnascope()
}

## =============================================================================
## SECTION 10 -- PROVENANCE
##
## A pipeline run should describe itself: which config, which input files, which
## acquisition and detection settings, which package versions, which commit.
## Parameters the config leaves blank are listed in the manifest as unrecorded,
## which is itself a finding.
## =============================================================================

hash_inputs <- function(manifest) {
  manifest |>
    dplyr::mutate(sha256 = vapply(path, function(p)
      digest::digest(file = p, algo = "sha256"), character(1))) |>
    dplyr::select(source_file, animal_id, size_mb, sha256)
}

git_commit <- function() {
  out <- suppressWarnings(try(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    silent = TRUE))
  if (inherits(out, "try-error") || !length(out)) return(NA_character_)
  status <- suppressWarnings(try(
    system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
    silent = TRUE))
  dirty <- !inherits(status, "try-error") && length(status) > 0
  paste0(out[1], if (dirty) " (working tree dirty)" else "")
}

EXPECTED_PARAMS <- list(
  acquisition = c("objective_magnification", "z_handling", "section_thickness_um",
                  "exposure_fixed_across_slides", "batch_randomised"),
  detection   = c("software", "software_version", "cell_expansion_radius_um",
                  "min_spot_size", "spot_intensity_threshold",
                  "stain_vectors_fixed", "roi_selection", "operator_blinded"))

write_run_manifest <- function(cfg, dirs, manifest = NULL, extra = list()) {

  unknown <- character(0)
  for (sec in names(EXPECTED_PARAMS)) {
    for (k in EXPECTED_PARAMS[[sec]]) {
      if (is.null(cfg[[sec]][[k]])) unknown <- c(unknown, paste0(sec, "$", k))
    }
  }

  man <- list(
    run = list(timestamp  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
               user       = Sys.info()[["user"]],
               hostname   = Sys.info()[["nodename"]],
               git_commit = git_commit(),
               r_version  = paste(R.version$major, R.version$minor, sep = ".")),
    config = list(path = cfg$.config_path, sha256 = cfg$.config_hash,
                  study = cfg$study),
    parameters = list(acquisition = cfg$acquisition,
                      detection   = cfg$detection,
                      readout     = cfg$readout,
                      design_type = cfg$design$type),
    unrecorded_parameters = if (length(unknown)) unknown else "none",
    output_dir = dirs$root,
    extra = extra)

  yaml::write_yaml(man, file.path(dirs$run, "run_manifest.yml"))

  if (!is.null(manifest)) {
    readr::write_csv(hash_inputs(manifest),
                     file.path(dirs$run, "input_file_hashes.csv"))
  }
  writeLines(capture.output(utils::sessionInfo()),
             file.path(dirs$run, "sessionInfo.txt"))
  file.copy(cfg$.config_path,
            file.path(dirs$run, basename(cfg$.config_path)), overwrite = TRUE)

  if (length(unknown)) {
    warning("Parameters not recorded in the config, listed as unknown in the ",
            "manifest:\n  ", paste(unknown, collapse = "\n  "),
            "\nThese are the settings that silently determine results.",
            call. = FALSE)
  }

  message("Run manifest written to ", dirs$run)
  invisible(man)
}

#' Methods stub populated from the config and the run.
#'
#' Anything it cannot fill is left as an explicit TODO rather than omitted.
write_methods_stub <- function(cfg, dirs, animal_tbl = NULL, regime_tbl = NULL,
                               threshold_info = NULL) {
  g <- function(x, default = "TODO") if (is.null(x)) default else as.character(x)

  probes <- if (!is.null(cfg$assay$probes))
    paste(vapply(cfg$assay$probes, function(p) g(p$name), character(1)),
          collapse = ", ") else "TODO"

  lines <- c(
    "# Methods (auto-generated stub -- edit before use)", "",
    sprintf("Study: %s", g(cfg$study$name)),
    sprintf("Assay: %s", g(cfg$assay$type)),
    sprintf("Probe(s): %s", probes),
    sprintf("Detection chemistry: %s", g(cfg$assay$chemistry)), "",
    "## Controls",
    sprintf("- Positive control probe: %s", g(cfg$assay$positive_control_probe)),
    sprintf("- Negative control probe: %s", g(cfg$assay$negative_control_probe)),
    sprintf("- Control group: %s", g(cfg$design$control_group)), "",
    "## Tissue and acquisition",
    sprintf("- Section thickness: %s um", g(cfg$acquisition$section_thickness_um)),
    sprintf("- Objective: %s", g(cfg$acquisition$objective_magnification)),
    sprintf("- Instrument: %s", g(cfg$acquisition$instrument)),
    sprintf("- z handling: %s", g(cfg$acquisition$z_handling)),
    sprintf("- Exposure and gain fixed across slides: %s",
            g(cfg$acquisition$exposure_fixed_across_slides)),
    sprintf("- Slide batches randomised: %s", g(cfg$acquisition$batch_randomised)),
    "",
    paste0("Section thickness is a first-order multiplier on puncta per cell and ",
           "is not corrected for here. Counts are not comparable across studies ",
           "at different thickness."), "",
    "## Image analysis",
    sprintf("- Software: %s %s", g(cfg$detection$software),
            g(cfg$detection$software_version, "")),
    sprintf("- Segmentation: %s", g(cfg$detection$segmentation)),
    sprintf("- Cell expansion radius: %s um",
            g(cfg$detection$cell_expansion_radius_um)),
    sprintf("- Minimum spot size: %s", g(cfg$detection$min_spot_size)),
    sprintf("- Spot intensity threshold: %s",
            g(cfg$detection$spot_intensity_threshold)),
    sprintf("- Single-dot calibration: %s", g(cfg$readout$single_dot_calibration)),
    sprintf("- Deconvolution / unmixing vectors fixed across batch: %s",
            g(cfg$detection$stain_vectors_fixed)),
    sprintf("- ROI selection: %s", g(cfg$detection$roi_selection)),
    sprintf("- Operator blinded to group: %s", g(cfg$detection$operator_blinded)),
    "", "## Readout",
    sprintf("- Primary readout: %s", g(cfg$readout$primary)),
    sprintf("- Positivity threshold: %s (%s)",
            g(cfg$readout$positivity_threshold), g(cfg$readout$threshold_source)),
    sprintf("- Scoring scale: %s", g(cfg$readout$bin_scale_name)), "")

  if (!is.null(threshold_info)) {
    lines <- c(lines,
      sprintf("The positivity threshold was derived as %s, giving %.3f.",
              g(threshold_info$rule), threshold_info$threshold), "")
  }
  if (!is.null(regime_tbl) && nrow(regime_tbl)) {
    lines <- c(lines, sprintf(paste0(
      "Expression regime across strata: %s. Where the regime is not ",
      "'countable', per-cell values are a cluster-derived index and are ",
      "reported as such, never as copy number."),
      paste(unique(regime_tbl$regime), collapse = ", ")), "")
  }

  lines <- c(lines,
    "## Statistics",
    sprintf("- Analysis unit: one animal (n = %d).", nrow(cfg$design$animals_tbl)),
    "- Cells are nested in sections within animals; animal is a random effect.",
    "- Counts modelled as negative binomial with a log link; Poisson, NB1 and",
    "  NB2 compared by AIC. Counts were not log(x+1) transformed.",
    "- The hurdle components (fraction of positive cells, and counts among",
    "  positive cells) are reported separately.",
    "- Unequal exposure handled with offset(log(exposure)), not a ratio.",
    "- Intraclass correlation reported.",
    "- Multiplicity: BH-FDR across the declared grid.", "",
    "## Reporting",
    paste0("No reporting standard exists for quantitative ISH. This report ",
           "follows the Schmied 2024 image-analysis checklist, QUAREP-LiMi and ",
           "ARRIVE 2.0. See docs/METHODS_CHECKLIST.md."), "")

  path <- file.path(dirs$root, "METHODS_STUB.md")
  writeLines(lines, path)
  message("Methods stub written to ", path)
  invisible(path)
}

#' Interpretive limits for oligonucleotide localisation by ISH.
#'
#' Written next to the numbers so they end up in the same folder.
evidence_boundaries_ont <- function() {
  c("INTERPRETIVE LIMITS -- oligonucleotide detection by in situ hybridisation",
    "",
    "Detected compound is not active compound. Two uptake pathways exist, one",
    "  ending in the endolysosome where the oligonucleotide is inert. Knockdown",
    "  does not correlate with bulk intracellular accumulation; measured",
    "  endosomal escape is 1-2%, while fewer than 2,000 cytosolic copies",
    "  suffice to silence. Report PRESENCE and PATTERN. Perinuclear versus",
    "  peripheral distribution has distinguished productive from non-productive",
    "  uptake where total amount did not.",
    "",
    "No limit of detection has been published for oligonucleotide detection by",
    "  ISH, in any units. A negative region means 'below an unknown threshold',",
    "  not 'no compound'.",
    "",
    "No correlation of ISH signal against LC-MS or hybridisation-ELISA tissue",
    "  concentration has been published. Do not phrase the readout as a",
    "  concentration.",
    "",
    "No test exists of whether these assays distinguish full-length compound",
    "  from n-1, n-2 or chain-shortened metabolites. Only LC-MS/MS answers that.",
    "",
    "For siRNA, state which strand the probe targets. Passenger-strand",
    "  persistence is a weaker proxy for activity than guide strand, since the",
    "  passenger is discarded on RISC loading.",
    "",
    "No validation exists for LNA gapmers, cEt, GalNAc-siRNA or CNS-conjugated",
    "  siRNA, and no assessment of whether a terminal conjugate sterically",
    "  blocks probe binding.")
}

## =============================================================================
## SECTION 11 -- MAIN
## =============================================================================

main <- function(config_path) {

  setup_session()
  cfg  <- load_study_config(config_path)
  dirs <- setup_output_dirs(cfg)
  qc   <- list()

  LEVEL       <- cfg$detection$export_level %||% "summary"
  marker_name <- cfg$assay$antibodies[[1]]$name %||% "Marker"
  cols        <- cfg$export_columns
  if (is.null(cols)) {
    stop("This pipeline needs an export_columns section in the config.",
         call. = FALSE)
  }

  message("\n== Pipeline 04: compound distribution with ", marker_name, " ==")
  message(cfg$study$name, "  |  export level: ", LEVEL, "\n")

  if (identical(LEVEL, "summary")) {
    message("Summary-level input. The following are NOT available and are ",
            "skipped rather than approximated:\n",
            "  - per-cell count distributions and the hurdle model\n",
            "  - counts among positive cells\n",
            "  - subcellular (perinuclear vs peripheral) pattern\n",
            "  - local-density check on the ", marker_name, " classification\n",
            "  - signal not associated with any segmented cell\n",
            "Read the matching object-level export to recover these.\n")
  }

  ## --- 1. Inputs -------------------------------------------------------------
  default_pattern <- if (identical(LEVEL, "summary")) "_summary\\.csv$"
                     else "_object_detection\\.csv$"
  manifest <- build_file_manifest(
    input_dir = cfg$paths$input_dir,
    pattern   = cfg$paths$file_pattern %||% default_pattern,
    id_regex  = cfg$paths$id_regex)
  assert_roster(manifest, cfg)
  qc$file_manifest <- manifest

  if (identical(LEVEL, "summary")) {
    res <- read_all(manifest, function(p)
      read_summary_export(p, columns = cols$summary_columns,
                          phenotypes = cols$phenotypes,
                          delim = cfg$paths$delimiter %||% ","))
  } else {
    res <- read_all(manifest, function(p)
      read_object_export(p, marker_cols = cols$markers, probe_cols = cols$probes,
                         extra_cols = unname(unlist(cols$object_extras)),
                         delim = cfg$paths$delimiter %||% ","))
  }
  qc$read_qc <- res$qc
  tbl <- attach_design(res$data, cfg)

  ## --- Annotation recodes ----------------------------------------------------
  raw_region <- if (identical(LEVEL, "summary")) tbl$region
                else tbl[[make_safe_name(cols$region)]]

  rec <- apply_annotation_recodes(raw_region, tbl$animal_id,
                                  cfg$design$annotation_recodes)
  qc$annotation_recodes <- rec$log
  if (nrow(rec$log)) {
    message("\nAnnotation recodes applied:")
    print(rec$log)
    message("Each asserts the tissue is not what the annotation said. Confirm ",
            "against the image before publishing.\n")
  }

  tbl$region <- resolve_region(rec$values, cfg$design$regions,
                               cfg$design$region_aliases)
  unresolved <- is.na(tbl$region)
  qc$unresolved_regions <- tibble::tibble(
    animal_id = tbl$animal_id[unresolved], raw = rec$values[unresolved]) |>
    dplyr::count(animal_id, raw, sort = TRUE)
  if (nrow(qc$unresolved_regions)) {
    warning(nrow(qc$unresolved_regions), " annotation(s) could not be resolved ",
            "to a declared region. See qc/unresolved_regions.csv.", call. = FALSE)
  }
  tbl <- dplyr::filter(tbl, !is.na(region))

  meta_cols <- intersect(c("animal_id", "sex", "dose"),
                         names(cfg$design$animals_tbl))
  if (length(meta_cols) > 1) {
    tbl <- dplyr::left_join(tbl, dplyr::select(cfg$design$animals_tbl,
                                               dplyr::all_of(meta_cols)),
                            by = "animal_id")
  }

  ## --- 2. Build the animal x region table ------------------------------------
  ph <- c("marker_with_compound", "marker_without_compound",
          "nonmarker_with_compound", "double_negative")

  if (identical(LEVEL, "summary")) {
    pn <- names(cols$phenotypes)
    region_tbl <- tbl |>
      dplyr::rename(marker_with_compound    = !!pn[1],
                    marker_without_compound = !!pn[2],
                    nonmarker_with_compound = !!pn[3],
                    double_negative         = !!pn[4]) |>
      dplyr::mutate(
        marker_total          = marker_with_compound + marker_without_compound,
        compound_positive     = marker_with_compound + nonmarker_with_compound,
        marker_negative_total = nonmarker_with_compound + double_negative,
        area_mm2              = as.numeric(area_analyzed_um2) / 1e6)
  } else {
    marker_col   <- paste0(names(cols$markers)[1], "_pos")
    compound_col <- names(cols$probes)[1]
    thr          <- cfg$readout$positivity_threshold

    region_tbl <- tbl |>
      dplyr::mutate(compound_pos = .data[[compound_col]] >= thr) |>
      dplyr::group_by(group, animal_id, region,
                      dplyr::across(dplyr::any_of(c("sex", "dose")))) |>
      dplyr::summarise(
        total_cells             = dplyr::n(),
        marker_total            = sum(.data[[marker_col]], na.rm = TRUE),
        marker_with_compound    = sum(.data[[marker_col]] & compound_pos, na.rm = TRUE),
        marker_without_compound = sum(.data[[marker_col]] & !compound_pos, na.rm = TRUE),
        nonmarker_with_compound = sum(!.data[[marker_col]] & compound_pos, na.rm = TRUE),
        double_negative         = sum(!.data[[marker_col]] & !compound_pos, na.rm = TRUE),
        compound_positive       = sum(compound_pos, na.rm = TRUE),
        marker_negative_total   = sum(!.data[[marker_col]], na.rm = TRUE),
        area_mm2                = NA_real_,
        .groups = "drop")
  }

  ## Do the phenotypes partition the detected cells?
  qc$phenotype_partition <- qc_phenotype_partition(
    region_tbl, phenotype_cols = ph, total_col = "total_cells")

  bad <- dplyr::filter(qc$phenotype_partition, !partition_ok)
  if (nrow(bad)) {
    print(bad)
    warning(nrow(bad), " stratum/strata where the phenotype counts do not sum ",
            "to the total cell count. Either a phenotype is missing from the ",
            "export or the phenotypes overlap. Every percentage below depends ",
            "on this.", call. = FALSE)
  }

  ## ===========================================================================
  ## 3. Marker failure
  ##
  ## A stratum where the antibody yielded zero positive cells while every other
  ## animal in that region yielded thousands is a staining or classification
  ## failure. Percentages already come out NA through division by zero and
  ## vanish from a figure with a "Removed n rows" warning; densities and raw
  ## fractions come out as a clean zero and are plotted, where they read as
  ## "this animal had none of that cell type". This makes the treatment
  ## consistent.
  ## ===========================================================================

  region_tbl <- qc_marker_failure(region_tbl, marker_col = "marker_total",
                                  by = "region")

  qc$marker_failure <- dplyr::select(region_tbl, dplyr::any_of(c(
    "animal_id", "group", "region", "sex", "dose", "total_cells",
    "marker_total", "stratum_median", "fraction_of_median", "marker_failed",
    "failure_note")))

  failed <- dplyr::filter(qc$marker_failure, marker_failed)
  if (nrow(failed)) {
    message("\n", marker_name, " FAILED in ", nrow(failed), " stratum/strata:")
    print(dplyr::select(failed, dplyr::any_of(
      c("animal_id", "group", "region", "marker_total", "stratum_median",
        "failure_note"))))
    message("")
  }

  ## --- 4. Derived quantities -------------------------------------------------
  region_tbl <- region_tbl |>
    dplyr::mutate(
      ## Marker-independent. Unaffected by an antibody failure, because the
      ## phenotype partition is complete: total compound-positive is the sum of
      ## the two compound-positive phenotypes however the marker split them.
      pct_compound_all_cells       = 100 * compound_positive / total_cells,
      compound_density_per_mm2     = compound_positive / area_mm2,
      cells_per_mm2                = total_cells / area_mm2,
      ## Marker-dependent.
      pct_marker_of_all_cells      = 100 * marker_total / total_cells,
      marker_density_per_mm2       = marker_total / area_mm2,
      pct_marker_with_compound     = 100 * marker_with_compound / marker_total,
      pct_nonmarker_with_compound  = 100 * nonmarker_with_compound /
                                     marker_negative_total,
      marker_compound_density_per_mm2    = marker_with_compound / area_mm2,
      nonmarker_compound_density_per_mm2 = nonmarker_with_compound / area_mm2,
      ## Preferential occupancy. A ratio of 1 means the compound is distributed
      ## in proportion to cell numbers. It does not mean more compound per cell.
      compound_enrichment_in_marker = pct_marker_with_compound /
                                      pct_nonmarker_with_compound)

  region_tbl <- mask_failed_marker(region_tbl, affected_cols = c(
    "pct_marker_of_all_cells", "marker_density_per_mm2",
    "pct_marker_with_compound", "pct_nonmarker_with_compound",
    "marker_compound_density_per_mm2", "nonmarker_compound_density_per_mm2",
    "compound_enrichment_in_marker"))

  readr::write_csv(region_tbl,
                   file.path(dirs$tables, "animal_region_summary.csv"))

  ## --- 5. Effective coverage -------------------------------------------------
  qc$effective_coverage <- qc_effective_coverage(
    region_tbl, by = "region", group_col = "group", sex_col = "sex",
    invalid_col = "marker_failed")

  message("Effective coverage for ", marker_name, "-dependent endpoints:")
  print(qc$effective_coverage, n = Inf)

  lost <- dplyr::filter(qc$effective_coverage, nzchar(flag))
  if (nrow(lost)) {
    warning("After exclusions and marker failures, ", nrow(lost),
            " group x region cell(s) rest on a single animal of a single sex. ",
            "See qc/effective_coverage.csv.", call. = FALSE)
  }

  ## --- 6. Empirical background from the control animals ----------------------
  ## Control animals carry no compound, so whatever the classifier calls
  ## positive in them is the false-positive rate. That number, not zero, is the
  ## floor against which every treated value is read.
  ctrl <- dplyr::filter(region_tbl,
                        as.character(group) == cfg$design$control_group)

  qc$background <- ctrl |>
    dplyr::group_by(region) |>
    dplyr::summarise(
      n_control_animals            = dplyr::n(),
      bg_pct_compound_all_cells    = mean(pct_compound_all_cells, na.rm = TRUE),
      bg_pct_marker_with_compound  = mean(pct_marker_with_compound, na.rm = TRUE),
      bg_pct_nonmarker_with_compound = mean(pct_nonmarker_with_compound, na.rm = TRUE),
      bg_max_pct_compound          = suppressWarnings(
                                       max(pct_compound_all_cells, na.rm = TRUE)),
      .groups = "drop")

  message("\nEmpirical false-positive rate in control animals:")
  print(qc$background)
  message("Treated values at or below these are not detection.\n")

  ## --- 7. Models -------------------------------------------------------------
  power_note <- design_power_note(cfg)
  qc$design_power <- power_note
  message(power_note$note)

  if (isTRUE(cfg$statistics$fit_models) &&
      identical(cfg$design$type, "dose_response")) {

    endpoints <- intersect(c("pct_compound_all_cells", "pct_marker_with_compound",
                             "pct_nonmarker_with_compound",
                             "compound_enrichment_in_marker",
                             "marker_density_per_mm2"), names(region_tbl))

    trend <- dplyr::bind_rows(lapply(endpoints, function(ep) {
      dplyr::bind_rows(lapply(unique(as.character(region_tbl$region)), function(rg) {
        d <- region_tbl[as.character(region_tbl$region) == rg &
                          is.finite(region_tbl[[ep]]), , drop = FALSE]
        if (nrow(d) < 4) return(NULL)
        fit <- try(fit_dose_trend(
          d, value_col = ep, dose_col = "dose",
          covariates = intersect(cfg$statistics$covariates, names(d)),
          random = "(1 | animal_id)",
          log_dose = isTRUE(cfg$statistics$dose_log)), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        tidy_model(fit, exponentiate = FALSE) |>
          dplyr::filter(grepl("dose", term)) |>
          dplyr::mutate(endpoint = ep, region = rg, n_animals = nrow(d))
      }))
    }))

    if (nrow(trend)) {
      trend <- adjust_multiplicity(
        trend, grid_cols = c("endpoint", "region"),
        method = cfg$statistics$multiplicity_method %||% "BH")
      readr::write_csv(trend, file.path(dirs$models, "dose_trend.csv"))
      print(dplyr::select(trend, endpoint, region, n_animals, estimate,
                          conf_low, conf_high, p_value, p_adjusted))
    }
  }

  ## Object level only: the hurdle separates two different claims that a
  ## per-compartment percentage conflates -- whether marker-positive cells are
  ## more likely to contain ANY detectable compound, and whether those that do
  ## contain MORE.
  if (identical(LEVEL, "object") && isTRUE(cfg$statistics$fit_models)) {
    treated <- dplyr::filter(tbl,
                             as.character(group) != cfg$design$control_group)
    compound_col <- names(cols$probes)[1]
    marker_col   <- paste0(names(cols$markers)[1], "_pos")

    if (nrow(treated) > 200) {
      if (any(abs(treated[[compound_col]] - round(treated[[compound_col]])) > 1e-8,
              na.rm = TRUE)) {
        treated[[compound_col]] <- floor(treated[[compound_col]])
      }
      fit <- try(fit_hurdle_nb(
        treated, compound_col, fixed = paste0(marker_col, " + region"),
        random = cfg$statistics$random_effects,
        zi = paste0("~ ", marker_col)), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        readr::write_csv(tidy_model(fit),
                         file.path(dirs$models, "compartment_hurdle.csv"))
        saveRDS(fit, file.path(dirs$models, "compartment_hurdle.rds"))
        icc <- report_icc(fit)
        if (!is.null(icc)) capture.output(print(icc),
          file = file.path(dirs$models, "compartment_icc.txt"))
      }
    }
  }

  ## --- 8. Figures ------------------------------------------------------------
  shape_by <- if ("sex" %in% names(region_tbl)) "sex" else NULL
  fw <- cfg$output$figure_width  %||% 13
  fh <- cfg$output$figure_height %||% 6

  bg_note <- function(col) {
    b <- qc$background[[col]]
    if (is.null(b) || all(is.na(b))) return(NULL)
    sprintf(paste0("Control-animal background: %.2f to %.2f%%. Values in that ",
                   "range are not detection."),
            min(b, na.rm = TRUE), max(b, na.rm = TRUE))
  }

  p1 <- plot_animal_points(
    region_tbl, y_col = "pct_compound_all_cells", x_col = "group",
    facet = "region", shape_by = shape_by,
    y_label = "% cells with detectable compound signal",
    title = "Compound distribution across all detected cells",
    subtitle = paste0("Marker-independent: unaffected by ", marker_name,
                      " staining failures"),
    caption = paste(c(bg_note("bg_pct_compound_all_cells"),
      paste0("Detected compound is not active compound, and no limit of ",
             "detection has been published for oligonucleotide ISH.")),
      collapse = "\n"))
  save_plot(p1, "01_compound_positivity_all_cells.png", dirs$plots,
            width = fw, height = fh)

  p2 <- plot_animal_points(
    dplyr::filter(region_tbl, !marker_failed),
    y_col = "pct_marker_with_compound", x_col = "group", facet = "region",
    shape_by = shape_by,
    y_label = paste0("% ", marker_name, "+ cells containing compound"),
    title = paste0("Compound uptake by ", marker_name, "-positive cells"),
    caption = paste(c(bg_note("bg_pct_marker_with_compound"),
      paste0("Strata where ", marker_name, " failed are omitted from this ",
             "panel AND from every other marker-dependent panel.")),
      collapse = "\n"))
  save_plot(p2, "02_uptake_in_marker_positive.png", dirs$plots,
            width = fw, height = fh)

  if ("compound_enrichment_in_marker" %in% names(region_tbl)) {
    p3 <- plot_animal_points(
      dplyr::filter(region_tbl, !marker_failed),
      y_col = "compound_enrichment_in_marker", x_col = "group",
      facet = "region", shape_by = shape_by,
      y_label = paste0("positive rate in ", marker_name, "+ / in ",
                       marker_name, "-"),
      title = "Preferential occupancy of the marker-positive compartment",
      caption = paste0(
        "A ratio of 1 means the compound is distributed in proportion to cell ",
        "numbers.\nAbove 1 means over-representation in ", marker_name,
        "+ cells; it does not mean more compound per cell.")) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")
    save_plot(p3, "03_compartment_enrichment.png", dirs$plots,
              width = fw, height = fh)
  }

  p4 <- plot_animal_points(
    dplyr::filter(region_tbl, !marker_failed),
    y_col = "marker_density_per_mm2", x_col = "group", facet = "region",
    shape_by = shape_by, y_label = paste0(marker_name, "+ cells / mm2"),
    title = paste0(marker_name, "-positive cell density"),
    subtitle = "A tissue-response readout in its own right, and the antibody QC",
    caption = paste0("Strata where ", marker_name, " yielded zero positive ",
                     "cells are masked, not plotted as zero.\nSee ",
                     "qc/marker_failure.csv."))
  save_plot(p4, "04_marker_density.png", dirs$plots, width = fw, height = fh)

  p5 <- ggplot2::ggplot(qc$marker_failure,
                        ggplot2::aes(x = animal_id, y = marker_total,
                                     fill = marker_failed)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~region, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey45", `TRUE` = "#C1272D"),
                               name = paste0(marker_name, " failed")) +
    ggplot2::labs(x = NULL, y = paste0(marker_name, "+ cells detected"),
                  title = paste0(marker_name, " detection per animal and region"),
      caption = paste0(
        "Red bars are strata where the antibody yielded zero or near-zero ",
        "positive cells. That is a\nstaining or classification failure, not a ",
        "biological zero, and it is masked everywhere downstream.")) +
    theme_rnascope_slide() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  save_plot(p5, "05_marker_detection_qc.png", dirs$plots, width = fw, height = 7)

  p6 <- plot_animal_points(
    region_tbl, y_col = "cells_per_mm2", x_col = "group", facet = "region",
    shape_by = shape_by, y_label = "Detected cells / mm2", title = "Cellularity",
    subtitle = "Reported alongside every per-area readout, always")
  save_plot(p6, "06_cellularity.png", dirs$plots, width = fw, height = fh)

  ## --- 9. Provenance ---------------------------------------------------------
  write_qc_tables(qc, dirs$qc)
  write_run_manifest(cfg, dirs, manifest,
    extra = list(pipeline = "04_compound_distribution_ihc",
                 export_level = LEVEL, marker = marker_name))
  write_methods_stub(cfg, dirs, region_tbl)
  writeLines(evidence_boundaries_ont(),
             file.path(dirs$root, "EVIDENCE_BOUNDARIES.txt"))

  message("\nDone. Results in ", dirs$root)
  message("Read qc/marker_failure.csv and qc/effective_coverage.csv first, ",
          "then EVIDENCE_BOUNDARIES.txt.\n")

  invisible(list(cfg = cfg, dirs = dirs, region_tbl = region_tbl, qc = qc))
}

if (!interactive() && !exists("SOURCED_FOR_INTERACTIVE_USE")) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (!exists("CONFIG_PATH")) {
    CONFIG_PATH <- if (length(.args)) .args[1] else
      stop("Usage: Rscript 04_compound_distribution_ihc.R <config.yml>",
           call. = FALSE)
  }
  RESULT <- main(CONFIG_PATH)
} else if (exists("CONFIG_PATH")) {
  RESULT <- main(CONFIG_PATH)
}
