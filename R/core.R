
\`%||%\` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

ish_required_packages <- function() c("yaml","digest","dplyr","readr","tidyr","purrr","tibble","ggplot2","scales")

ish_check_packages <- function() {
  miss <- ish_required_packages()[!vapply(ish_required_packages(), requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required packages: ", paste(miss, collapse=", "), call.=FALSE)
  invisible(TRUE)
}

ish_validate_enum <- function(x, allowed, name) {
  if (is.null(x) || !length(x)) stop(name, " is required.", call.=FALSE)
  if (!all(x %in% allowed)) stop(name, " must be one of: ", paste(allowed, collapse=", "), call.=FALSE)
  invisible(TRUE)
}

ish_load_config <- function(path) {
  ish_check_packages()
  if (!file.exists(path)) stop("Config not found: ", path, call.=FALSE)
  cfg <- yaml::read_yaml(path)
  req <- c("schema_version","study","paths","design","panel","anatomy","endpoints","output")
  miss <- setdiff(req, names(cfg))
  if (length(miss)) stop("Config missing sections: ", paste(miss, collapse=", "), call.=FALSE)
  if (!identical(as.character(cfg$schema_version), "1.0")) stop("Unsupported schema_version.", call.=FALSE)

  animals <- dplyr::bind_rows(lapply(cfg$design$animals, tibble::as_tibble))
  if (!all(c("animal_id","group") %in% names(animals))) stop("design.animals requires animal_id and group.", call.=FALSE)
  if (anyDuplicated(as.character(animals$animal_id))) stop("Duplicate animal_id in design.animals.", call.=FALSE)

  cfg$design$group_levels <- cfg$design$group_levels %||% unique(as.character(animals$group))
  cfg$design$control_group <- cfg$design$control_group %||% cfg$design$group_levels[[1]]
  cfg$design$animals_tbl <- animals |>
    dplyr::mutate(animal_id=as.character(animal_id),
                  group=factor(as.character(group), levels=cfg$design$group_levels))

  cfg$panel$input_dir <- cfg$panel$input_dir %||% cfg$paths$input_dir
  cfg$panel$file_pattern <- cfg$panel$file_pattern %||% "\\.csv$"
  cfg$panel$id_regex <- cfg$panel$id_regex %||% cfg$paths$id_regex %||% "^([A-Za-z0-9.]+?)_"
  cfg$panel$delimiter <- cfg$panel$delimiter %||% ","
  cfg$phenotypes <- cfg$phenotypes %||% list()
  cfg$markers <- cfg$markers %||% list()
  cfg$qc <- cfg$qc %||% list()
  cfg$spatial <- cfg$spatial %||% list(enabled=FALSE)
  cfg$presentation <- cfg$presentation %||% list()
  cfg$output$dpi <- cfg$output$dpi %||% 300

  allowed_kind <- c("integer_count","estimated_spot_index","area_fraction","integrated_intensity","continuous","proportion")
  allowed_role <- c("compound","sirna","target_mrna","marker","other")
  allowed_missing <- c("error","preserve_na","zero_if_export_blank")
  for (ep in cfg$endpoints) {
    if (is.null(ep$name) || is.null(ep$source_column)) stop("Every endpoint requires name and source_column.", call.=FALSE)
    ish_validate_enum(ep$kind, allowed_kind, paste0(ep$name,".kind"))
    ish_validate_enum(ep$role %||% "other", allowed_role, paste0(ep$name,".role"))
    ish_validate_enum(ep$missing_policy %||% "error", allowed_missing, paste0(ep$name,".missing_policy"))
  }

  cfg$.config_path <- normalizePath(path, mustWork=TRUE)
  cfg$.config_dir <- dirname(cfg$.config_path)
  cfg$.config_hash <- digest::digest(file=path, algo="sha256")
  cfg
}

ish_resolve_path <- function(path, cfg) {
  path <- path.expand(path)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(normalizePath(path, mustWork=FALSE))
  normalizePath(file.path(cfg$.config_dir, path), mustWork=FALSE)
}

ish_output_dirs <- function(cfg) {
  root <- ish_resolve_path(cfg$paths$output_dir, cfg)
  dirs <- list(root=root, tables=file.path(root,"tables"), qc=file.path(root,"qc"),
               models=file.path(root,"models"),
               scientific=file.path(root,"figures","FULL_SCIENTIFIC"),
               presentation=file.path(root,"figures","PRESENTATION_COMPACT"),
               spatial=file.path(root,"figures","SPATIAL"),
               run=file.path(root,"run_manifest"))
  for (d in dirs) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  dirs
}

ish_endpoint_spec <- function(cfg, name) {
  i <- which(vapply(cfg$endpoints, function(x) identical(x$name,name), logical(1)))
  if (!length(i)) stop("Unknown endpoint: ", name, call.=FALSE)
  cfg$endpoints[[i[1]]]
}

ish_safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
ish_safe_median <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else stats::median(x) }
