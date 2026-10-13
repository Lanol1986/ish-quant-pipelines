
ish_as_positive_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x > 0)
  z <- tolower(trimws(as.character(x)))
  yes <- c("1","true","positive","pos","yes","y")
  no  <- c("0","false","negative","neg","no","n","")
  bad <- setdiff(unique(z), c(yes,no,NA_character_))
  if (length(bad)) stop("Unrecognized positivity values: ", paste(bad, collapse=", "), call.=FALSE)
  ifelse(is.na(z), NA, z %in% yes)
}

ish_apply_missing_policy <- function(x, policy, endpoint) {
  raw <- x
  val <- suppressWarnings(as.numeric(raw))
  nonempty_bad <- is.na(val) & !is.na(raw) & trimws(as.character(raw)) != ""
  if (any(nonempty_bad)) stop(endpoint, ": non-empty values could not be parsed as numeric.", call.=FALSE)
  n_missing <- sum(is.na(val))
  if (policy == "error" && n_missing) stop(endpoint, ": missing values found with missing_policy='error'.", call.=FALSE)
  if (policy == "zero_if_export_blank") val[is.na(val)] <- 0
  list(values=val, n_missing_input=n_missing, policy=policy)
}

ish_build_manifest <- function(cfg) {
  input_dir <- ish_resolve_path(cfg$panel$input_dir, cfg)
  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir, call.=FALSE)
  files <- list.files(input_dir, pattern=cfg$panel$file_pattern, full.names=TRUE)
  files <- files[!grepl("^\\._", basename(files))]
  if (!length(files)) stop("No files matched panel.file_pattern.", call.=FALSE)
  ids <- stringr::str_match(basename(files), cfg$panel$id_regex)[,2]
  if (anyNA(ids)) stop("Could not extract animal ID from all filenames.", call.=FALSE)
  tibble::tibble(panel=cfg$panel$name, source_file=basename(files), path=files,
                 animal_id=as.character(ids), size_bytes=file.info(files)$size)
}

ish_assert_roster <- function(manifest, cfg) {
  expected <- as.character(cfg$design$animals_tbl$animal_id)
  absent <- cfg$design$expected_absent %||% list()
  absent_ids <- if (is.list(absent)) names(absent) else as.character(absent)
  found <- unique(manifest$animal_id)
  missing <- setdiff(expected, c(found, absent_ids))
  extra <- setdiff(found, expected)
  if (length(missing)) stop("Missing declared animals: ", paste(missing,collapse=", "), call.=FALSE)
  if (length(extra)) stop("Undeclared animals found: ", paste(extra,collapse=", "), call.=FALSE)
  invisible(TRUE)
}

ish_read_delimited_strict <- function(path, delim=",") {
  x <- readr::read_delim(path, delim=delim, show_col_types=FALSE, progress=FALSE,
                         name_repair="minimal", guess_max=100000)
  probs <- readr::problems(x)
  if (nrow(probs)) stop(nrow(probs), " parsing problem(s) in ", basename(path), call.=FALSE)
  x
}

ish_required_source_columns <- function(cfg) {
  ep <- vapply(cfg$endpoints, function(x) x$source_column, character(1))
  ph <- unlist(lapply(cfg$phenotypes, function(x) x$source_column), use.names=FALSE)
  mk <- unlist(lapply(cfg$markers, function(x) c(x$positivity_column,x$intensity_column)), use.names=FALSE)
  qm <- unlist(lapply(cfg$qc_metrics %||% list(), function(x) x$source_column), use.names=FALSE)
  an <- unlist(lapply(cfg$anatomy$fields %||% list(), function(x) x$source_column), use.names=FALSE)
  sp <- unlist(cfg$spatial$columns %||% list(), use.names=FALSE)
  unique(stats::na.omit(c(ep,ph,mk,qm,an,sp)))
}

ish_read_one_export <- function(path, cfg) {
  df <- ish_read_delimited_strict(path, cfg$panel$delimiter)
  req <- ish_required_source_columns(cfg)
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing required columns in ", basename(path), ": ", paste(miss,collapse=", "), call.=FALSE)

  out <- tibble::tibble(.row_id=seq_len(nrow(df)))
  qc <- list()
  for (ep in cfg$endpoints) {
    z <- ish_apply_missing_policy(df[[ep$source_column]], ep$missing_policy %||% "error", ep$name)
    out[[ep$name]] <- z$values
    qc[[ep$name]] <- tibble::tibble(endpoint=ep$name,n_rows=nrow(df),
      n_missing_input=z$n_missing_input, missing_policy=z$policy)
  }
  for (nm in names(cfg$phenotypes)) {
    out[[paste0("phenotype__",nm)]] <- ish_as_positive_flag(df[[cfg$phenotypes[[nm]]$source_column]])
  }
  for (nm in names(cfg$markers)) {
    m <- cfg$markers[[nm]]
    if (!is.null(m$positivity_column)) out[[paste0(nm,"_positive")]] <- ish_as_positive_flag(df[[m$positivity_column]])
    if (!is.null(m$intensity_column)) out[[paste0(nm,"_intensity")]] <- suppressWarnings(as.numeric(df[[m$intensity_column]]))
  }
  for (nm in names(cfg$qc_metrics %||% list())) {
    q <- cfg$qc_metrics[[nm]]
    out[[paste0("qc_metric__",nm)]] <- suppressWarnings(as.numeric(df[[q$source_column]]))
  }
  for (nm in names(cfg$anatomy$fields %||% list())) {
    out[[paste0("anatomy_raw__",nm)]] <- as.character(df[[cfg$anatomy$fields[[nm]]$source_column]])
  }
  for (nm in names(cfg$spatial$columns %||% list())) {
    out[[paste0("spatial__",nm)]] <- suppressWarnings(as.numeric(df[[cfg$spatial$columns[[nm]]]]))
  }
  list(data=out, qc=dplyr::bind_rows(qc))
}

ish_read_panel <- function(cfg) {
  manifest <- ish_build_manifest(cfg)
  ish_assert_roster(manifest,cfg)
  design <- cfg$design$animals_tbl
  items <- lapply(seq_len(nrow(manifest)), function(i) {
    m <- manifest[i,]
    z <- ish_read_one_export(m$path,cfg)
    d <- z$data |>
      dplyr::mutate(panel=cfg$panel$name, animal_id=m$animal_id,
                    source_file=m$source_file,
                    section_id=paste(cfg$panel$name,m$animal_id,m$source_file,sep=":")) |>
      dplyr::left_join(design,by="animal_id")
    q <- z$qc |> dplyr::mutate(panel=cfg$panel$name,animal_id=m$animal_id,source_file=m$source_file)
    list(data=d,qc=q)
  })
  list(data=dplyr::bind_rows(lapply(items,\`[[\`,"data")),
       qc=dplyr::bind_rows(lapply(items,\`[[\`,"qc")),
       manifest=manifest)
}
