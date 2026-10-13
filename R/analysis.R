
ish_norm_annotation <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("[^a-z0-9]+","_",z)
  gsub("^_|_$","",gsub("_+","_",z))
}

ish_apply_anatomy <- function(data,cfg) {
  out <- data
  for (nm in names(cfg$anatomy$fields %||% list())) {
    spec <- cfg$anatomy$fields[[nm]]
    raw <- out[[paste0("anatomy_raw__",nm)]]
    if (is.null(spec$levels)) { out[[nm]] <- raw; next }
    rn <- ish_norm_annotation(raw)
    val <- rep(NA_character_,length(raw))
    aliases <- spec$aliases %||% list()
    for (level in spec$levels) {
      pats <- ish_norm_annotation(c(level, aliases[[level]] %||% character(0)))
      for (p in pats) {
        hit <- is.na(val) & grepl(p,rn,fixed=TRUE)
        val[hit] <- level
      }
    }
    out[[nm]] <- val
  }
  out
}

ish_assign_phenotypes <- function(data,cfg) {
  if (!length(cfg$phenotypes)) {
    data$cell_type <- "All cells"
    data$phenotype_memberships <- 1L
    return(data)
  }
  mat <- matrix(FALSE,nrow(data),length(cfg$phenotypes),dimnames=list(NULL,names(cfg$phenotypes)))
  for (nm in names(cfg$phenotypes)) {
    rule <- cfg$phenotypes[[nm]]
    if (!is.null(rule$source_column)) {
      mat[,nm] <- data[[paste0("phenotype__",nm)]] %in% TRUE
    } else {
      ans <- rep(TRUE,nrow(data))
      pa <- rule$positive_all %||% rule$positive %||% character(0)
      py <- rule$positive_any %||% character(0)
      ng <- rule$negative %||% character(0)
      if (length(pa)) ans <- ans & apply(data[paste0(pa,"_positive")],1,function(z) all(z %in% TRUE))
      if (length(py)) ans <- ans & apply(data[paste0(py,"_positive")],1,function(z) any(z %in% TRUE))
      if (length(ng)) ans <- ans & apply(data[paste0(ng,"_positive")],1,function(z) all(!(z %in% TRUE)))
      mat[,nm] <- ans
    }
  }
  n <- rowSums(mat)
  cell_type <- rep("Unclassified",nrow(data))
  one <- n==1
  if (any(one)) cell_type[one] <- colnames(mat)[max.col(mat[one,,drop=FALSE],ties.method="first")]
  cell_type[n>1] <- "Ambiguous"
  if ((cfg$phenotype_policy %||% "report_ambiguous")=="exclude_ambiguous") cell_type[n>1] <- NA_character_
  data$cell_type <- cell_type
  data$phenotype_memberships <- n
  data
}

ish_add_all_cells <- function(data) {
  all <- data; all$cell_type <- "All cells"
  dplyr::bind_rows(data,all)
}

ish_anatomy_columns <- function(data) intersect(c("organ","region","subregion","layer"),names(data))

ish_summarise_endpoint <- function(data,cfg,endpoint_name) {
  ep <- ish_endpoint_spec(cfg,endpoint_name)
  by <- intersect(c("study_id","animal_id","group","panel","section_id",
                    ish_anatomy_columns(data),"cell_type"),names(data))
  threshold <- ep$positivity_threshold %||% NA_real_
  d <- data |> dplyr::mutate(.value=.data[[endpoint_name]],
                             .positive=if (is.finite(threshold)) .value>=threshold else NA)
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(n_cells=dplyr::n(), n_measured=sum(is.finite(.value)),
      n_positive=if (is.finite(threshold)) sum(.positive,na.rm=TRUE) else NA_integer_,
      pct_positive=if (is.finite(threshold)) 100*mean(.positive,na.rm=TRUE) else NA_real_,
      mean_all=ish_safe_mean(.value), median_all=ish_safe_median(.value),
      mean_positive=if (is.finite(threshold)) ish_safe_mean(.value[.positive]) else NA_real_,
      median_positive=if (is.finite(threshold)) ish_safe_median(.value[.positive]) else NA_real_,
      .groups="drop") |>
    dplyr::mutate(endpoint=endpoint_name,endpoint_kind=ep$kind,
                  endpoint_role=ep$role %||% "other",unit=ep$unit %||% "a.u.",
                  positivity_threshold=threshold)
}

ish_standardize <- function(summary_tbl) {
  ids <- intersect(c("study_id","animal_id","group","panel","section_id","organ","region","subregion","layer",
                     "cell_type","endpoint","endpoint_kind","endpoint_role","unit","n_cells","n_measured",
                     "positivity_threshold"),names(summary_tbl))
  metrics <- intersect(c("pct_positive","mean_all","median_all","mean_positive","median_positive"),names(summary_tbl))
  tidyr::pivot_longer(summary_tbl,cols=dplyr::all_of(metrics),names_to="metric",values_to="raw_value") |>
    dplyr::select(dplyr::all_of(ids),metric,raw_value,dplyr::everything()) |>
    dplyr::mutate(normalized_value=NA_real_,normalization_reference=NA_character_)
}

ish_normalize_target <- function(std,cfg) {
  ctrl <- cfg$design$control_group
  out <- std |> dplyr::mutate(.rowid=dplyr::row_number())
  d <- out |> dplyr::filter(endpoint_role=="target_mrna",metric=="mean_all")
  if (!nrow(d)) return(dplyr::select(out,-.rowid))
  strata <- intersect(c("panel","organ","region","subregion","layer","cell_type","endpoint","metric"),names(d))
  ref <- d |> dplyr::filter(as.character(group)==ctrl) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata))) |>
    dplyr::summarise(control_mean=ish_safe_mean(raw_value),.groups="drop")
  z <- dplyr::left_join(d,ref,by=strata) |>
    dplyr::mutate(.norm=100*raw_value/control_mean,
                  .ref=paste0("matched ",ctrl," mean = 100%")) |>
    dplyr::select(.rowid,.norm,.ref)
  dplyr::left_join(out,z,by=".rowid") |>
    dplyr::mutate(normalized_value=dplyr::coalesce(.norm,normalized_value),
                  normalization_reference=dplyr::coalesce(.ref,normalization_reference)) |>
    dplyr::select(-.rowid,-.norm,-.ref)
}

ish_global_equal_region <- function(std) {
  if (!"region" %in% names(std)) return(std)
  grp <- intersect(c("study_id","animal_id","group","panel","cell_type","endpoint","endpoint_kind",
                     "endpoint_role","unit","metric"),names(std))
  std |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(raw_value=ish_safe_mean(raw_value),
      normalized_value=ish_safe_mean(normalized_value),
      n_cells=sum(n_cells,na.rm=TRUE), n_regions=dplyr::n_distinct(region),.groups="drop") |>
    dplyr::mutate(region="All analyzed regions (equal-region mean)",
                  normalization_reference=ifelse(is.finite(normalized_value),
                    "equal-region mean of stratum-normalized values",NA_character_))
}

ish_assert_integer_endpoint <- function(x,name) {
  if (any(is.finite(x) & abs(x-round(x))>1e-8))
    stop("Endpoint '",name,"' is not integer-valued. Do not round/floor an estimated index to run a count model.",call.=FALSE)
  invisible(TRUE)
}

ish_fit_two_part_count_model <- function(data,endpoint,threshold=1,fixed="group",
                                         random="(1 | animal_id/section_id)") {
  if (!requireNamespace("glmmTMB",quietly=TRUE))
    return(list(status="NOT_RUN",reason="glmmTMB not installed"))
  ish_assert_integer_endpoint(data[[endpoint]],endpoint)
  d <- data |> dplyr::mutate(.positive=.data[[endpoint]]>=threshold)
  rhs <- paste(fixed,random,sep=" + ")
  occ <- glmmTMB::glmmTMB(stats::as.formula(paste(".positive ~",rhs)),
                          family=stats::binomial(),data=d)
  pos <- d |> dplyr::filter(.positive)
  burden <- glmmTMB::glmmTMB(stats::as.formula(paste(endpoint,"~",rhs)),
                             family=glmmTMB::truncated_nbinom2(),data=pos)
  list(status="PASS",occurrence=occ,positive_count=burden)
}
