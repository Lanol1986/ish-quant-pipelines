
ish_qc_row <- function(check,status,message,panel=NA_character_,animal_id=NA_character_,
                       metric=NA_real_,threshold=NA_real_) {
  tibble::tibble(check=check,status=status,message=message,panel=panel,
                 animal_id=animal_id,metric=metric,threshold=threshold)
}

ish_qc_status <- function(data,cfg,manifest) {
  roster_ok <- setequal(unique(manifest$animal_id),
                        setdiff(as.character(cfg$design$animals_tbl$animal_id),
                                names(cfg$design$expected_absent %||% list())))
  q <- list(
    roster=ish_qc_row("animal_roster",if (roster_ok) "PASS" else "FAIL",
                      paste0("files=",length(unique(manifest$animal_id))),cfg$panel$name)
  )
  for (ep in cfg$endpoints) {
    x <- data[[ep$name]]
    bad <- sum(is.finite(x) & x<0)
    q[[paste0("ep_",ep$name)]] <- ish_qc_row(paste0("endpoint_validity__",ep$name),
      if (bad) "REVIEW" else "PASS",paste0("negative values=",bad),cfg$panel$name)
  }
  if (length(cfg$phenotypes)) {
    amb <- 100*mean(data$phenotype_memberships>1,na.rm=TRUE)
    alarm <- cfg$qc$ambiguous_phenotype_alarm_pct %||% 5
    q$phenotype <- ish_qc_row("phenotype_partition",if (amb>alarm) "REVIEW" else "PASS",
                              sprintf("ambiguous=%.2f%%",amb),cfg$panel$name,metric=amb,threshold=alarm)
  }
  if (isTRUE(cfg$spatial$enabled)) {
    cols <- paste0("spatial__",c("xmin","xmax","ymin","ymax"))
    ok <- all(cols %in% names(data))
    bad <- if (ok) sum(!stats::complete.cases(data[,cols])) else nrow(data)
    q$coords <- ish_qc_row("coordinates",if (!ok) "FAIL" else if (bad) "REVIEW" else "PASS",
                           paste0("rows with incomplete coordinates=",bad),cfg$panel$name,metric=bad)
  }
  dplyr::bind_rows(q)
}

ish_cell_proportions <- function(data) {
  by <- intersect(c("panel","animal_id","group",ish_anatomy_columns(data)),names(data))
  data |>
    dplyr::filter(!is.na(cell_type),cell_type!="All cells") |>
    dplyr::count(dplyr::across(dplyr::all_of(c(by,"cell_type"))),name="n_cells") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::mutate(pct_cells=100*n_cells/sum(n_cells)) |>
    dplyr::ungroup()
}

ish_marker_separation <- function(data,cfg) {
  out <- list()
  for (nm in names(cfg$markers)) {
    pc <- paste0(nm,"_positive"); ic <- paste0(nm,"_intensity")
    if (!all(c(pc,ic) %in% names(data))) next
    out[[nm]] <- data |>
      dplyr::group_by(panel,animal_id,group) |>
      dplyr::summarise(pos_median=ish_safe_median(.data[[ic]][.data[[pc]] %in% TRUE]),
                       neg_median=ish_safe_median(.data[[ic]][.data[[pc]] %in% FALSE]),
                       ratio=pos_median/neg_median,.groups="drop") |>
      dplyr::mutate(marker=nm)
  }
  dplyr::bind_rows(out)
}

ish_threshold_sweep <- function(data,cfg,endpoint_name,thresholds=NULL) {
  ep <- ish_endpoint_spec(cfg,endpoint_name)
  base <- ep$positivity_threshold
  if (is.null(base)) return(tibble::tibble())
  thresholds <- thresholds %||% unique(pmax(0,base+c(-1,0,1,2)))
  by <- intersect(c("panel","animal_id","group",ish_anatomy_columns(data),"cell_type"),names(data))
  purrr::map_dfr(thresholds,function(k) {
    data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
      dplyr::summarise(threshold=k,n_cells=dplyr::n(),
                       pct_positive=100*mean(.data[[endpoint_name]]>=k,na.rm=TRUE),
                       .groups="drop")
  }) |> dplyr::mutate(endpoint=endpoint_name)
}

ish_panel_concordance <- function(data,endpoint,metric,cell_type="All cells") {
  d <- data |> dplyr::filter(.data$endpoint==endpoint,.data$metric==metric,
                             .data$cell_type==cell_type) |>
    dplyr::mutate(value=dplyr::coalesce(normalized_value,raw_value)) |>
    dplyr::select(animal_id,panel,value) |> dplyr::distinct()
  panels <- unique(as.character(d$panel))
  if (length(panels)!=2) return(tibble::tibble())
  w <- tidyr::pivot_wider(d,names_from=panel,values_from=value)
  ok <- stats::complete.cases(w[,panels,drop=FALSE])
  tibble::tibble(endpoint=endpoint,metric=metric,cell_type=cell_type,
                 panel_1=panels[1],panel_2=panels[2],n_pairs=sum(ok),
                 pearson_r=if (sum(ok)>=3) stats::cor(w[[panels[1]]][ok],w[[panels[2]]][ok]) else NA_real_,
                 spearman_r=if (sum(ok)>=3) stats::cor(w[[panels[1]]][ok],w[[panels[2]]][ok],method="spearman") else NA_real_)
}


ish_specimen_qc_metrics <- function(data,cfg) {
  cols <- grep("^qc_metric__",names(data),value=TRUE)
  if (!length(cols)) return(tibble::tibble())
  data |>
    dplyr::group_by(panel,animal_id,group) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cols),
      list(median=~ish_safe_median(.x),
           p10=~if (sum(is.finite(.x))) stats::quantile(.x,.10,na.rm=TRUE,names=FALSE) else NA_real_,
           p90=~if (sum(is.finite(.x))) stats::quantile(.x,.90,na.rm=TRUE,names=FALSE) else NA_real_),
      .names="{.col}__{.fn}"),.groups="drop")
}

ish_vehicle_background <- function(std,cfg) {
  ctrl <- cfg$design$control_group
  std |>
    dplyr::filter(as.character(group)==ctrl,
                  endpoint_role %in% c("sirna","compound"),
                  metric=="pct_positive") |>
    dplyr::select(dplyr::any_of(c("panel","animal_id","group","region","subregion","layer",
      "cell_type","endpoint","metric","raw_value","n_cells")))
}

ish_vehicle_target_baseline <- function(std,cfg) {
  ctrl <- cfg$design$control_group
  std |>
    dplyr::filter(as.character(group)==ctrl,
                  endpoint_role=="target_mrna",
                  metric=="mean_all") |>
    dplyr::select(dplyr::any_of(c("panel","animal_id","group","region","subregion","layer",
      "cell_type","endpoint","metric","raw_value","n_cells")))
}
