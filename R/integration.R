
ish_load_integration_config <- function(path) {
  if (!file.exists(path)) stop("Integration config not found.",call.=FALSE)
  cfg <- yaml::read_yaml(path)
  req <- c("schema_version","study","design","panels","integration","output")
  miss <- setdiff(req,names(cfg))
  if (length(miss)) stop("Integration config missing: ",paste(miss,collapse=", "),call.=FALSE)
  if (!identical(as.character(cfg$schema_version),"1.0")) stop("Unsupported schema_version.",call.=FALSE)
  if (length(cfg$panels)<2) stop("At least two panels are required.",call.=FALSE)
  if (isTRUE(cfg$design$panels_independent_replicates))
    stop("Consecutive-section integration requires panels_independent_replicates: false.",call.=FALSE)
  cfg$.config_path <- normalizePath(path,mustWork=TRUE)
  cfg$.config_dir <- dirname(cfg$.config_path)
  cfg$.config_hash <- digest::digest(file=path,algo="sha256")
  cfg$output$dpi <- cfg$output$dpi %||% 300
  cfg
}

ish_int_path <- function(path,cfg) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])",path)) return(normalizePath(path,mustWork=FALSE))
  normalizePath(file.path(cfg$.config_dir,path),mustWork=FALSE)
}

ish_read_integrated_panels <- function(cfg) {
  res <- list(); comp <- list()
  for (nm in names(cfg$panels)) {
    p <- cfg$panels[[nm]]
    f <- ish_int_path(p$standardized_results,cfg)
    if (!file.exists(f)) stop("Panel results not found: ",f,call.=FALSE)
    d <- readr::read_csv(f,show_col_types=FALSE)
    d$panel <- p$display_name %||% nm
    d$panel_key <- nm
    res[[nm]] <- d
    if (!is.null(p$cell_proportions)) {
      cf <- ish_int_path(p$cell_proportions,cfg)
      if (file.exists(cf)) {
        z <- readr::read_csv(cf,show_col_types=FALSE)
        z$panel <- p$display_name %||% nm
        comp[[nm]] <- z
      }
    }
  }
  list(results=dplyr::bind_rows(res),cell_proportions=dplyr::bind_rows(comp))
}

ish_validate_integration <- function(data,cfg) {
  req <- c("animal_id","group","panel","endpoint","metric","raw_value")
  miss <- setdiff(req,names(data))
  if (length(miss)) stop("Standardized tables missing: ",paste(miss,collapse=", "),call.=FALSE)
  roster <- dplyr::bind_rows(lapply(cfg$design$animals,tibble::as_tibble)) |>
    dplyr::mutate(animal_id=as.character(animal_id),expected_group=as.character(group)) |>
    dplyr::select(animal_id,expected_group)
  chk <- data |> dplyr::distinct(animal_id,group) |> dplyr::left_join(roster,by="animal_id")
  if (any(is.na(chk$expected_group)) || any(as.character(chk$group)!=chk$expected_group))
    stop("Animal/treatment metadata are inconsistent across panels.",call.=FALSE)
  invisible(TRUE)
}

ish_integrated_global <- function(data) {
  grp <- intersect(c("study_id","animal_id","group","panel","panel_key","cell_type","endpoint",
                     "endpoint_kind","endpoint_role","unit","metric"),names(data))
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(raw_value=ish_safe_mean(raw_value),
      normalized_value=ish_safe_mean(normalized_value),
      n_cells=sum(n_cells,na.rm=TRUE),
      n_regions=if ("region" %in% names(data)) dplyr::n_distinct(region) else NA_integer_,
      .groups="drop") |>
    dplyr::mutate(region="All analyzed regions (equal-region mean)")
}

ish_run_integration <- function(path) {
  cfg <- ish_load_integration_config(path)
  root <- ish_int_path(cfg$output$root,cfg)
  dirs <- list(root=root,tables=file.path(root,"tables"),qc=file.path(root,"qc"),
    scientific=file.path(root,"figures","FULL_SCIENTIFIC"),
    presentation=file.path(root,"figures","PRESENTATION_COMPACT"),
    run=file.path(root,"run_manifest"))
  for (d in dirs) dir.create(d,recursive=TRUE,showWarnings=FALSE)

  inp <- ish_read_integrated_panels(cfg)
  data <- inp$results
  ish_validate_integration(data,cfg)
  global <- ish_integrated_global(data)

  readr::write_csv(data,file.path(dirs$tables,"integrated_standardized_results.csv"))
  readr::write_csv(global,file.path(dirs$tables,"integrated_global_equal_region.csv"))

  if (nrow(inp$cell_proportions)) {
    readr::write_csv(inp$cell_proportions,file.path(dirs$qc,"integrated_cell_proportions.csv"))
    ish_save_plot(ish_plot_cell_proportions(inp$cell_proportions),
      file.path(dirs$scientific,"QC_cell_proportions_by_animal.png"),12,6,cfg$output$dpi)
  }

  cc <- list()
  for (nm in names(cfg$integration$shared_measurements %||% list())) {
    s <- cfg$integration$shared_measurements[[nm]]
    cc[[nm]] <- ish_panel_concordance(data,s$endpoint,s$metric,s$cell_type %||% "All cells")
  }
  concord <- dplyr::bind_rows(cc)
  if (nrow(concord)) readr::write_csv(concord,file.path(dirs$qc,"panel_concordance.csv"))

  for (nm in names(cfg$integration$main_figures %||% list())) {
    f <- cfg$integration$main_figures[[nm]]
    source <- if (f$scope %in% c("global_shared","celltype")) global else data
    d <- source |> dplyr::filter(endpoint==f$endpoint,metric==f$metric)
    if (!is.null(f$cell_type)) d <- dplyr::filter(d,cell_type==f$cell_type)
    if (f$scope=="celltype") d <- dplyr::filter(d,cell_type!="All cells")
    d <- d |> dplyr::mutate(plot_value=dplyr::coalesce(normalized_value,raw_value))
    if (!nrow(d)) next
    if (f$scope=="global_shared") {
      ps <- ish_plot_panel_pair(d,"plot_value",f$y_label,f$title,f$reference_line %||% NULL,FALSE)
      pp <- ish_plot_panel_pair(d,"plot_value",f$y_label,NULL,f$reference_line %||% NULL,TRUE)
    } else {
      xcol <- if (f$scope=="regional") "region" else "cell_type"
      ps <- ish_plot_panel_facets(d,"plot_value",xcol,f$y_label,f$title,f$reference_line %||% NULL,FALSE)
      pp <- ish_plot_panel_facets(d,"plot_value",xcol,f$y_label,NULL,f$reference_line %||% NULL,TRUE)
    }
    ish_save_plot(ps,file.path(dirs$scientific,paste0(nm,".png")),10,5.5,cfg$output$dpi)
    ish_save_plot(pp,file.path(dirs$presentation,paste0(nm,".png")),8.5,4.8,cfg$output$dpi)
  }

  mrna <- global |> dplyr::filter(endpoint_role=="target_mrna",metric=="mean_all",cell_type!="All cells")
  if (nrow(mrna)) {
    p <- ish_presentation_mrna(mrna,cfg$design$control_group)
    ish_save_plot(p,file.path(dirs$presentation,"target_mRNA_remaining_by_celltype_0_100.png"),
                  9,4.8,cfg$output$dpi)
    flagged <- mrna |> dplyr::filter(is.finite(normalized_value),
                                     normalized_value<0|normalized_value>100)
    readr::write_csv(flagged,file.path(dirs$qc,"presentation_values_outside_0_100.csv"))
  }

  n <- data |> dplyr::distinct(animal_id,group) |> dplyr::count(group,name="n_animals")
  readr::write_csv(n,file.path(dirs$qc,"biological_n.csv"))
  q <- dplyr::bind_rows(
    ish_qc_row("panel_replication_rule","PASS",
      "Panels are repeated consecutive-section measurements and do not increase biological n."),
    ish_qc_row("animal_count","PASS",paste0(dplyr::n_distinct(data$animal_id)," unique animals")))
  readr::write_csv(q,file.path(dirs$qc,"QC_STATUS.csv"))

  yaml::write_yaml(list(schema_version=cfg$schema_version,timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
    config_sha256=cfg$.config_hash,git_commit=ish_git_commit(),
    panel_relationship=cfg$integration$panel_relationship,
    biological_unit=cfg$design$biological_unit,
    panels_independent_replicates=cfg$design$panels_independent_replicates),
    file.path(dirs$run,"run_manifest.yml"))
  file.copy(cfg$.config_path,file.path(dirs$run,basename(cfg$.config_path)),overwrite=TRUE)
  ish_write_output_manifest(root)
  invisible(list(cfg=cfg,dirs=dirs,data=data,global=global,concordance=concord,qc=q))
}
