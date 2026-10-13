
ish_run_single_panel <- function(config_path,assay_mode="multiplex") {
  cfg <- ish_load_config(config_path)
  dirs <- ish_output_dirs(cfg)
  rr <- ish_read_panel(cfg)

  data <- rr$data |>
    dplyr::mutate(study_id=cfg$study$id %||% cfg$study$name) |>
    ish_apply_anatomy(cfg) |>
    ish_assign_phenotypes(cfg)

  analysis_data <- ish_add_all_cells(data)

  summaries <- dplyr::bind_rows(lapply(vapply(cfg$endpoints,\`[[\`,character(1),"name"),
    function(ep) ish_summarise_endpoint(analysis_data,cfg,ep)))
  std <- ish_standardize(summaries) |> ish_normalize_target(cfg)
  global <- ish_global_equal_region(std)

  readr::write_csv(summaries,file.path(dirs$tables,"animal_endpoint_summary.csv"))
  readr::write_csv(std,file.path(dirs$tables,"standardized_results.csv"))
  readr::write_csv(global,file.path(dirs$tables,"standardized_global_equal_region.csv"))
  readr::write_csv(rr$manifest,file.path(dirs$qc,"file_manifest.csv"))
  readr::write_csv(rr$qc,file.path(dirs$qc,"read_qc.csv"))

  comp <- ish_cell_proportions(analysis_data)
  if (nrow(comp)) {
    readr::write_csv(comp,file.path(dirs$qc,"cell_proportions.csv"))
    ish_save_plot(ish_plot_cell_proportions(comp),file.path(dirs$scientific,"QC_cell_proportions.png"),
                  11,6,cfg$output$dpi)
  }
  marker_qc <- ish_marker_separation(data,cfg)
  if (nrow(marker_qc)) readr::write_csv(marker_qc,file.path(dirs$qc,"marker_separation.csv"))

  qstat <- ish_qc_status(data,cfg,rr$manifest)
  readr::write_csv(qstat,file.path(dirs$qc,"QC_STATUS.csv"))

  for (ep in cfg$endpoints) {
    if (!is.null(ep$positivity_threshold)) {
      sw <- ish_threshold_sweep(analysis_data,cfg,ep$name)
      readr::write_csv(sw,file.path(dirs$qc,paste0("threshold_sweep__",ep$name,".csv")))
    }
  }

  if (isTRUE(cfg$spatial$enabled)) {
    meta <- list()
    for (ep in cfg$endpoints) {
      if (!((ep$role %||% "other") %in% c("sirna","compound","target_mrna"))) next
      for (a in unique(data$animal_id)) {
        d <- dplyr::filter(data,animal_id==a)
        sp <- ish_spatial_plot(d,ep$name,cfg,paste(cfg$panel$name,a,ep$name))
        ish_save_plot(sp$plot,file.path(dirs$spatial,paste0("spatial__",ep$name,"__",a,".png")),
                      8,7,cfg$output$dpi)
        meta[[paste(ep$name,a)]] <- dplyr::mutate(sp$metadata,animal_id=a,panel=cfg$panel$name)
      }
    }
    if (length(meta)) readr::write_csv(dplyr::bind_rows(meta),file.path(dirs$qc,"spatial_scale_metadata.csv"))
  }

  run_record <- list(assay_mode=assay_mode,
    analyses_executed=list(animal_level_descriptive=TRUE,
      control_normalization=any(vapply(cfg$endpoints,function(x) identical(x$role,"target_mrna"),logical(1))),
      automatic_count_model=FALSE,spatial=isTRUE(cfg$spatial$enabled)),
    note="Count models are never enabled by rounding/flooring continuous estimated indices.")
  ish_write_manifest(cfg,dirs,rr$manifest,run_record)
  ish_write_output_manifest(dirs$root)

  invisible(list(cfg=cfg,dirs=dirs,data=data,summary=summaries,standardized=std,qc=qstat))
}
