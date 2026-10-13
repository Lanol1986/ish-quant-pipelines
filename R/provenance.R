
ish_sha256 <- function(path) digest::digest(file=path,algo="sha256")

ish_repo_root <- function(start=getwd()) {
  p <- normalizePath(start,mustWork=FALSE)
  repeat {
    if (dir.exists(file.path(p,".git"))) return(p)
    q <- dirname(p); if (identical(q,p)) return(NA_character_); p <- q
  }
}

ish_git_commit <- function(root=ish_repo_root()) {
  if (is.na(root)) return(NA_character_)
  z <- suppressWarnings(try(system2("git",c("-C",root,"rev-parse","HEAD"),stdout=TRUE,stderr=FALSE),silent=TRUE))
  if (inherits(z,"try-error")||!length(z)) return(NA_character_)
  z[1]
}

ish_write_manifest <- function(cfg,dirs,manifest,run_record=list()) {
  hashes <- manifest |> dplyr::mutate(sha256=vapply(path,ish_sha256,character(1))) |>
    dplyr::select(panel,animal_id,source_file,size_bytes,sha256)
  readr::write_csv(hashes,file.path(dirs$run,"input_file_hashes.csv"))
  file.copy(cfg$.config_path,file.path(dirs$run,basename(cfg$.config_path)),overwrite=TRUE)
  writeLines(capture.output(utils::sessionInfo()),file.path(dirs$run,"sessionInfo.txt"))
  yaml::write_yaml(list(schema_version=cfg$schema_version,study=cfg$study,panel=cfg$panel$name,
    timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),git_commit=ish_git_commit(),
    config_sha256=cfg$.config_hash,run_record=run_record),
    file.path(dirs$run,"run_manifest.yml"))
}

ish_write_output_manifest <- function(root) {
  f <- list.files(root,recursive=TRUE,full.names=TRUE)
  f <- f[file.info(f)$isdir %in% FALSE]
  readr::write_csv(tibble::tibble(relative_path=sub(paste0("^",normalizePath(root),"/?"),"",normalizePath(f)),
    size_bytes=file.info(f)$size,sha256=vapply(f,ish_sha256,character(1))),
    file.path(root,"OUTPUT_MANIFEST.csv"))
}
