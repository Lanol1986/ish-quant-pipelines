
ish_source_core <- function(repo_root=getwd()) {
  files <- c("core.R","io.R","analysis.R","qc.R","plots.R","provenance.R",
             "single_panel.R","integration.R")
  for (f in files) source(file.path(repo_root,"R",f),local=.GlobalEnv)
  invisible(TRUE)
}
