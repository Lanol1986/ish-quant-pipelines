#!/usr/bin/env Rscript
library(testthat)
args <- commandArgs(trailingOnly=FALSE)
fa <- grep("^--file=",args,value=TRUE)
script <- if (length(fa)) normalizePath(sub("^--file=","",fa[1])) else normalizePath("tests/run_tests.R")
repo_root <- dirname(dirname(script))
source(file.path(repo_root,"R","load.R"))
ish_source_core(repo_root)
test_dir(file.path(repo_root,"tests","testthat"),reporter="summary")
