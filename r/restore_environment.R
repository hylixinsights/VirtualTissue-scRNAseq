#!/usr/bin/env Rscript
# Run in the repository root, optionally passing a fresh library directory.
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1]] else ".Rlib"
dir.create(target, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlib"), .Library))
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("Bootstrap renv into .Rlib first; see README.md")
}
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
renv::restore(lockfile = "renv.lock", library = normalizePath(target), prompt = FALSE)
.libPaths(c(normalizePath(target), .Library))
pkgs <- c("Seurat", "SeuratObject", "Matrix", "arrow", "jsonlite")
for (p in pkgs) stopifnot(requireNamespace(p, quietly = TRUE))
cat("Restored R library from renv.lock\n")
