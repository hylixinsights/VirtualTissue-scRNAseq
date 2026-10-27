#!/usr/bin/env Rscript
# Run from the repository root. Install into a project-specific library only.
local_lib <- normalizePath(".Rlib", mustWork = FALSE)
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .Library))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
packages <- c("renv", "Seurat", "SeuratObject", "Matrix", "arrow", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = local_lib, Ncpus = 2)
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
renv::snapshot(packages = packages, lockfile = "renv.lock", prompt = FALSE)
sink("sessionInfo.txt")
cat("Library paths:\n"); print(.libPaths())
for (p in packages) cat(p, as.character(packageVersion(p)), "\n")
print(sessionInfo())
sink()
