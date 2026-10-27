# Restore only the reader extension; the original renv.lock remains unchanged.
args <- commandArgs(trailingOnly=TRUE)
target <- if(length(args)) args[[1]] else '.Rlib-expansion'
dir.create(target, recursive=TRUE, showWarnings=FALSE)
.libPaths(c('.Rlib', .Library))
options(repos=c(CRAN='https://cloud.r-project.org'), timeout=600)
renv::restore(lockfile='renv-expansion.lock', library=normalizePath(target), prompt=FALSE)
.libPaths(c(normalizePath(target), '.Rlib', .Library))
stopifnot(as.character(packageVersion('qs2')) == '0.3.1')
x <- list(patient_id=c('synthetic_a', 'synthetic_b'), count=c(30L, 45L))
f <- tempfile(fileext='.qs2')
qs2::qs_save(x, f)
stopifnot(identical(x, qs2::qs_read(f, validate_checksum=TRUE)))
cat('Pinned reader extension restored; synthetic checksum round-trip passed.\n')
