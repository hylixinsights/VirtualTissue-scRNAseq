# Technical inventory only. Never print cell-by-state tables or model outcomes.
.libPaths(c('.Rlib-expansion', '.Rlib', .libPaths()))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1L)
input <- normalizePath(args[[1]], mustWork=TRUE)
record <- jsonlite::read_json('reports/evidence/cohort-expansion/lin-object-integrity.json')
stopifnot(file.info(input)$size == record$bytes)
stopifnot(unname(tools::md5sum(input)) == record$md5)
cat('Verified source. Reading qs2 with internal checksum verification.\n')
flush.console()
object <- qs2::qs_read(input, validate_checksum=TRUE, nthreads=2L)
stopifnot(inherits(object, 'Seurat'))
metadata <- object@meta.data
stopifnot(!anyDuplicated(rownames(metadata)))
stopifnot(identical(colnames(object), rownames(metadata)))
clinical <- read.delim('tables/cohort_expansion_lin_samples.tsv', stringsAsFactors=FALSE)
matches <- list()
for (column in names(metadata)) {
  values <- metadata[[column]]
  if (is.character(values) || is.factor(values)) {
    vocabulary <- unique(as.character(values))
    exact <- intersect(vocabulary, clinical$sample_id)
    if (length(exact)) matches[[column]] <- exact
  }
}
assays <- lapply(object@assays, function(a) {
  list(class=class(a), dimensions=dim(a), slots=slotNames(a),
       layers=if (inherits(a, 'Assay5')) names(a@layers) else NULL)
})
schema <- list(source_md5=record$md5, class=class(object),
               metadata_rows=nrow(metadata), metadata_columns=names(metadata),
               metadata_classes=lapply(metadata, class),
               exact_clinical_sample_id_matches=matches, assays=assays,
               reductions=names(object@reductions),
               object_size_bytes=as.numeric(object.size(object)),
               state_counts_computed=FALSE)
jsonlite::write_json(schema, 'reports/evidence/cohort-expansion/lin-object-schema.json',
                     pretty=TRUE, auto_unbox=TRUE, null='null')
target <- '../data/Lin2026/derived/audit-quarantine'
dir.create(target, recursive=TRUE, showWarnings=FALSE)
# Preserve source metadata for repeated whitelist-only audits. This is not a
# training-ready export and may contain authors' state labels and response.
saveRDS(metadata, file.path(target, 'source_metadata.rds'))
capture.output(sessionInfo(), file='reports/evidence/cohort-expansion/lin-inventory-sessionInfo.txt')
cat('Metadata alignment verified; schema and quarantined metadata saved.\n')
