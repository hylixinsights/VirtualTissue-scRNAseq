#!/usr/bin/env Rscript
# Structural inspection only. No state composition, change, embeddings or outcomes.
.libPaths(c(normalizePath(".Rlib"), .Library))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript r/inventory_seurat.R INPUT.rds OUTPUT.json")
obj <- readRDS(args[[1]])
if (!inherits(obj, "Seurat")) stop("Expected Seurat; got ", paste(class(obj), collapse = ", "))
meta <- obj@meta.data
assays <- lapply(names(obj@assays), function(name) {
  assay <- obj@assays[[name]]
  list(name = name, class = class(assay), dimensions = dim(assay),
       layers = if (inherits(assay, "Assay5")) Layers(assay) else intersect(c("counts", "data", "scale.data"), slotNames(assay)))
})
columns <- lapply(names(meta), function(name) {
  x <- meta[[name]]
  list(name = name, class = class(x), missing = sum(is.na(x)), unique_values = length(unique(x)))
})
inventory <- list(source = args[[1]], object_class = class(obj),
                  stored_seurat_version = as.character(obj@version), cells = nrow(meta),
                  assays = assays, reductions = names(obj@reductions), metadata_columns = columns,
                  time_field_candidates = grep("date|day|interval|time", names(meta), ignore.case = TRUE, value = TRUE),
                  sample_field_candidates = grep("sample|patient|orig.ident", names(meta), ignore.case = TRUE, value = TRUE),
                  inspected_state_compositions = FALSE)
write_json(inventory, args[[2]], auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("Saved structural inventory; no state compositions computed.\n")
