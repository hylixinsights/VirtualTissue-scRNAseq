#!/usr/bin/env Rscript
# Independent count-only check against the original RDS, not the exported matrix.
.libPaths(c(normalizePath(".Rlib"), .Library))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(jsonlite))
obj <- readRDS("../data/SPECIAL/raw/GC_all_immune_sharing.rds")
meta <- obj@meta.data
is_cd8 <- grepl("CD8_Tcells$", as.character(meta$high_resolution_clusters))
totals <- aggregate(as.integer(is_cd8), by = list(sample_id = meta$sample_ID), FUN = sum)
names(totals)[[2]] <- "r_cd8_total"
python <- read.delim("tables/samples.tsv", check.names = FALSE)
joined <- merge(totals, python[, c("sample_id", "author_annotated_cd8_cells")], by = "sample_id", all.x = TRUE)
stopifnot(nrow(joined) == length(unique(meta$sample_ID)), !anyNA(joined),
          all(joined$r_cd8_total == joined$author_annotated_cd8_cells))
write_json(list(status = "passed", source = "original RDS metadata", method = "suffix-based CD8 lineage total",
                biopsies_compared = nrow(joined), all_totals_identical = TRUE, per_biopsy = joined,
                state_compositions_computed = FALSE),
           "reports/evidence/cohort-counts-r-validation.json", auto_unbox = TRUE, pretty = TRUE)
cat("Independent RDS count check passed.\n")
