#!/usr/bin/env Rscript
# T1.1: technical export only; no state composition or pre/on comparison.
.libPaths(c(normalizePath(".Rlib"), .Library))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript r/export_seurat.R INPUT.rds OUTPUT_DIR")
input <- normalizePath(args[[1]])
out <- args[[2]]
if (dir.exists(out) && length(list.files(out))) stop("Output directory must be empty")
manifest <- fromJSON("data/MANIFEST.json", simplifyVector = FALSE)
record <- Filter(function(x) identical(x$filename, basename(input)), manifest$files)[[1]]
stopifnot(identical(record$local$status, "verified"))
input_sha <- digest::digest(file = input, algo = "sha256", serialize = FALSE)
stopifnot(identical(input_sha, record$local$sha256))
obj <- readRDS(input)
stopifnot(inherits(obj, "Seurat"))
meta <- obj@meta.data
required <- c("sample_ID", "patient_ID", "Timepoint")
stopifnot(all(required %in% names(meta)), !anyNA(meta[, required]))
assay <- obj@assays[["RNA"]]
counts <- if (inherits(assay, "Assay5")) LayerData(assay, layer = "counts") else assay@counts
stopifnot(inherits(counts, "sparseMatrix"), !anyDuplicated(colnames(counts)),
          !anyDuplicated(rownames(counts)), identical(colnames(counts), rownames(meta)),
          all(is.finite(counts@x)), all(counts@x >= 0), all(counts@x == floor(counts@x)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
samples <- sort(unique(as.character(meta$sample_ID)))
stopifnot(all(grepl("^[A-Za-z0-9_-]+$", samples)))
checks <- list(); exported <- character()
for (sample in samples) {
  cells <- rownames(meta)[meta$sample_ID == sample]
  m <- counts[, cells, drop = FALSE]
  dest <- file.path(out, sample); dir.create(dest)
  matrix_file <- file.path(dest, "counts.mtx")
  writeMM(m, matrix_file)
  obs <- data.frame(cell_id = cells, meta[cells, , drop = FALSE], check.names = FALSE, row.names = NULL)
  obs[] <- lapply(obs, function(x) if (is.factor(x)) as.character(x) else x)
  obs$qc_total_counts <- as.numeric(colSums(m))
  obs$qc_detected_genes <- as.integer(colSums(m > 0))
  mt <- grepl("^MT-", rownames(m))
  obs$qc_mito_fraction <- ifelse(obs$qc_total_counts > 0, colSums(m[mt, , drop = FALSE]) / obs$qc_total_counts, NA_real_)
  write_parquet(obs, file.path(dest, "obs.parquet"))
  write.table(data.frame(gene_id = rownames(m)), file.path(dest, "var.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  for (reduction in names(obj@reductions)) {
    embedding <- obj@reductions[[reduction]]@cell.embeddings
    if (!all(cells %in% rownames(embedding))) stop("Incomplete embedding: ", reduction)
    write_parquet(data.frame(cell_id = cells, embedding[cells, , drop = FALSE], row.names = NULL),
                  file.path(dest, paste0("author_", reduction, ".parquet")))
  }
  restored <- readMM(matrix_file)
  obs_restored <- as.data.frame(read_parquet(file.path(dest, "obs.parquet")))
  roundtrip_ok <- identical(dim(m), dim(restored)) &&
    isTRUE(all.equal(as.numeric(colSums(m)), as.numeric(colSums(restored)), tolerance = 0)) &&
    isTRUE(all.equal(as.numeric(rowSums(m)), as.numeric(rowSums(restored)), tolerance = 0)) &&
    Matrix::nnzero(m - restored) == 0 && identical(obs_restored$cell_id, cells) &&
    isTRUE(all.equal(obs, obs_restored, check.attributes = FALSE, tolerance = 0))
  stopifnot(roundtrip_ok)
  files <- list.files(dest, full.names = TRUE)
  checks[[sample]] <- list(sample_id = sample, cells = length(cells), genes = nrow(m),
                           count_total = sum(m), exact_roundtrip = roundtrip_ok,
                           files = lapply(files, function(p) list(path = basename(p), bytes = unname(file.info(p)$size),
                             sha256 = digest::digest(file = p, algo = "sha256", serialize = FALSE))))
  exported <- c(exported, cells)
}
stopifnot(length(exported) == ncol(counts), !anyDuplicated(exported), setequal(exported, colnames(counts)))
report <- list(input_sha256 = input_sha, total_cells = ncol(counts), exported_cells = length(exported),
               genes = nrow(counts), samples = length(samples), all_cells_exported = TRUE,
               exact_count_roundtrip = TRUE, matrix_orientation = "genes_by_cells",
               per_sample = checks, state_compositions_computed = FALSE,
               embeddings_policy = "Archival author embeddings only; refit every learned representation on training patients")
write_json(report, "reports/evidence/immune-export-validation.json", auto_unbox = TRUE, pretty = TRUE)
writeLines(c("# SPECIAL neutral-format export", "",
             "MatrixMarket: genes by cells. var.tsv order matches matrix rows; obs.parquet order matches columns.",
             "Source: doi:10.48804/GSAXBN, version 8; CC BY-NC-SA 4.0.",
             "Changes: per-sample RNA raw-count export and derived technical QC fields; no filtering.",
             "Author embeddings are archival and must not be used as a fitted training/test representation.",
             "State/outcome fields are quarantined until eligibility and lockbox assignment."), file.path(out, "README.md"))
cat("Export complete; exact count and metadata round trips passed.\n")
