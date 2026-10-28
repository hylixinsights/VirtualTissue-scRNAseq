#!/usr/bin/env Rscript
# A05 technical export. Explicitly excludes all clinical outcomes, CD8 states and embeddings.
.libPaths(c(normalizePath('.Rlib'), .Library))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4L) stop('Usage: INPUT.rds OUTPUT_DIRECTORY SAMPLE_COLUMN REPORT.json')
input<-normalizePath(args[1]);out<-args[2];sample_col<-args[3];report_path<-args[4]
if(dir.exists(out) && length(list.files(out))) stop('Refuse nonempty output directory')
if(file.exists(report_path)) stop('Refuse existing accepted report')
mf<-fromJSON('data/MANIFEST.json',simplifyVector=FALSE)
records<-Filter(function(x)identical(x$filename,basename(input)),mf$files)
stopifnot(length(records)==1L,identical(records[[1]]$local$status,'verified'))
sha<-digest::digest(file=input,algo='sha256',serialize=FALSE)
stopifnot(identical(sha,records[[1]]$local$sha256))
obj<-readRDS(input);stopifnot(inherits(obj,'Seurat'))
a<-obj@assays[['RNA']];counts<-if(inherits(a,'Assay5'))LayerData(a,layer='counts') else a@counts
meta<-obj@meta.data
stopifnot(sample_col %in% names(meta),identical(colnames(counts),rownames(meta)),!anyDuplicated(colnames(counts)),!anyDuplicated(rownames(counts)),all(is.finite(counts@x)),all(counts@x>=0),all(counts@x==floor(counts@x)))
# Drop unneeded assay data and quarantined metadata immediately; never print their values.
sample_ids<-as.character(meta[[sample_col]]);stopifnot(!anyNA(sample_ids),all(nzchar(sample_ids)))
rm(obj,a,meta);gc()
dir.create(out,recursive=TRUE,showWarnings=FALSE)
samples<-sort(unique(sample_ids));checks<-list();exported<-character()
for(k in seq_along(samples)){
 sample<-samples[k];cells<-colnames(counts)[sample_ids==sample];m<-counts[,cells,drop=FALSE]
 # Safe opaque output directories; sample labels retained only in technical metadata.
 dest<-file.path(out,sprintf('sample_%03d',k));dir.create(dest)
 file<-file.path(dest,'counts.mtx');writeMM(m,file)
 obs<-data.frame(cell_id=cells,sample_id=sample,qc_total_counts=as.numeric(colSums(m)),qc_detected_genes=as.integer(colSums(m>0)))
 mt<-grepl('^MT-',rownames(m));obs$qc_mito_fraction<-as.numeric(colSums(m[mt,,drop=FALSE]))/pmax(1,obs$qc_total_counts)
 # Fixed lineage screen; not state assignment or outcome analysis.
 tmarkers<-intersect(c('CD3D','CD3E','TRAC'),rownames(m));cd8markers<-intersect(c('CD8A','CD8B'),rownames(m))
 obs$provisional_cd8_screen<-as.integer(colSums(m[tmarkers,,drop=FALSE]>0)>=2 & colSums(m[cd8markers,,drop=FALSE]>0)>=1)
 write_parquet(obs,file.path(dest,'obs.parquet'))
 write.table(data.frame(gene_id=rownames(m)),file.path(dest,'var.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 restored<-as(readMM(file),'CsparseMatrix');readobs<-as.data.frame(read_parquet(file.path(dest,'obs.parquet')))
 stopifnot(identical(dim(m),dim(restored)),Matrix::nnzero(m-restored)==0,identical(readobs$cell_id,cells),identical(read.table(file.path(dest,'var.tsv'),header=TRUE,sep='\t',stringsAsFactors=FALSE)$gene_id,rownames(m)))
 checks[[length(checks)+1L]]<-list(directory=basename(dest),cells=ncol(m),genes=nrow(m),exact_roundtrip=TRUE,files=lapply(list.files(dest,full.names=TRUE),function(p)list(name=basename(p),bytes=unname(file.info(p)$size),sha256=digest::digest(file=p,algo='sha256',serialize=FALSE))))
 exported<-c(exported,cells)
}
stopifnot(length(exported)==ncol(counts),!anyDuplicated(exported))
write_json(list(input_sha256=sha,cells=ncol(counts),genes=nrow(counts),samples=length(samples),matrix_orientation='genes_by_cells',exact_roundtrip=TRUE,clinical_metadata_exported=FALSE,states_inspected=FALSE,paired_changes_computed=FALSE,licence='CC-BY-NC-SA-4.0',per_sample=checks),report_path,pretty=TRUE,auto_unbox=TRUE)
writeLines(c('# A05 outcome-blind technical export','Source: SPECIAL doi:10.48804/GSAXBN v8. CC BY-NC-SA 4.0; no commercial use.','Only raw RNA, sample/barcode IDs and technical QC exported. No outcomes, state labels or integrated embeddings.','Counts: genes x cells, aligned exactly with var.tsv rows and obs.parquet rows.','Provisional lineage screen is fixed and does not establish clinical eligibility.','The original source retains metadata requiring a future preregistered access decision.'),file.path(out,'README.md'))
cat('Technical export complete; exact sparse-matrix round trips verified.\n')
