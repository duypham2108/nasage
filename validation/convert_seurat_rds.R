# Export raw counts and cell metadata from a Seurat v3/v4 .rds without needing Seurat itself
# (SeuratObject + Matrix are enough). Python rebuilds an AnnData with
# validation/rds_export_to_h5ad.py.
#
# Usage: Rscript validation/convert_seurat_rds.R <in.rds> <out_dir> [assay]
# Writes to <out_dir>: counts_{i,p,x}.bin (CSC of genes x cells), counts_dim.txt,
#                      genes.txt, cells.txt, obs.csv, summary.txt
suppressPackageStartupMessages({ library(Matrix); library(SeuratObject) })
args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; out <- args[2]; assay <- if (length(args) > 2) args[3] else NULL
dir.create(out, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

msg("reading", inp)
obj <- readRDS(inp)
msg("class", class(obj), "assays:", paste(names(obj@assays), collapse = ","),
    "default:", obj@active.assay)
if (is.null(assay)) assay <- if ("RNA" %in% names(obj@assays)) "RNA" else obj@active.assay
a <- obj@assays[[assay]]
m <- a@counts
if (length(m) == 0 || ncol(m) == 0) stop("assay ", assay, " has no counts slot")
m <- as(m, "CsparseMatrix")
msg("counts", nrow(m), "genes x", ncol(m), "cells; nnz", length(m@x),
    "; integer-valued:", all(m@x[seq_len(min(1e6, length(m@x)))] %% 1 == 0))

writeBin(as.integer(m@i), file.path(out, "counts_i.bin"), size = 4)
writeBin(as.integer(m@p), file.path(out, "counts_p.bin"), size = 4)
writeBin(as.numeric(m@x), file.path(out, "counts_x.bin"), size = 4)   # float32
writeLines(as.character(dim(m)), file.path(out, "counts_dim.txt"))
writeLines(rownames(m), file.path(out, "genes.txt"))
writeLines(colnames(m), file.path(out, "cells.txt"))
md <- obj@meta.data
stopifnot(identical(rownames(md), colnames(m)))
write.csv(md, file.path(out, "obs.csv"), row.names = TRUE)

sink(file.path(out, "summary.txt"))
str(md, list.len = 200)
for (cn in colnames(md)) if (!is.numeric(md[[cn]]) || length(unique(md[[cn]])) < 60) {
  cat("\n##", cn, "\n"); print(table(md[[cn]], useNA = "ifany"))
}
sink()
msg("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
