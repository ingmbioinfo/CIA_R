## Setup for running the tests - creating object to be used suite-wide

suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
  library("sparseMatrixStats")
  library("SingleCellExperiment")
  library("scRNAseq")
})

sce <- ReprocessedAllenData()
sce

# Some simple reprocessing:

library("scater")
sce <- logNormCounts(sce, exprs_values="tophat_counts")
sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)
sce
rownames(sce) <- toupper(rownames(sce))

gmt <- load_signatures(system.file("extdata", "azimuth_human_motor_cortex.gmt",
                                   package = "CIA"))
gmt

mat_logcounts <- logcounts(sce)
