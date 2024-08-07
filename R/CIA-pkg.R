#' CIA
#'
#' Cluster Independent Annotation
#'
#' `CIA` (Cluster Independent Annotation) is a cutting-edge computational tool
#' designed to accurately classify cells in scRNA-seq datasets using gene
#' signatures.
#' This tool operates without the need for a fully annotated reference dataset
#' or complex machine learning processes, providing a highly user-friendly and
#' practical solution for cell type annotation.
#'
#' `CIA` summarizes the information of each signature expression into a single
#' score value for each cell.
#' By comparing these score values, `CIA` assigns labels to each cell based
#' on the top-scored signature.
#' `CIA` can filter scores by their distribution or significance, allowing
#' comparison of genesets with lengths spanning tens to thousands of genes.
#'
#' `CIA` is implemented in both R and Python, making it compatible with all
#' major single-cell analysis frameworks like `SingleCellExperiment`, `Seurat`,
#' and `Scanpy`.
#' This compatibility ensures a seamless integration into existing workflows.
#'
#' @name CIA-pkg
#' @docType package
#' @keywords internal
"_PACKAGE"

globalVariables(c("Cluster", "Count"))
