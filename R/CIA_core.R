#' Load Gene Signatures
#'
#' This function loads gene signatures from a given input, which can either be a file path to a TSV file or a named list.
#' Each signature in the file should be in a separate row with the first column being the signature name and
#' the remaining columns containing the genes.
#'
#' @param signatures_input A character string representing the file path to the TSV file containing the gene
#' signatures, or a list where each element is a vector of gene names with names of the list elements
#' representing signature names.
#' @return A list where each element is a character vector of gene names, named by the signature names.
#' @importFrom data.table fread
#' @export

load_signatures <- function(signatures_input) {
  if (is.character(signatures_input)) {
    df <- fread(signatures_input, sep = '\t', header = FALSE)
    signatures <- split(df[, -1], df[[1]])
    signatures <- lapply(signatures, na.omit)
    signatures <- lapply(signatures, as.character)
  } else if (is.list(signatures_input)) {
    signatures <- signatures_input
  } else {
    stop("signatures_input must be either a list or a string path to a TSV file.")
  }
  return(signatures)
}

#' Compute Signature Scores
#'
#' Calculates signature scores for a given set of genes by comparing their expression levels
#' in the provided data. The score for each gene set is calculated as the product of the
#' count of genes expressed (presence) in a given cell and their relative expression levels (expression),
#' normalized by the total expression of the cell.
#'
#' @param data Numeric matrix, data frame, SeuratObject, or SingleCellExperiment with genes in rows and cells in columns.
#' @param geneset Vector of gene names to compute the scores for.
#' @param seurat_assay Assay to use for SeuratObject objects (default "RNA").
#' @param matrix Slot to use for Seurat or SingleCellExperiment objects (default "data" for SO, 'logcounts' for SCE).
#' @param total_col_sums Optional precomputed column sums for normalization.
#' @return Vector of signature scores for each column (cell) in the data.
#' @importFrom sparseMatrixStats colSums2
#' @export

compute_signature_scores <- function(data, geneset, seurat_assay="RNA", matrix="data", total_col_sums = NULL) {
  if (inherits(data, "Seurat")) {
    datam <- slot(data[[seurat_assay]],matrix)
  } else if (inherits(data, "SingleCellExperiment")) {
    if(matrix=='data'){matrix<-'logcounts'}
    datam <- assay(data, matrix)
  } else {datam <- data}
  geneset <- intersect(geneset, rownames(datam))
  if (length(geneset) == 0) {
    return(rep(0, ncol(datam)))
  }

  if (is.null(total_col_sums)) {
    total_col_sums <- colSums2(datam)
  }

  subdata <- datam[geneset, , drop = FALSE]
  count <- colSums2(subdata > 0)
  exp <- colSums2(subdata) / total_col_sums

  return(count * exp)
}


## ' Calculate Signature Scores
#'
#' This function calculates the signature scores for each gene set in the provided data.
#' It loads the gene signatures from an input source and computes the scores based on
#' the expression levels of these genes in the data.
#'
#' @param data A numeric matrix, data frame, SeuratObject, or SingleCellExperiment with genes in rows and cells in columns,
#' representing the expression level of each gene in each cell.
#' @param signatures_input A character string representing the file path to the TSV file containing the gene
#' signatures, or a list where each element is a vector of gene names with names of the list elements
#' representing signature names.
#' @param return_score Boolean to return scores directly (default FALSE).
#' @param seurat_assay Assay for Seurat objects (default "RNA").
#' @param matrix Slot to use for Seurat or SingleCellExperiment objects (default "data" for SO, 'logcounts' for SCE).
#' @param score_mode Calculation mode for scores: 'raw', 'scaled', or log-transformed ('log', 'log2', 'log10').
#' @param n_cpus Number of CPU cores for parallel computation (default uses a quarter of available cores).
#' @return Matrix of signature scores if return_score=TRUE or if data is a matrix/data,frame.Otherwise, updates the input object's metadata with scores.
#' @importFrom future future_lapply plan
#' @importFrom data.table fread
#' @importFrom sparseMatrixStats colSums2
#' @export

signature_score <- function(data, signatures_input, return_score=FALSE, seurat_assay="RNA",matrix="data", score_mode = 'raw', n_cpus = NULL) {
  signatures <- load_signatures(signatures_input)
    # Check the type of data and extract expression matrix accordingly
  if (inherits(data, "Seurat")) {
    datam <- slot(data[[seurat_assay]],matrix)
  } else if (inherits(data, "SingleCellExperiment")) {
    if(matrix=='data'){matrix<-'logcounts'}
    datam <- assay(data, matrix)
  } else {datam <- data}

  cat('Checking if genes are in the dataset matrix...', "\n")
  result <- sapply(names(signatures), function(x) {
    lt <- length(signatures[[x]])
    l <- sum(signatures[[x]] %in% rownames(data))
    message <- sprintf("%s: %d / %d", x, l, lt)
    cat(message, "\n")
    invisible(message)  # Use invisible to avoid printing the return value of cat
})

  n_cpus <- if (is.null(n_cpus)) { ceiling(availableCores() / 4) } else { n_cpus }
  tot <- colSums2(datam)
  plan(multicore, workers = n_cpus)

  scores <- future_lapply(names(signatures), future.seed=TRUE, FUN =function(name) {
    geneset <- signatures[[name]]
    compute_signature_scores(datam, geneset, tot)
  })

  scores_df <- do.call(cbind, scores)
  colnames(scores_df) <- names(signatures)

  if (score_mode == 'scaled') {
    scores_df <- sweep(scores_df, 2, apply(scores_df, 2, max), `/`)
  } else if (score_mode %in% c('log', 'log2', 'log10')) {
    scores_df <- match.fun(score_mode)(scores_df + .Machine$double.xmin)
  } else if (score_mode != 'raw') {
    stop("Invalid score_mode. Must be one of: ['raw', 'scaled', 'log', 'log2', 'log10']")
  }


  if (return_score){
    return(scores_df)
  }

if (inherits(data, "Seurat")) {
    data@meta.data[, colnames(scores_df)] <- scores_df
    cat('Scores have been added in data@meta.data', "\n")
    return(data)

  } else if (inherits(data, "SingleCellExperiment")) {
    colData(data)[, colnames(scores_df)] <- scores_df
    cat('Scores have been added in colData(data)', "\n")
    return(data)
  } else {return(scores_df)}

}


#' Classification Based on Signature Scores
#'
#' Performs classification of cells in `data` based on the computed signature scores. Each sample/cell
#' is classified to the signature with the highest score, provided the difference between the top two scores
#' exceeds a given threshold, otherwise, it is marked as unassigned.
#'
#' @param data A numeric matrix, data frame, SeuratObject, or SingleCellExperiment with genes in rows and cells in columns,
#' representing the expression level of each gene in each cell.
#' @param signatures_input A character string representing the file path to the TSV file containing the gene
#' signatures, or a list where each element is a vector of gene names with names of the list elements
#' representing signature names.
#' @param seurat_assay Assay for Seurat objects (default "RNA").
#' @param matrix Slot to use for Seurat or SingleCellExperiment objects (default "data" for SO, 'logcounts' for SCE).
#' @param n_cpus An optional integer indicating the number of CPU cores to use for parallel computation. If NULL, the function will decide the number based on
#' available system resources.
#' @param similarity_threshold A numeric threshold used to decide if the highest score is significantly higher than the second-highest to
#' assign a specific class to a cell. If the difference is less than or equal to this threshold, the cell
#' is labeled as unassigned.
#' @param column_name The name of the column to be added to the metadata storing the classification labels, default is 'CIA_prediction'.
#' @param unassigned_label The label to assign to the cells where no clear majority signature is identified, default is 'Unassigned'.
#' @return Modifies the input data object by adding the classification labels to its metadata and returns the modified data object.
#' If input data is a matrix or data frame, returns a vector of classification labels.
#' @export

signature_based_classification <- function(data, signatures_input, n_cpus = NULL, similarity_threshold = 0.1,
                                           seurat_assay="RNA",matrix="data",column_name='CIA_prediction',
                                           unassigned_label='Unassigned') {
  start_time <- Sys.time()  # Capture start time

  get_label <- function(row) {
    order_scores <- order(row, decreasing = TRUE)
    if (row[order_scores[1]] - row[order_scores[2]] <= similarity_threshold) {
      return(unassigned_label)
    } else {
      return(names(row)[order_scores[1]])
    }
  }

  score_matrix <- signature_score(data, signatures_input, score_mode = 'scaled',seurat_assay=seurat_assay,
                                  matrix=matrix, return_score=T, n_cpus = n_cpus)
  labels <- apply(score_matrix, 1, get_label)

  end_time <- Sys.time()  # Capture end time

  # Format and print the message correctly
  cat("\nClassification complete!    Start:", format(start_time, "%H:%M:%S"),
      "    End:", format(end_time, "%H:%M:%S"), "\n")

if (inherits(data, "Seurat")) {
    data@meta.data[,column_name] <-as.factor(labels)
    cat(column_name,'has been added in data@meta.data', "\n")
    return(data)
  } else if (inherits(data, "SingleCellExperiment")) {
    colData(data)[,column_name]<- as.data.frame(labels)
    cat(column_name,'has been added in colData(data)', "\n")
    return(data)
  }
  else {return(as.factor(labels))}
}

