#' Load Gene Signatures
#'
#' This function loads gene signatures from a given input, which can either be a
#' file path to a TSV file or a named list.
#' Each signature in the file should be in a separate row with the first column
#' being the signature name and the remaining columns containing the genes.
#'
#' @param signatures_input A character string representing the file path to the
#' TSV file containing the gene signatures.
#' @param description_field_available Logical value, to accommodate for the use
#' of "custom" gmt file formats, where the description field might not be
#' provided. Defaults to TRUE.
#'
#' @return A list where each element is a character vector of gene names, named
#' by the signature names.
#'
#' @export
#'
#' @importFrom utils head tail
#' @importFrom methods is
#'
#' @examples
#' gmt_file <- system.file("extdata", "azimuth_human_motor_cortex.gmt",
#'                         package = "CIA")
#' signatures <- load_signatures(gmt_file)
#' signatures
#'
#' gmt_url <-
#'   "https://data.wikipathways.org/20240710/gmt/wikipathways-20240710-gmt-Homo_sapiens.gmt"
#'
#' sig_wikipathways <- load_signatures(gmt_url)
#' head(names(sig_wikipathways))
load_signatures <- function(signatures_input, description_field_available = TRUE) {
  if (!is.character(signatures_input)) {
    stop("signatures_input must be either a URL or a string path to a local TSV file.")
  }

  if (!grepl("^http", signatures_input)) {
    stopifnot(file.exists(signatures_input))
  }

  input_lines <- strsplit(readLines(signatures_input), "\t")

  if (description_field_available) {
    signatures <- lapply(input_lines, tail, -2)
  } else {
    # strip only the first field
    signatures <- lapply(input_lines, tail, -1)
    message(
      "You are loading a gmt file which is not entirely conform to the ",
      "official file specification. Please consider reformatting that ",
      "if possible."
    )
  }

  names(signatures) <- lapply(input_lines, head, 1)

  # if some fields were left empty for a conversion issue...
  for (i in names(signatures)) {
    signatures[[i]] <- signatures[[i]][signatures[[i]] != ""]
  }

  return(signatures)
}


#' Compute Signature Scores
#'
#' Calculates signature scores for a given set of genes by comparing their
#' expression levels in the provided data. The score for each gene set is
#' calculated as the product of the count of genes expressed (presence) in a
#' given cell and their relative expression levels (expression), normalized by
#' the total expression of the cell.
#'
#' @param data Numeric matrix, data frame, SeuratObject, or SingleCellExperiment
#' with genes in rows and cells in columns.
#' @param geneset Vector of gene names to compute the scores for.
#' @param seurat_assay Assay to use for SeuratObject objects (default "RNA").
#' @param seurat_layer Character string, indicating which layer to use of the
#' Seurat object. Defaults to `data`.
#' @param sce_assay Character string, indicating which assay to use of the
#' SingleCellExperiment object. Defaults to `logcounts`.
#' @param total_col_sums Optional precomputed column sums for normalization.
#'
#' @return Vector of signature scores for each column (cell) in the data.
#'
#' @importFrom sparseMatrixStats colSums2
#' @importFrom methods slot is
#' @importFrom SeuratObject Layers LayerData Version
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Seurat GetAssay
#'
#' @export
#'
#' @examples
#' ## TODO example
#'
score_signature <- function(data,
                            geneset,
                            seurat_assay = "RNA",
                            seurat_layer = "data",
                            sce_assay = "logcounts",
                            total_col_sums = NULL) {
  ## Checks on arguments
  allowed_formats <- c("SingleCellExperiment", "Seurat", "matrix", "Matrix", "DelayedMatrix")
  if (!any(unlist((lapply(allowed_formats, function(arg) is(data, arg)))))) {
    stop(
      "The data provided should be in one of the following formats: ",
      paste(allowed_formats, collapse = "|")
    )
  }

  stopifnot(is.character(geneset))
  stopifnot(length(geneset) > 0)
  stopifnot(is.character(seurat_assay))
  stopifnot(length(seurat_assay) == 1)
  stopifnot(is.character(seurat_layer))
  stopifnot(length(seurat_layer) == 1)
  stopifnot(is.character(sce_assay))
  stopifnot(length(sce_assay) == 1)

  if (!is.null(total_col_sums)) {
    stopifnot(length(total_col_sums) == ncol(data))
  }

  if (is(data, "Seurat")) {
    stopifnot(is(data[[seurat_assay]], "Assay") | is(data[[seurat_assay]], "Assay5"))
    stopifnot(seurat_layer %in% Layers(data))

    obj_version_major <- as.numeric(
      strsplit(as.character(Version(data)), ".", fixed = TRUE)[[1]][1]
    )

    if (obj_version_major < 5) {
      message("Seurat object version is < 5.0.0 . It's suggested to run Seurat::UpdateSeuratObject first.")
      datam <- slot(data[[seurat_assay]], seurat_layer)
    } else {
      datam <- LayerData(data[[seurat_assay]], seurat_layer)
    }
  } else if (is(data, "SingleCellExperiment")) {
    stopifnot(sce_assay %in% assayNames(data))

    datam <- assay(data, sce_assay)
  } else if (is(data, "matrix") | is(data, "Matrix") | is(data, "DelayedMatrix")) {
    stopifnot(!is.null(dim(data)) & all(dim(data) > 0))
    if (is(data, "matrix")) {
      stopifnot(is.numeric(data))
    }
    if (is(data, "Matrix")) {
      stopifnot(is.numeric(data[1, 1]))
    }

    datam <- data
  }

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

  sigscore <- count * exp

  return(sigscore)
}


#' Calculate Signature Scores
#'
#' This function calculates the signature scores for each gene set in the
#' provided data. It loads the gene signatures from an input source and computes
#' the scores based on the expression levels of these genes in the data.
#'
#' @param data A numeric matrix, data frame, `SeuratObject`, or
#' `SingleCellExperiment` with genes in rows and cells in columns,
#' representing the expression level of each gene in each cell.
#' @param signatures_input A character string representing the file path to the
#' TSV file containing the gene signatures, or a list where each element is a
#' vector of gene names with names of the list elements representing signature
#' names.
#' @param return_score Boolean to return scores directly (default FALSE).
#' @param seurat_assay Assay for Seurat objects (default "RNA").
#' @param seurat_layer Character string, indicating which layer to use of the
#' Seurat object. Defaults to `data`.
#' @param sce_assay Character string, indicating which assay to use of the
#' SingleCellExperiment object. Defaults to `logcounts`.
#' @param score_mode Character string, specifying the calculation mode to be
#' used for scores: "raw", "scaled", or log-transformed ("log", "log2", "log10").
#' In the scope of this function, this defaults to "raw".
#'
#' @param n_cpus Number of CPU cores for parallel computation (default uses a
#' quarter of available cores).
#'
#' @return Matrix of signature scores if return_score=TRUE or if data is a
#' matrix/data,frame. Otherwise, updates the input object's metadata with scores.
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom parallel detectCores
#' @importFrom sparseMatrixStats colSums2
#' @importFrom SummarizedExperiment assay colData<- colData assayNames
#' @importFrom SingleCellExperiment SingleCellExperiment counts logcounts
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' ## TODO example
#'
score_all_signatures <- function(data,
                                 signatures_input,
                                 return_score = FALSE,
                                 seurat_assay = "RNA",
                                 seurat_layer = "data",
                                 sce_assay = "logcounts",
                                 score_mode = "raw",
                                 n_cpus = NULL) {
  ## Checks on arguments
  allowed_formats <- c("SingleCellExperiment", "Seurat", "matrix", "Matrix")
  if (!any(unlist((lapply(allowed_formats, function(arg) is(data, arg)))))) {
    stop(
      "The data provided should be in one of the following formats: ",
      paste(allowed_formats, collapse = "|")
    )
  }

  if (is.character(signatures_input)) {
    signatures <- load_signatures(signatures_input)
  } else if (is.list(signatures_input)) {
    signatures <- signatures_input
  }

  stopifnot(is.logical(return_score))
  stopifnot(length(return_score) == 1)

  stopifnot(is.character(seurat_assay))
  stopifnot(length(seurat_assay) == 1)
  stopifnot(is.character(seurat_layer))
  stopifnot(length(seurat_layer) == 1)
  stopifnot(is.character(sce_assay))
  stopifnot(length(sce_assay) == 1)

  stopifnot(is.character(score_mode))
  score_mode <- match.arg(score_mode, c("raw", "scaled", "log", "log2", "log10"))

  if (!is.null(n_cpus)) {
    stopifnot(is.numeric(n_cpus))
    stopifnot(n_cpus > 0)
    stopifnot(length(n_cpus) == 1)
  }

  # Check the type of data and extract expression matrix accordingly


  if (is(data, "Seurat")) {
    stopifnot(is(data[[seurat_assay]], "Assay") | is(data[[seurat_assay]], "Assay5"))
    stopifnot(seurat_layer %in% Layers(data))

    obj_version_major <- as.numeric(
      strsplit(as.character(Version(data)), ".", fixed = TRUE)[[1]][1]
    )

    if (obj_version_major < 5) {
      message("Seurat object version is < 5.0.0 . It's suggested to run Seurat::UpdateSeuratObject first.")
      datam <- slot(data[[seurat_assay]], seurat_layer)
    } else {
      datam <- LayerData(data[[seurat_assay]], seurat_layer)
    }
  } else if (is(data, "SingleCellExperiment")) {
    stopifnot(sce_assay %in% assayNames(data))

    datam <- assay(data, sce_assay)
  } else if (is(data, "matrix") | is(data, "Matrix") | is(data, "DelayedMatrix")) {
    stopifnot(!is.null(dim(data)) & all(dim(data) > 0))
    if (is(data, "matrix")) {
      stopifnot(is.numeric(data))
    }
    if (is(data, "Matrix")) {
      stopifnot(is.numeric(data[1, 1]))
    }

    datam <- data
  }

  message("Checking if genes are in the dataset matrix...", "\n")
  result <- sapply(names(signatures), function(x) {
    lt <- length(signatures[[x]])
    l <- sum(signatures[[x]] %in% rownames(data))
    message <- sprintf("%s: %d / %d", x, l, lt)
    message(message)
    invisible(message) # Use invisible to avoid printing the return value of cat
  })

  n_cpus <- if (is.null(n_cpus)) {
    ceiling(parallel::detectCores() / 4)
  } else {
    n_cpus
  }
  tot <- colSums2(datam)
  # plan(multicore, workers = n_cpus)

  scores <- BiocParallel::bplapply(
    names(signatures),
    function(name) {
      geneset <- signatures[[name]]
      score_signature(
        data = datam,
        geneset = geneset,
        total_col_sums = tot
      )
    },
    BPPARAM = BiocParallel::MulticoreParam(n_cpus)
  )

  scores_df <- do.call(cbind, scores)
  colnames(scores_df) <- names(signatures)

  if (score_mode == "scaled") {
    scores_df <- sweep(scores_df, 2, apply(scores_df, 2, max), `/`)
  } else if (score_mode %in% c("log", "log2", "log10")) {
    scores_df <- match.fun(score_mode)(scores_df + .Machine$double.xmin)
  } else if (score_mode != "raw") {
    stop("Invalid score_mode. Must be one of: ['raw', 'scaled', 'log', 'log2', 'log10']")
  }


  if (return_score) {
    return(scores_df)
  }


  # TODO: call the column names in a clever way to make people notice where
  # these come from

  # ideally: CIA_signatureid_celltype
  # CIA_azimuth_Vip, for example


  if (is(data, "Seurat")) {
    data@meta.data[, colnames(scores_df)] <- scores_df

    message("Scores have been added in data@meta.data", "\n")

    return(data)
  } else if (is(data, "SingleCellExperiment")) {
    colData(data)[, colnames(scores_df)] <- scores_df

    message("Scores have been added in colData(data)", "\n")

    return(data)
  } else {
    return(scores_df)
  }
}


#' Classification Based on Signature Scores
#'
#' Performs classification of cells in `data` based on the computed signature
#' scores. Each sample/cell is classified to the signature with the highest score,
#' provided the difference between the top two scores exceeds a given threshold,
#' otherwise, it is marked as unassigned.
#'
#' @param data A numeric matrix, data frame, SeuratObject, or
#' SingleCellExperiment with genes in rows and cells in columns,
#' representing the expression level of each gene in each cell.
#' @param signatures_input A character string representing the file path to the
#' TSV file containing the gene signatures, or a list where each element is a
#' vector of gene names with names of the list elements representing signature names.
#' @param seurat_assay Assay for Seurat objects (default "RNA").
#' @param seurat_layer Character string, indicating which layer to use of the
#' Seurat object. Defaults to `data`.
#' @param sce_assay Character string, indicating which assay to use of the
#' SingleCellExperiment object. Defaults to `logcounts`.
#' @param score_mode Character string, specifying the calculation mode to be
#' used for scores: "raw", "scaled", or log-transformed ("log", "log2", "log10").
#' In the scope of this function, this defaults to "scaled".
#' @param n_cpus An optional integer indicating the number of CPU cores to use
#' for parallel computation. If NULL, the function will decide the number based on
#' available system resources.
#' @param similarity_threshold A numeric threshold used to decide if the highest
#' score is significantly higher than the second-highest to
#' assign a specific class to a cell. If the difference is less than or equal to
#' this threshold, the cell is labeled as unassigned.
#' @param column_name The name of the column to be added to the metadata storing
#' the classification labels, default is "CIA_prediction".
#' @param unassigned_label The label to assign to the cells where no clear majority
#' signature is identified, default is "Unassigned".
#'
#' @return Modifies the input data object by adding the classification labels to
#' its metadata and returns the modified data object.
#' If input data is a matrix or data frame, returns a vector of classification labels.
#'
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' ## TODO example
#'
#' ## A thought: rename it to CIA_classify? it is somehow catchy?
CIA_classify <- function(data,
                         signatures_input,
                         n_cpus = NULL,
                         similarity_threshold = 0.1,
                         seurat_assay = "RNA",
                         seurat_layer = "data",
                         sce_assay = "logcounts",
                         score_mode = "scaled",
                         column_name = "CIA_prediction",
                         unassigned_label = "Unassigned") {
  ## Checks on arguments
  allowed_formats <- c("SingleCellExperiment", "Seurat", "matrix", "Matrix")
  if (!any(unlist((lapply(allowed_formats, function(arg) is(data, arg)))))) {
    stop(
      "The data provided should be in one of the following formats: ",
      paste(allowed_formats, collapse = "|")
    )
  }

  if (is.character(signatures_input)) {
    signatures <- load_signatures(signatures_input)
  } else if (is.list(signatures_input)) {
    signatures <- signatures_input
  }

  if (!is.null(n_cpus)) {
    stopifnot(is.numeric(n_cpus))
    stopifnot(n_cpus > 0)
    stopifnot(length(n_cpus) == 1)
  }

  stopifnot(is.numeric(similarity_threshold))
  stopifnot(similarity_threshold >= 0 & similarity_threshold <= 1)

  stopifnot(is.character(seurat_assay))
  stopifnot(length(seurat_assay) == 1)
  stopifnot(is.character(seurat_layer))
  stopifnot(length(seurat_layer) == 1)

  stopifnot(is.character(sce_assay))
  stopifnot(length(sce_assay) == 1)

  stopifnot(is.character(score_mode))
  score_mode <- match.arg(score_mode, c("raw", "scaled", "log", "log2", "log10"))

  stopifnot(is.character(column_name))
  stopifnot(length(column_name) > 0)

  stopifnot(is.character(unassigned_label))

  start_time <- Sys.time() # Capture start time

  score_matrix <- score_all_signatures(data,
    signatures_input,
    score_mode = score_mode,
    seurat_assay = seurat_assay,
    seurat_layer = seurat_layer,
    sce_assay = sce_assay,
    return_score = TRUE,
    n_cpus = n_cpus
  )

  labels <- apply(score_matrix, 1,
    get_label,
    similarity_threshold = similarity_threshold,
    unassigned_label = unassigned_label
  )

  end_time <- Sys.time() # Capture end time

  # Format and print the message correctly
  message(
    "\nClassification complete!    Start:", format(start_time, "%H:%M:%S"),
    "    End:", format(end_time, "%H:%M:%S"), "\n"
  )

  if (is(data, "Seurat")) {
    data@meta.data[, column_name] <- as.factor(labels)

    message(column_name, "has been added in data@meta.data", "\n")

    return(data)
  } else if (is(data, "SingleCellExperiment")) {
    colData(data)[, column_name] <- as.data.frame(labels)

    message(column_name, "has been added in colData(data)", "\n")

    return(data)
  } else {
    return(as.factor(labels))
  }
}

#' Label extraction from the CIA score matrix
#'
#' @param row The vector of values, by row
#' @param similarity_threshold Numeric vector, to control the behavior of the
#' assignment of the labels.
#' @param unassigned_label The label to assign to the cells where no clear majority
#' signature is identified, default is "Unassigned" - handled via `CIA_classify`.
#'
#' @return The string label for the cell of interest, defined with the simple
#' heuristics on the similarity threshold.
get_label <- function(row,
                      similarity_threshold,
                      unassigned_label) {
  order_scores <- order(row, decreasing = TRUE)
  if (row[order_scores[1]] - row[order_scores[2]] <= similarity_threshold) {
    return(unassigned_label)
  } else {
    return(names(row)[order_scores[1]])
  }
}
