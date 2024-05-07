#' Compute Classification Metrics
#'
#' Calculates classification metrics such as sensitivity (SE), specificity (SP),
#' precision (PR), accuracy (ACC), F1-score (F1), and the percentage of
#' 'Unassigned' cells (%UN) for given classification results, based on the
#' specified unassigned label. Cells labeled as specified by the value of
#' `unassigned_label` can be excluded from the metrics calculation.
#' Additionally, the percentage of unassigned cells for each classification
#' method is reported.
#'
#' @param cells_info Data frame containing the observed data including the labels
#' assigned by the classification algorithm(s) and the true group labels.
#' Typically, this is the output of `colData()` for a SingleCellExperiment
#' object, and the `meta.data` slot for a Seurat object.
#' @param classification_cols A vector of character strings specifying the column
#' names in `cells_info` where the labels assigned by the methods of interest are
#' stored.
#' @param ref_labels Character string specifying the column in `cells_info` that
#' contains the true group labels (or the ones assigned by the reference method).
#' @param unassigned_label The label used in `classification_cols` to denote
#' unclassified or unassigned samples.
#' If this parameter is specified and detected in the `cells_info` data, %UN
#' will be calculated and included in the output metrics.
#'
#' @return A matrix with classification metrics (SE, SP, PR, ACC, F1) for each
#' classifier. If unassigned_label is specified and present in the `cells_info```
#' data, an additional '%UN' column is included in the output.
#'
#' @details
#' The classification metrics reported in the output object are;
#' * the overall sensitivity (SE)
#' * the specificity (SP)
#' * the precision (PR)
#' * the accuracy (ACC)
#' * the F1-score (F1)
#' * the percentage of unassigned cells (%UN) for each classification method,
#' compared to the reference method.
#'
#' @export
#'
#' @examples
#'
#' ## TODO example
#'
#' ## # Assuming cells_info is a data frame with true labels and classification results
#' ## cells_info <- read.csv('your_data_file.csv')
#' ## classification_cols <- c('classifier1', 'classifier2')
#' ## ref_labels <- 'true_group'
#' ## unassigned_label <- 'Unassigned' # Specify the label that denotes unassigned samples
#' ## metrics_report <-
#' ##   compute_classification_metrics(cells_info, classification_cols, ref_labels, unassigned_label)
#' ## print(metrics_report)
compute_classification_metrics <- function(cells_info,
                                           classification_cols,
                                           ref_labels,
                                           unassigned_label = "") {
  stopifnot(is(cells_info, "DataFrame") || is.data.frame(cells_info))
  stopifnot(is.character(classification_cols),
            length(classification_cols) > 0)
  stopifnot(is.character(ref_labels),
            length(ref_labels) == 1)
  stopifnot(is.character(unassigned_label),
            length(unassigned_label) == 1)

  if (!all(classification_cols %in% colnames(cells_info)))
    stop("Some column names for the classification were not found in the `cells_info`")

  if (!ref_labels %in% colnames(cells_info))
    stop("The value specified for `ref_labels` was not found as a column name in the `cells_info`")

  report <- list()

  for (m in classification_cols) {
    TP_l <- TN_l <- FP_l <- FN_l <- numeric()
    UN <- round(100 * sum(cells_info[[m]] == unassigned_label) / nrow(cells_info), 2)
    datam <- cells_info[cells_info[, m] != unassigned_label, ]

    for (i in unique(datam[[ref_labels]])) {
      TP_l <- c(TP_l, sum(datam[[m]] == i & datam[[ref_labels]] == i, na.rm = TRUE))
      TN_l <- c(TN_l, sum(datam[[m]] != i & datam[[ref_labels]] != i, na.rm = TRUE))
      FP_l <- c(FP_l, sum(datam[[m]] == i & datam[[ref_labels]] != i, na.rm = TRUE))
      FN_l <- c(FN_l, sum(datam[[m]] != i & datam[[ref_labels]] == i, na.rm = TRUE))
    }

    TP <- sum(TP_l, na.rm = TRUE)
    TN <- sum(TN_l, na.rm = TRUE)
    FP <- sum(FP_l, na.rm = TRUE)
    FN <- sum(FN_l, na.rm = TRUE)

    SE <- TP / (TP + FN)
    SP <- TN / (TN + FP)
    PR <- TP / (TP + FP)
    ACC <- (TP + TN) / (TP + TN + FP + FN)
    F1 <- 2 * TP / (2 * TP + FP + FN)

    report[[m]] <- c(SE, SP, PR, ACC, F1, UN)
  }

  report <- t(as.data.frame(report))
  colnames(report) <- c("SE", "SP", "PR", "ACC", "F1", "%UN")
  if (sum(report[, "%UN"]) == 0) {
    report <- report[, 1:5, drop = FALSE]
  }
  return(report)
}
