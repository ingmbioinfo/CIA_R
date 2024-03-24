#' Compute Classification Metrics
#'
#' @description
#' Computes the main metrics of classification by comparing the labels of cells classified
#' by given methods to labels assigned by a reference method. Cells labeled as `unassigned_label`
#' are excluded from the metrics calculation. Additionally, the percentage of unassigned cells
#' for each classification method is reported.
#'
#' @param data A data frame or matrix representing the cell data, typically an AnnData object in R.
#' @param classification_col A vector of strings specifying the column names in `data` where the
#' labels assigned by the methods of interest are stored.
#' @param groups_col A string specifying the column name in `data` where the labels assigned by
#' the reference method are stored.
#' @param unassigned_label A string representing the label used to mark unassigned cells in the
#' classification columns. Default is an empty string, which means no cells are excluded based on
#' their label.
#'
#' @return A data frame containing the overall sensitivity (SE), specificity (SP), precision (PR),
#' accuracy (ACC), F1-score (F1), and percentage of unassigned cells (%UN) for each classification
#' method compared to the reference method.
#'
#' @export
#'
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows
#'
#' @examples
#' ## TODO
classification_metrics <- function(data, classification_col, groups_col, unassigned_label = '') {
  report <- purrr::map_dfr(classification_col, function(m) {
    total_cells <- nrow(data)
    unassigned_count <- sum(data[[m]] == unassigned_label)
    percent_unassigned <- (unassigned_count / total_cells) * 100

    filtered_data <- data[data[[m]] != unassigned_label, ]

    metrics <- lapply(unique(filtered_data[[groups_col]]), function(group) {
      TP <- sum(filtered_data[[m]] == group & filtered_data[[groups_col]] == group)
      TN <- sum(filtered_data[[m]] != group & filtered_data[[groups_col]] != group)
      FP <- sum(filtered_data[[m]] == group & filtered_data[[groups_col]] != group)
      FN <- sum(filtered_data[[m]] != group & filtered_data[[groups_col]] == group)

      SE <- TP / (TP + FN)
      SP <- TN / (TN + FP)
      PR <- TP / (TP + FP)
      ACC <- (TP + TN) / (TP + TN + FP + FN)
      F1 <- 2 * TP / (2 * TP + FP + FN)

      c(SE = SE, SP = SP, PR = PR, ACC = ACC, F1 = F1)
    }) %>% dplyr::bind_rows(.f = mean)

    c(metrics, `%UN` = percent_unassigned)
  }, .id = "Method")

  rownames(report) <- classification_col
  return(report)
}

#' Compute Classification Metrics
#'
#' Calculates classification metrics such as sensitivity (SE), specificity (SP),
#' precision (PR), accuracy (ACC), F1-score (F1), and the percentage of 'Unassigned' cells (%UN)
#' for given classification results, based on the specified unassigned label.
#'
#' @param data Data frame containing the observed data including classification and true group labels.
#' @param classification_col Character vector of column names in `data` that contain the classification results to evaluate.
#' @param groups_cols Character string specifying the column in `data` that contains the true group labels.
#' @param unassigned_label The label used in `classification_col` to denote unclassified or unassigned samples.
#' If this label is present, %UN will be calculated and included in the metrics.
#'
#' @return A matrix with classification metrics (SE, SP, PR, ACC, F1) for each classifier.
#' If unassigned_label is specified and present in the data, an additional '%UN' column is included in the output.
#'
#' @export
#'
#' @examples
#' ## TODO: probably have a pre-computed set of labels
#' ## TODO: define-describe the process to get there in inst/scripts or so?
#'
#' # Assuming data is a data frame with true labels and classification results
#' data <- read.csv('your_data_file.csv')
#' classification_col <- c('classifier1', 'classifier2')
#' groups_cols <- 'true_group'
#' unassigned_label <- 'Unassigned' # Specify the label that denotes unassigned samples
#' metrics_report <- compute_metrics(data, classification_col, groups_cols, unassigned_label)
#' print(metrics_report)
compute_metrics <- function(data, classification_col, groups_cols, unassigned_label = '') {
  report <- list()

  for (m in classification_col) {
    TP_l <- TN_l <- FP_l <- FN_l <- numeric()
    UN <- round(100 * sum(data[[m]] == unassigned_label) / nrow(data), 2)
    datam <- data[data[,m] != unassigned_label,]

    for (i in unique(datam[[groups_cols]])) {
      TP_l <- c(TP_l, sum(datam[[m]] == i & datam[[groups_cols]] == i))
      TN_l <- c(TN_l, sum(datam[[m]] != i & datam[[groups_cols]] != i))
      FP_l <- c(FP_l, sum(datam[[m]] == i & datam[[groups_cols]] != i))
      FN_l <- c(FN_l, sum(datam[[m]] != i & datam[[groups_cols]] == i))
    }

    TP <- sum(TP_l)
    TN <- sum(TN_l)
    FP <- sum(FP_l)
    FN <- sum(FN_l)

    SE <- TP / (TP + FN)
    SP <- TN / (TN + FP)
    PR <- TP / (TP + FP)
    ACC <- (TP + TN) / (TP + TN + FP + FN)
    F1 <- 2 * TP / (2 * TP + FP + FN)

    report[[m]] <- c(SE, SP, PR, ACC, F1, UN)
  }

  report <- t(as.data.frame(report))
  colnames(report) <- c('SE', 'SP', 'PR', 'ACC', 'F1', '%UN')
  if (sum(report[,'%UN'])==0){
    report <- report[,1:5]
  }
  return(report)
}
