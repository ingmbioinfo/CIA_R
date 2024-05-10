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

#' Compute Classification Metrics Per Reference Group
#'
#' Calculates classification metrics such as sensitivity (SE), specificity (SP),
#' precision (PR), accuracy (ACC), F1-score (F1), and the percentage of
#' 'Unassigned' cells (%UN),based on the
#' specified unassigned label, for each reference group. Cells labeled as specified by the value of
#' `unassigned_label` can be excluded from the metrics calculation.
#' Additionally, the percentage of unassigned cells for each classification
#' method is reported.
#'
#' @param cells_info Data frame containing the observed data including the labels
#' assigned by the classification algorithm(s) and the true group labels.
#' Typically, this is the output of `colData()` for a SingleCellExperiment
#' object, and the `meta.data` slot for a Seurat object.
#' @param classification_col A strings specifying the column
#' name in `cells_info` where the labels assigned by the method of interest are
#' stored.
#' @param ref_labels Character string specifying the column in `cells_info` that
#' contains the true group labels (or the ones assigned by the reference method).
#' @param unassigned_label The label used in `classification_col` to denote
#' unclassified or unassigned samples.
#' If this parameter is specified and detected in the `cells_info` data, %UN
#' will be calculated and included in the output metrics.
#'
#' @return A matrix with classification metrics (SE, SP, PR, ACC, F1) for each
#' reference group. If unassigned_label is specified and present in the `cells_info```
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
grouped_classification_metrics <- function(cells_info,
                                           classification_col,
                                           ref_labels,
                                           unassigned_label = "") {
  stopifnot(is(cells_info, "DataFrame") || is.data.frame(cells_info))
  stopifnot(is.character(classification_col),
            length(classification_col) > 0)
  stopifnot(is.character(ref_labels),
            length(ref_labels) == 1)
  stopifnot(is.character(unassigned_label),
            length(unassigned_label) == 1)

  if (!all(classification_col %in% colnames(cells_info)))
    stop("Some column names for the classification were not found in the `cells_info`")

  if (!ref_labels %in% colnames(cells_info))
    stop("The value specified for `ref_labels` was not found as a column name in the `cells_info`")


    TP_l <- TN_l <- FP_l <- FN_l <- UN_l<- numeric()
    datam <- cells_info[cells_info[, classification_col] != unassigned_label, ]

    for (i in unique(datam[[ref_labels]])) {
      TP_l <- c(TP_l, sum(datam[[classification_col]] == i  & datam[[ref_labels]] == i, na.rm = TRUE))
      TN_l <- c(TN_l, sum(datam[[classification_col]] != i & datam[[ref_labels]] != i, na.rm = TRUE))
      FP_l <- c(FP_l, sum(datam[[classification_col]] == i & datam[[ref_labels]] != i, na.rm = TRUE))
      FN_l <- c(FN_l, sum(datam[[classification_col]] != i & datam[[ref_labels]] == i, na.rm = TRUE))
      UN_l <- c(UN_l, round(100 * sum(cells_info[[classification_col]]==unassigned_label & cells_info[[ref_labels]]==i)
                    / sum(cells_info[[ref_labels]]==i), 2))
    }



    SE <- TP_l / (TP_l + FN_l)
    SP <- TN_l/ (TN_l + FP_l)
    PR <- TP_l / (TP_l + FP_l)
    ACC <- (TP_l + TN_l) / (TP_l + TN_l + FP_l + FN_l)
    F1 <- 2 * TP_l / (2 * TP_l + FP_l + FN_l)
  
       
  report <- cbind(SE, SP, PR, ACC, F1, UN_l)
  report <- as.data.frame(report)
  rownames(report)<- unique(datam[[ref_labels]])
  colnames(report) <- c("SE", "SP", "PR", "ACC", "F1", "%UN")
  if (sum(report[, "%UN"]) == 0) {
    report <- report[, 1:5, drop = FALSE]
  }
  return(report)
}

#' Plot Group Composition
#'
#' This function plots the composition of each reference group as a horizontal stacked bar plot. 
#' The composition can be shown either as raw counts or as percentages.
#'
#' @param df DataFrame containing the data to be plotted.
#' @param ref_col Character. The name of the column representing the reference grouping variable.
#' @param comp_col Character. The name of the column representing the grouping to be compared.
#' @param plot_type Character. Indicates whether to plot 'percentage' or 'raw' counts. Defaults to 'percentage'.
#' @param palette Character vector. Specifies the color palette to use. Defaults to RColorBrewer::brewer.pal(8, "Dark2").
#' @param show_legend Logical. Whether to display the legend on the plot. Defaults to TRUE.
#'
#' @return A horizontal stacked bar plot showing the composition of each reference group.
#' @export
plot_group_composition <- function(df, ref_col, comp_col,
                                   plot_type = "percentage", palette = RColorBrewer::brewer.pal(8, "Set3"),
                                   show_legend = TRUE) {
  # Check if the specified columns exist in the DataFrame
  if (!(comp_col %in% names(df) && ref_col %in% names(df))) {
    stop("Specified columns are not in the DataFrame")
  }

  # Create a contingency table of counts
  pivot_table <- table(df[[ref_col]], df[[comp_col]])
  
  # Convert to DataFrame for plotting
  plot_data <- as.data.frame.matrix(pivot_table)
  plot_data <- cbind(Cluster = rownames(plot_data), plot_data)
  plot_data <- tidyr::pivot_longer(plot_data, cols = -Cluster, names_to = "Group", values_to = "Count")

  # Calculate percentages if required
  if (plot_type == "percentage") {
    plot_data <- plot_data %>%
      group_by(Cluster) %>%
      mutate(Value = Count / sum(Count) * 100) %>%
      ungroup()
    y_label <- "Percentage"
  } else if (plot_type == "raw") {
    plot_data <- plot_data %>%
      group_by(Cluster) %>%
      mutate(Value = Count) %>%
      ungroup()
    y_label <- "Count"
  } else {
    stop("plot_type must be 'percentage' or 'raw'")
  }
  
  # Plotting
  p <- ggplot(plot_data, aes(x = Cluster, y = Value, fill = Group)) +
    geom_bar(stat = 'identity', position = 'stack') +
    coord_flip() +
    scale_fill_manual(values = palette) +
    labs(y = y_label, x = ref_col, fill = comp_col) +
    theme_minimal()

  # Add legend if required
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  # Show plot
  print(p)
}

