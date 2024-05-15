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
#'
#' @importFrom dplyr group_by mutate ungroup
#'
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
    plot_data <- plot_data |>
      group_by(Cluster) |>
      mutate(Value = Count / sum(Count) * 100) |>
      ungroup()
    y_label <- "Percentage"
  } else if (plot_type == "raw") {
    plot_data <- plot_data |>
      group_by(Cluster) |>
      mutate(Value = Count) |>
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
#' Group Composition Heatmap
#'
#' This function generates a heatmap showing the percentage composition of each classification
#' within reference groups. It is used for visualizing how different classifications distribute
#' across predefined groups.
#'
#' @param data Data frame or matrix containing the classification and reference data.
#' @param classification_obs Character string, the name of the column in \code{data} that contains the classification data.
#' @param ref_obs Character string, the name of the column in \code{data} that contains the reference group data.
#' @param columns_order Optional; a character vector specifying the order of columns in the heatmap.
#' @param color_map Optional; character string specifying the color palette to use, defaults to 'Greens'.
#'                  It should be one of the color palettes supported by \code{RColorBrewer}.
#'
#' @return A ggplot object representing the heatmap.
#'
#' @examples
#' ## Assuming 'data' is a data frame with columns 'Group' and 'Category'
#' # group_composition(data, 'Category', 'Group')
#' # TODO
#'
#' @import ggplot2 tidyr RColorBrewer
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats reshape
#'
#' @export

group_composition <- function(data, classification_obs, ref_obs, columns_order=NULL, color_map='Greens') {

  # Extract data
  ref_data <- data[[ref_obs]]
  class_data <- data[[classification_obs]]

  # Create a contingency table
  contingency_table <- table(ref_data, class_data)

  # Convert counts to percentage
  percentage_table <- prop.table(contingency_table, margin = 1) * 100
  percentage_table <- round(percentage_table, 2)

  # Convert the table to a data frame for ggplot
  df <- as.data.frame.matrix(percentage_table)
  if (!is.null(columns_order)) {
    df <- df[, columns_order]
  }
    df$rowname <- rownames(df)
    df <- reshape(df, varying = list(names(df)[1:(dim(df)[2]-1)]), v.names = "Value",
                  timevar = "Column", times = names(df)[1:(dim(df)[2]-1)], direction = "long")

  # Setting up the color palette
  my_palette <- colorRampPalette(brewer.pal(9, color_map))  # Customize this part for different palettes

  # Plot heatmap
  p <- ggplot(df, aes(x = rowname, y = Column, fill = Value)) +
    geom_tile() +
    scale_fill_gradientn(name = "", colors = my_palette(100)) +  # Use gradient n for multiple color transitions
    scale_x_discrete(name = "") +
    scale_y_discrete(name = "")

  print(p)
}

#' Plot Grouped Distributions and Perform Statistical Tests
#'
#' This function plots a heatmap of median values for selected columns in a dataframe across cell groups and performs statistical tests to evaluate the differences in distributions. The Wilcoxon test checks if each group's signature score is significantly higher than others in the same group. The Mann-Whitney U test checks if each signature has the highest score values in the corresponding group compared to all other groups.
#'
#' @param data DataFrame containing the cell data.
#' @param columns_obs Character vector. Column names in the dataframe where the values of interest are stored.
#' @param ref_obs Character. Column name in the dataframe where the cell group labels are stored.
#' @param color_map Character. Colormap for the heatmap. Defaults to 'Reds'.
#' @param scale_medians Character. How to scale the median values in the heatmap. Options: 'row-wise', 'column-wise', or NULL. Defaults to NULL.
#' @param save Character. Filename to save the heatmap. If provided, saves the heatmap in 'figures' directory with 'CIA_' prefix. Defaults to NULL.
#' 
#' @return None. The function either saves the heatmap to a file or prints it.
#' 
#' @examples
#' grouped_distributions(data, columns_obs = c('feature1', 'feature2'), ref_obs = 'group_column')
#' 
#' @import ggplot2
#' @import reshape2
#' @importFrom stats wilcox.test
grouped_distributions <- function(data, columns_obs, ref_obs, color_map = 'Reds', scale_medians = NULL, save = NULL) {
  # Compute median values for each group
  unique_groups <- levels(data[[ref_obs]])
  grouped_df <- do.call(rbind, lapply(unique_groups, function(group) {
    group_data <- data[data[[ref_obs]] == group, columns_obs, drop = FALSE]
    medians <- apply(group_data, 2, median, na.rm = TRUE)
    c(group, medians)
  }))
  colnames(grouped_df) <- c(ref_obs, columns_obs)
  grouped_df <- as.data.frame(grouped_df)
  
  if (!is.null(scale_medians)) {
    if (scale_medians == 'row-wise') {
      grouped_df[columns_obs] <- t(apply(grouped_df[columns_obs], 1, function(x) x / sum(x, na.rm = TRUE)))
    } else if (scale_medians == 'column-wise') {
      grouped_df[columns_obs] <- sweep(grouped_df[columns_obs], 2, colSums(grouped_df[columns_obs], na.rm = TRUE), `/`)
    }
  }
  # Convert to long format for plotting
  df_long <- melt(grouped_df, id.vars = ref_obs, variable.name = "Column", value.name = "Value")
  # Wilcoxon test
  print('Performing Wilcoxon test on each cell group ...')
  count <- 0
  subsets <- list()
  
  for (group in unique_groups) {
    subset_data <- data[data[[ref_obs]] == group, columns_obs, drop = FALSE]
    subsets[[group]] <- subset_data
    medians <- apply(subset_data, 2, median, na.rm = TRUE)
    pos <- which.max(medians)
    
    combs <- combn(columns_obs, 2, simplify = FALSE)
    for (comb in combs) {
      if (all(subset_data[[comb[1]]] == 0) || all(subset_data[[comb[2]]] == 0)) {
        next
      }
      result <- tryCatch({
        wilcox.test(subset_data[[comb[1]]], subset_data[[comb[2]]], paired = TRUE, alternative = "two.sided")
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(result) && !is.na(result$p.value) && (result$p.value >= 0.01) && (comb[1] == names(medians)[pos])) {
        count <- count + 1
        print(paste('WARNING in cell group', group, ':', comb[1], 'values are not significantly different from', comb[2], 'values.'))
        print(paste('(p=', result$p.value, ')'))
      }
    }
  }
  
  if (count == 0) {
    print('For each cell group there is a distribution significantly higher than the others (p<0.01)')
  }
  
  # Mann-Whitney U test
  print('Performing Mann-Whitney U test on each selected AnnData.obs column ...')
  count <- 0
  for (column in columns_obs) {
    sign <- lapply(unique_groups, function(group) data[data[[ref_obs]] == group, column])
    names(sign) <- unique_groups
    sign_medians <- sapply(sign, median, na.rm = TRUE)
    pos <- which.max(sign_medians)
    
    combs <- combn(unique_groups, 2, simplify = FALSE)
    for (comb in combs) {
      group1_values <- sign[[comb[1]]]
      group2_values <- sign[[comb[2]]]
      
      if (all(is.na(group1_values)) || all(is.na(group2_values)) || all(group1_values == 0) || all(group2_values == 0)) {
        next
      }
      
      result <- tryCatch({
        wilcox.test(group1_values, group2_values, alternative = "two.sided", paired = FALSE)
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(result) && !is.na(result$p.value) && (result$p.value >= 0.01) && (comb[1] == names(sign_medians)[pos])) {
        count <- count + 1
        print(paste('WARNING in', column, 'distribution: values in', comb[1], 'group are not significantly different from values in', comb[2], 'group'))
        print(paste('(p=', result$p.value, ')'))
      }
    }
  }
  
  if (count == 0) {
    print('For each distribution, there is only a cell group in which values are higher with respect to all the other groups (p<0.01)')
  }
  
  # Plotting
# If the columns are not already factors, convert them
df_long$ref <- as.factor(df_long[,ref_obs])
df_long$Column <- as.factor(df_long$Column)
df_long$Value <- as.numeric(df_long$Value)

  p <- ggplot(df_long, aes(x = ref, y = Column, fill = Value)) +
    geom_tile() +
    scale_fill_gradientn(name = "", colors = colorRampPalette(brewer.pal(9, color_map))(100)) +
    theme_minimal() +
    labs(title = "Median score values", x = ref_obs, y = "Signatures")
  
  # Save or display plot
  if (!is.null(save)) {
    if (!dir.exists("figures")) dir.create("figures")
    ggsave(paste0("figures/CIA_", save), plot = p, width = 10, height = 8)
  } else {
    print(p)
  }
}
