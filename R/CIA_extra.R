# TO BE TESTED
# TO BE TESTED
# TO BE TESTED
#' Celltypist Majority Vote
#'
#' Assigns cell group labels based on the majority voting of cell type predictions within each group.
#' If no reference cell groups are provided, an over-clustering step is performed using the Leiden algorithm.
#'
#' @param data A `Seurat` or `SingleCellExperiment` object containing the cell data.
#' @param classification_obs A string or list of strings specifying the column(s) in metadata where the cell type predictions (labels) are stored.
#' @param groups_obs A string specifying the column in metadata where the reference group labels are stored. If NULL, an over-clustering with the Leiden algorithm is performed based on the dataset size.
#' @param graph A string specifying the graph to use for clustering in `Seurat` objects. Default is "RNA_snn".
#' @param min_prop A numeric value specifying the minimum proportion of cells required to assign a majority vote label to a group. Default is 0.
#' @param unassigned_label A string specifying the label to assign to cell groups where no cell type reaches the minimum proportion. Default is 'Unassigned'.
#'
#' @return The modified `Seurat` or `SingleCellExperiment` object with new majority voting classification labels stored in metadata.
#'
#' @details
#' This function computes the majority voting for cell type labels within each group of cells.
#' If `groups_obs` is not provided, the function performs over-clustering using the Leiden algorithm with a resolution adjusted based on the dataset size.
#' The results of the majority voting are stored back in the metadata of the object, adding a column for each classification considered.
#'
#' @importFrom clusterExperiment ClusterExperiment
#' @importFrom Seurat FindClusters
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with Seurat object
#' seurat_obj <- celltypist_majority_vote(
#'   seurat_obj,
#'   classification_obs = "predicted_labels",
#'   graph = "RNA_snn",
#'   min_prop = 0.5
#' )
#'
#' # Example usage with SingleCellExperiment object
#' sce_obj <- celltypist_majority_vote(
#'   sce_obj,
#'   classification_obs = "predicted_labels",
#'   min_prop = 0.5
#' )
#' }
celltypist_majority_vote <- function(data, classification_obs, groups_obs = NULL, graph = "RNA_snn", min_prop = 0, unassigned_label = 'Unassigned') {
  if (!is(data, "Seurat") && !is(data, "SingleCellExperiment")) {
    stop("Data must be a Seurat or SingleCellExperiment object")
  }

  if (is(data, "Seurat")) {
    obs <- data@meta.data
  } else if (is(data, "SingleCellExperiment")) {
    obs <- as.data.frame(colData(data))
  }

  # Determine resolution for clustering if groups_obs is not provided
  if (is.null(groups_obs)) {
    if (is(data, "Seurat")) {
      resolution <- 5 + 5 * (nrow(data@meta.data) %/% 20000)
      # Check if the specified graph is available for clustering
      if (!(graph %in% names(data@graphs))) {
        stop(paste("Provided graph.name", graph, "not present in Seurat object. Please run FindNeighbors() first."))
      }
      data <- FindClusters(data, graph.name = graph, resolution = resolution)
      groups_obs <- paste0("seurat_clusters_", resolution)
    } else if (is(data, "SingleCellExperiment")) {
      resolution <- 5 + 5 * (nrow(colData(data)) %/% 20000)
      data <- clusterExperiment::ClusterExperiment(data, method = "Leiden", resolution = resolution)
      groups_obs <- paste0("leiden_", resolution)
    }
    message(paste("Reference annotation not selected. Computing over-clustering with Leiden algorithm (resolution=", resolution, ") ...", sep = ""))
    message(paste("Dataset has been divided into", length(unique(obs[[groups_obs]])), "groups according to transcriptional similarities."))
    message(paste0("Over-clustering result saved in metadata as leiden_", as.character(resolution) ))
  } else {
    message(paste(groups_obs, " in metadata has been selected as reference annotation."))
  }

  message("Extending the more represented cell type label to each cell group...\n")

  groups <- obs[[groups_obs]]

  # Ensure classification_obs is a list
  if (!is.list(classification_obs)) {
    classification_obs <- list(classification_obs)
  }

  for (classification in classification_obs) {
    votes <- table(obs[[classification]], groups)
    majority <- apply(votes, 2, function(x) names(which.max(x)))
    freqs <- apply(votes, 2, max) / colSums(votes)

    majority_labels <- ifelse(freqs >= min_prop, majority, unassigned_label)
    obs[[paste0(classification, "_majority_voting")]] <- factor(groups, levels = names(majority_labels), labels = majority_labels)

    message(paste0("New classification labels have been stored in metadata as " , classification, "_majority_voting."))
  }

  if (is(data, "Seurat")) {
    data@meta.data <- obs
  } else if (is(data, "SingleCellExperiment")) {
    colData(data) <- DataFrame(obs)
  }

  return(data)
}

# Example usage with Seurat object
# seurat_obj <- celltypist_majority_vote(seurat_obj, classification_obs = "predicted_labels", graph = "RNA_snn", min_prop = 0.5)

# Example usage with SingleCellExperiment object
# sce_obj <- celltypist_majority_vote(sce_obj, classification_obs = "predicted_labels", min_prop = 0.5)
