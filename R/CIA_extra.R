#' Perform Majority Voting on Cell Type Classifications
#'
#' A function that wraps Celltypist majority voting (DOI: 10.1126/science.abl5197). 
#' Assigns cell group labels based on the majority voting of cell type predictions within each group.
#' If no reference cell groups are provided, an over-clustering step is performed using the Leiden algorithm.
#'
#' @param data A Seurat or SingleCellExperiment object containing the single-cell data.
#' @param classification_obs A character vector specifying the column names in the metadata that contain cell type classifications.
#' @param groups_obs A character string specifying the column name in the metadata that contains the grouping information. If NULL, over-clustering will be performed.
#' @param graph A character string specifying the graph name used for clustering in Seurat. Default is "RNA_snn".
#' @param min_prop A numeric value specifying the minimum proportion threshold for assigning a majority vote label. Default is 0.
#' @param unassigned_label A character string specifying the label for groups that do not meet the minimum proportion threshold. Default is 'Unassigned'.
#' @param res A numeric value specifying the resolution parameter for the Leiden clustering algorithm. Default is 0.05.
#' 
#' @return A Seurat or SingleCellExperiment object with updated metadata containing the majority vote labels.
#'
#' @details
#' The function can handle both Seurat and SingleCellExperiment objects. If the grouping information (`groups_obs`) is not provided, the function will perform over-clustering using the Leiden algorithm on the specified graph. The majority vote label for each group is determined by the most frequent cell type label within the group. Groups that do not meet the `min_prop` threshold are assigned the `unassigned_label`.
#'
#' @examples
#' \dontrun{
#' # Example usage with Seurat object
#' seurat_obj <- celltypist_majority_vote(seurat_obj, classification_obs = "predicted_labels", graph = "RNA_snn", min_prop = 0.5)
#'
#' # Example usage with SingleCellExperiment object
#' sce_obj <- celltypist_majority_vote(sce_obj, classification_obs = "predicted_labels", min_prop = 0.5)
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom igraph graph_from_adjacency_matrix cluster_leiden
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom methods is
#'
#' @export
celltypist_majority_vote <- function(data, classification_obs, groups_obs = NULL, graph = "RNA_snn", min_prop = 0, unassigned_label = 'Unassigned', res=0.05) {
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
      if (!(graph %in% names(data@graphs))) {
        stop(paste("Provided graph name", graph, "not present in Seurat object. Please run FindNeighbors() first."))
      }
       adj<- as.matrix(data@graphs[[graph]])
       gr <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
       leid<- cluster_leiden(
              gr,
              weights = NULL,
              resolution_parameter = res,
              beta = 0.01,
              n_iterations = -1,
              vertex_weights = NULL
            )
      data@meta.data$overclustering <- leid$membership
      groups_obs<- 'overclustering'
      obs<- data@meta.data
    } else if (is(data, "SingleCellExperiment")) {
         if (!(graph %in% names(metadata(data)$graphs))) {
        stop(paste("Provided graph name", graph, "not present in Seurat object. Please run FindNeighbors() first."))
      }
   adj<- as.matrix(data@metadata$graphs[[graph]])
       gr <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
       leid<- cluster_leiden(
              gr,
              weights = NULL,
              resolution_parameter = res,
              beta = 0.01,
              n_iterations = -1,
              vertex_weights = NULL
            )
      colData(data)$overclustering <- leid$membership
      groups_obs<- 'overclustering'
      obs <- colData(data)
    }
    message("Reference annotation not selected. Computing over-clustering with Leiden algorithm")
  }

  # Check if groups_obs column exists
  if (!(groups_obs %in% colnames(obs))) {
    stop(paste("Column", groups_obs, "not found in metadata. Please check the clustering step."))
  }

  message(paste("Dataset has been divided into", length(unique(obs[[groups_obs]])), "groups according to transcriptional similarities."))
  message(paste0("Over-clustering result saved in metadata as ", groups_obs))

  message("Extending the more represented cell type label to each cell group...\n")

  groups <- obs[[groups_obs]]

  #Check if classification_obs columns exist
  missing_columns <- classification_obs[!(classification_obs %in% colnames(obs))]
  if (length(missing_columns) > 0) {
    stop(paste("The following classification_obs columns are missing in metadata:", paste(missing_columns, collapse = ", ")))
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
    colData(data) <- S4Vectors::DataFrame(obs)
  }

  return(data)
}

# Example usage with Seurat object
# seurat_obj <- celltypist_majority_vote(seurat_obj, classification_obs = "predicted_labels", graph = "RNA_snn", min_prop = 0.5)

# Example usage with SingleCellExperiment object
# sce_obj <- celltypist_majority_vote(sce_obj, classification_obs = "predicted_labels", min_prop = 0.5)
