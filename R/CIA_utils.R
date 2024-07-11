#' Compute Similarity Between Gene Signatures
#'
#' Computes the similarity between gene signatures using either the Jaccard
#' index or the percentage of intersection.
#'
#' @param signatures_dict A list where the keys are signature names and the
#' values are lists of gene names (gene signatures).
#' @param metric_name A character string specifying the metric for showing
#' similarities: "jaccard" for Jaccard index or "percentage" for percentages of
#' intersection. Defaults to "jaccard".
#'
#' @return A data frame containing the similarity of each pair of signatures,
#' with signatures as both rows and columns.
#'
#' @export
#'
#' @examples
#' signatures <- list(
#'   "signature1" = c("gene1", "gene2", "gene3"),
#'   "signature2" = c("gene2", "gene3", "gene4"),
#'   "signature3" = c("gene1", "gene5")
#' )
#' similarity <- signatures_similarity(signatures, metric_name = "jaccard")
#' print(similarity)
signatures_similarity <- function(signatures_dict,
                                  metric_name = "jaccard") {
  metric_name <- match.arg(metric_name, c("jaccard", "percentage"))

  signature_names <- names(signatures_dict)
  n <- length(signature_names)
  similarity_matrix <- matrix(0, n, n,
                              dimnames = list(signature_names, signature_names))

  for (i in 1:n) {
    for (j in i:n) {
      intersec <- length(
        intersect(
          signatures_dict[[signature_names[i]]],
          signatures_dict[[signature_names[j]]]
        )
      )
      if (metric_name == "jaccard") {
        union <- length(
          union(
            signatures_dict[[signature_names[i]]],
            signatures_dict[[signature_names[j]]]
          )
        )
        similarity <- intersec / union
      } else if (metric_name == "percentage") {
        similarity <-
          round(100 * intersec / length(signatures_dict[[signature_names[i]]]), 2)
      }
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }

  similarity <- as.data.frame(similarity_matrix)

  return(similarity)
}
