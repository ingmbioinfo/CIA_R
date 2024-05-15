#' Compute Similarity Between Gene Signatures
#'
#' Computes the similarity between gene signatures using either the Jaccard index or the percentage of intersection.
#'
#' @param signatures_dict A list where the keys are signature names and the values are lists of gene names (gene signatures).
#' @param show A character string specifying the metric for showing similarities: 'J' for Jaccard index or '%' for percentages of intersection. Default is 'J'.
#' 
#' @return A data frame containing the similarity of each pair of signatures, with signatures as both rows and columns.
#' 
#' @examples
#' signatures <- list(
#'   'signature1' = c('gene1', 'gene2', 'gene3'),
#'   'signature2' = c('gene2', 'gene3', 'gene4'),
#'   'signature3' = c('gene1', 'gene5')
#' )
#' similarity <- signatures_similarity(signatures, show = 'J')
#' print(similarity)
#' 
#' @export
signatures_similarity <- function(signatures_dict, show = 'J') {
  if (!show %in% c('J', '%')) {
    stop('show must be "J" or "%".')
  }
  
  signature_names <- names(signatures_dict)
  n <- length(signature_names)
  similarity_matrix <- matrix(0, n, n, dimnames = list(signature_names, signature_names))
  
  for (i in 1:n) {
    for (j in i:n) {
      intersec <- length(intersect(signatures_dict[[signature_names[i]]], signatures_dict[[signature_names[j]]]))
      if (show == 'J') {
        union <- length(union(signatures_dict[[signature_names[i]]], signatures_dict[[signature_names[j]]]))
        similarity <- intersec / union
      } else if (show == '%') {
        similarity <- round(100 * intersec / length(signatures_dict[[signature_names[i]]]), 2)
      }
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }
  
  similarity <- as.data.frame(similarity_matrix)
  return(similarity)
}
