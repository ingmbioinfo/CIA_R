% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIA_extra.R
\name{celltypist_majority_vote}
\alias{celltypist_majority_vote}
\title{Perform Majority Voting on Cell Type Classifications}
\usage{
celltypist_majority_vote(
  data,
  classification_obs,
  groups_obs = NULL,
  graph_name = "RNA_snn",
  min_prop = 0,
  unassigned_label = "Unassigned",
  leiden_res = 0.05
)
}
\arguments{
\item{data}{A \code{Seurat} or \code{SingleCellExperiment} object containing the
single-cell data.}

\item{classification_obs}{A character vector specifying the column names in
the metadata that contain the current cell type classifications.}

\item{groups_obs}{A character string specifying the column name in the
metadata that contains the grouping information. If \code{NULL}, over-clustering
will be performed.}

\item{graph_name}{A character string specifying the graph name used for
clustering in Seurat. Default is "RNA_snn".}

\item{min_prop}{A numeric value specifying the minimum proportion threshold
for assigning a majority vote label. Default is 0.}

\item{unassigned_label}{A character string specifying the label for groups
that do not meet the minimum proportion threshold. Default is "Unassigned".}

\item{leiden_res}{A numeric value specifying the resolution parameter for the
Leiden clustering algorithm. Default is 0.05.}
}
\value{
A \code{Seurat} or \code{SingleCellExperiment} object, with updated metadata
containing the majority vote labels (with the information stored in vectors
with suffix \verb{_majority_voting}).
}
\description{
A function that wraps Celltypist majority voting (DOI: 10.1126/science.abl5197).
Assigns cell group labels based on the majority voting of cell type
predictions within each group.
If no reference cell groups are provided, an over-clustering step is
performed using the Leiden algorithm.
}
\details{
The function can handle both \code{Seurat} and \code{SingleCellExperiment} objects.
If the grouping information (\code{groups_obs}) is not provided, the function will
perform over-clustering using the Leiden algorithm on the specified graph.
The majority vote label for each group is determined by the most frequent
cell type label within the group.
Groups that do not meet the \code{min_prop} threshold are assigned the
\code{unassigned_label}.
}
\examples{
# TODO we need an example to run fully
\dontrun{
# Example usage with Seurat object
seurat_obj <- celltypist_majority_vote(
  seurat_obj,
  classification_obs = "predicted_labels",
  graph_name = "RNA_snn",
  min_prop = 0.5
)

# Example usage with SingleCellExperiment object
sce_obj <- celltypist_majority_vote(
  sce_obj,
  classification_obs = "predicted_labels",
  min_prop = 0.5
)
}
}
\references{
DOI: 10.1126/science.abl5197
}
