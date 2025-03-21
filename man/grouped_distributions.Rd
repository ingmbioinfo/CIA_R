% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIA_metrics.R
\name{grouped_distributions}
\alias{grouped_distributions}
\title{Plot Grouped Distributions and Perform Statistical Tests}
\usage{
grouped_distributions(
  data,
  columns_obs,
  ref_obs,
  color_map = "Reds",
  scale_medians = NULL,
  save = NULL
)
}
\arguments{
\item{data}{A \code{data.frame} object containing the cell data.}

\item{columns_obs}{Character vector. Column names in the data.frame where the
values of interest are stored.}

\item{ref_obs}{Character. Column name in the data.frame where the cell group
labels are stored.}

\item{color_map}{Character. Colormap for the heatmap. Defaults to "Reds".}

\item{scale_medians}{Character. How to scale the median values in the heatmap.
Options: 'row-wise', 'column-wise', or NULL. Defaults to NULL.}

\item{save}{Character. Filename to save the heatmap. If provided, saves the
heatmap in 'figures' directory with 'CIA_' prefix. Defaults to NULL.}
}
\value{
None. The function either saves the heatmap to a file or prints it.
}
\description{
This function plots a heatmap of median values for selected columns in a
dataframe across cell groups and performs statistical tests to evaluate the
differences in distributions. The Wilcoxon test checks if each group's
signature score is significantly higher than others in the same group.
The Mann-Whitney U test checks if each signature has the highest score values
in the corresponding group compared to all other groups.
}
\examples{
# TODO
# grouped_distributions(data, columns_obs = c('feature1', 'feature2'), ref_obs = 'group_column')


## TODO: we need to warn right away that the scores have to be in!
# grouped_distributions_heatmap?
}
