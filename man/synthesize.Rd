% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synthesize.R
\name{synthesize}
\alias{synthesize}
\title{Creates a synthetic population based on ipu results}
\usage{
synthesize(weight_tbl, group_by = NULL, primary_id = "id")
}
\arguments{
\item{weight_tbl}{the \code{data.frame} of the same name output by
\code{\link{ipu}}.}

\item{group_by}{if provided, the \code{data.frame} will be grouped by this
variable before sampling. If not provided, tidyverse/dplyr groupings will
be respected. If no grouping info is present, samples are drawn from the
entire table.}

\item{primary_id}{The field used to join the primary and secondary seed
tables. Only necessary if \code{secondary_seed} is provided.}
}
\value{
A \code{data.frame} with one record for each synthesized member of
  the population (e.g. household). A \code{new_id} column is created, but
  the previous \code{primary_id} column is maintained to facilitate joining
  back to other data sources (e.g. a person attribute table).
}
\description{
A simple function that takes the \code{weight_tbl} output from
\code{\link{ipu}} and randomly samples based on the weight.
}
\examples{
hh_seed <- dplyr::tibble(
id = c(1, 2, 3, 4),
siz = c(1, 2, 2, 1),
weight = c(1, 1, 1, 1),
geo_cluster = c(1, 1, 2, 2)
)
hh_targets <- list()
hh_targets$siz <- dplyr::tibble(
  geo_cluster = c(1, 2),
  `1` = c(75, 100),
  `2` = c(25, 150)
)
result <- ipu(hh_seed, hh_targets, max_iterations = 5)
synthesize(result$weight_tbl, "geo_cluster")
}
