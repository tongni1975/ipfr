#' Creates a synthetic population based on ipu results
#' 
#' A simple function that takes the \code{weight_tbl} output from
#' \code{\link{ipu}} and randomly samples based on the weight.
#' 
#' @inheritParams ipu
#' @param weight_tbl the \code{data.frame} of the same name output by
#'   \code{\link{ipu}}.
#' @param group_by if provided, the random sampling will be done within each
#'   group.
#' @return A \code{data.frame} with one record for each synthesized member of
#'   the population (e.g. household). A \code{new_id} column is created, but
#'   the previous \code{primary_id} column is maintained to facilitate joining
#'   back to other data sources (e.g. a person table).
#' @export

synthesize <- function(weight_tbl, group_by = NULL, primary_id = "id") {
  
  if (!primary_id %in% colnames(weight_tbl)) {
    stop("primary_id not found in weight_tbl") # nocov
  }
  
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(weight_tbl)) {
      stop("group_by not found in weight_tbl") # nocov
    }
    weight_tbl <- weight_tbl %>%
      group_by(!!as.name(group_by))
  }
  
  synthetic_table <- weight_tbl %>%
    sample_n(round(sum(weight), 0), replace = TRUE, weight = weight) %>%
    ungroup() %>%
    select(-weight, -avg_weight, -weight_factor) %>%
    mutate(new_id = seq(1, n())) %>%
    select(new_id, one_of(primary_id), everything())
  
  return(synthetic_table)
}
