#' ipfr: A package to perform iterative proportional fitting
#' 
#' The two functions are \code{\link{ipf}} and
#'    \code{\link{ipf_multi}}
#' 
#' @docType package
#' @name ipfr
NULL
#> NULL

#' Reweight a Seed Table to Marginal Controls
#' 
#' @param targets A \code{named list} of data frames.  Each name in the list defines a
#'    marginal dimension and must match a column from the seed table.  The data
#'    frame associated with each name must start with an identical \code{ID} column.
#'    The other column names define the marginal categories that targets are
#'    provided for.  Each row of the data frames represents a unique set of
#'    targets to fit.  Every data frame must have the same number of rows.
#'    
#' @param seed A \code{data frame} including a \code{weight} field and necessary
#'    colums for matching to marginal targets.
#'
#' @param relative_gap target for convergence.  Maximum percent change to allow
#'    any seed weight to move by while considering the process converged.  By 
#'    default, if no weights change by more than 1%, the process has converged.
#'    The process is said to be converged if either \code{relative_gap} or 
#'    \code{absolute_gap} parameters have been met.
#'
#' @param absolute_gap target for convergence.  Maximum absolute change to allow
#'    any seed weight to move by while considering the process converged.  By 
#'    default, if no weights change by more than 10, the process has converged.
#'    The process is said to be converged if either \code{relative_gap} or 
#'    \code{absolute_gap} parameters have been met.
#'
#' @param max_iterations maximimum number of iterations to perform, even if 
#'    convergence is not reached.
#'
#' @param min_weight Minimum weight to allow in any cell to prevent zero weights.
#'    Set to .0001 by default.  Should be arbitrarily small compared to your 
#'    seed table weights.
#'   
#' @param verbose Print details on the maximum expansion factor with each 
#'    iteration? Default \code{FALSE}. 
#'   
#' @return the seed \code{data frame} with a column of weights appended for each
#'    row in the target dataframes
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
#' 
ipf <- function(seed, targets,
                relative_gap = 0.01, absolute_gap = 10, max_iterations = 50,
                min_weight = .0001, verbose = FALSE){

  # Check check that seed and target are provided
  if (is.null(seed)) {
    stop("Seed table not provided")
  }
  if (is.null(targets)) {
    stop("Targets not provided")
  }
  
  # Check that at least one observation of each marginal category exists
  # in the seed table.  Otherwise, the process produces wrong answers without
  # throwing errors.
  for (name in names(targets)) {
    col_names <- colnames(targets[[name]])
    col_names <- as.numeric(col_names[!col_names == "ID"])
    
    test <- match(col_names, seed[[name]])
    if (any(is.na(test))) {
      prob_cat <- col_names[which(is.na(test))]
      stop(paste0(
        "Marginal ", name, "; category ", prob_cat,
        " is missing from seed table"
      ))
    }
  }
  
  # Create df of totals from the first marginal table
  # (e.g. total households, persons, etc.)
  totals <- targets[[1]] %>%
    tidyr::gather(key = category, value = count, -ID) %>%
    group_by(ID) %>%
    summarize(total = sum(count))
  
  # Create a long data frame by repeating the seed table for each
  # row in the target tables.
  seed_long <- merge(totals$ID, seed) %>%
    dplyr::rename(ID = x) %>%
    dplyr::arrange(ID)
  
  # Convert the weights into percents.  The percents will sum to 1 for each ID
  seed_long <- seed_long %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(weight = weight / sum(weight)) %>%
    dplyr::ungroup()
  
  # IPF ---
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations){
    
    # In the following loop, track the maximum gap in this vector
    rel_gap <- vector("numeric", length(targets))
    abs_gap <- vector("numeric", length(targets))
    
    # For each table in targets
    for (i in 1:length(targets)) {
      mName <- names(targets)[i]
      
      # Prepare the target table
      target <- targets[[mName]] %>%
        tidyr::gather(key = marg, value = target, -ID) %>%
        dplyr::mutate(marg = as.numeric(marg)) %>%
        dplyr::group_by(ID) %>%
        dplyr::mutate(target = target / sum(target)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ID)
      colnames(target) <- c("ID", mName, "target")
      
      # Prepare a summary of the seed table
      seed_summary <- seed_long %>%
        dplyr::group_by_("ID", mName) %>%
        dplyr::summarize(weight = sum(weight)) %>%
        dplyr::ungroup()
      
      # Join target to seed to calculate factor
      fac_tbl <- seed_summary %>%
        dplyr::left_join(target, by = setNames(c("ID", mName), c("ID", mName))) %>%
        dplyr::mutate(factor = ifelse(weight == 0, 1, target / weight))
      
      # Handle missing targets.  For example, if the seed table had size
      # categories 1-5, but the targets only had 1-4, leave the 5s alone without
      # causing an error
      fac_tbl <- fac_tbl %>%
        dplyr::mutate(factor = ifelse(is.na(target), 1, factor ))
      
      # Prepare fac_tbl for joining to seed table
      fac_tbl <- fac_tbl %>%
        dplyr::select(dplyr::one_of("ID", mName, "factor"))
      
      # Join to the seed_long table and calculate new weight
      seed_long <- seed_long %>%
        dplyr::left_join(fac_tbl, by = setNames(c("ID", mName), c("ID", mName))) %>%
        dplyr::mutate(new_weight = weight * factor)
      
      # Because closure will depend on relative and absolute gap, add the
      # totals vector back to determine the absolute difference.
      gap_tbl <- seed_long %>%
        dplyr::left_join(totals, by = "ID") %>%
        dplyr::mutate(abs_diff = abs(factor - 1) * total)
      
      rel_gap[i] <- max(abs(gap_tbl$factor - 1))
      abs_gap[i] <- max(abs(gap_tbl$abs_diff))
      
      # Clean up seed_long for next iteration
      seed_long <- seed_long %>%
        dplyr::mutate(weight = new_weight) %>%
        dplyr::select(-c(factor, new_weight))
        
    }
    
    # Check for convergence and increment iter
    if(verbose){
      message("Finished iteration ", iter)
    }
    converged <- all(rel_gap <= relative_gap) | all(abs_gap <= absolute_gap)
    iter = iter + 1
  }
  
  # After the loop, scale up the weights to match the totals
  seed_long <- seed_long %>%
    dplyr::left_join(totals, by = "ID") %>%
    dplyr::mutate(weight = weight * total) %>%
    dplyr::select(-total)
  
  # Change seed_long into final format to return
  seed_long <- seed_long %>%
    tidyr::spread(key = ID, value = weight)
  
  # if iterations exceeded, throw a warning.
  if(iter == max_iterations & !converged){
    warning("Failed to converge after ", iter, " iterations")
    utils::flush.console()
  }
  
  if (verbose) {
    position <- which(rel_gap == max(rel_gap))[1]
    cat("\n", "Max Rel Gap:", rel_gap[position])
    cat("\n", "Absolute Gap:", abs_gap[position])
    cat("\n", "ID:", gap_tbl$ID[position])
    cat("\n", "Marg:", names(targets)[position])
    utils::flush.console()
  }
  
  # return the table with new weights
  return(seed_long)
}



