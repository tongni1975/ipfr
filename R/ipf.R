#' Reweight a Seed Table to Marginal Controls
#' 
#' @param seed a \code{data.frame} including necessary colums for matching to 
#' marginals and an optional weight field.
#'   
#' @param weight_var existing weights column in seed table. 
#' If \code{NULL}, defaults to 1.
#' 
#' @param marginals a two-column \code{data.frame}.  Column names must be 
#'    \code{marginal} and \code{value}.  The \code{marginal} column is filled  
#'    with combinations of marginal names and category numbers (e.g. "persons1",
#'    or "wrk2").
#'    
#'    The character portion (\code{persons} or \code{wrk}) must 
#'    correspond to a column name in the \code{seed} table.
#'    
#'    The numeric portion must correspond to the values found in the column of
#'    the \code{seed} table.
#'
#' @param relative_gap defines convergence.  If if no cells are factored by more 
#'    than this amount, the process is said to have converged.  For example, the
#'    default value of .01 imlies that convergence is reached when no cell in 
#'    seed table is changed by more than 10 percent.
#'
#' @param max_iterations maximimum number of iterations to perform, even if 
#'    convergence is not reached.
#'
#' @param min_weight Minimum weight to allow in any cell to prevent zero weights.
#'    Set to .0001 by default.  Should be arbitrarily small compared to your 
#'    seed table weights.
#'   
#' @param verbose Print the maximum expansion factor with each iteration? 
#'   Default \code{FALSE}. 
#'   
#' @return a vector of weights for each row in \code{seed}
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
#' 
ipf <- function(seed, weight_var = NULL, marginals, relative_gap = 0.01,
                max_iterations = 50, min_weight = .0001, verbose = FALSE){

  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    warning("weight_var not specified.  Initializing with equal weights.")
    seed <- dplyr::mutate(seed, weight = 1)
  } else {
    seed <- dplyr::rename_(seed, weight = weight_var)
  }
  
  # Split the marginal column into marginal and category columns
  marginals <- marginals %>%
    dplyr::mutate(
      category = as.numeric(gsub("[A-z]", "", marginal)),
      marginal = gsub("[0-9]", "", marginal)
    )
  
  # Check to see if the marginal totals match
  equal <- marginals %>%
    dplyr::group_by(marginal) %>%
    summarize(total = sum(value)) %>%
    .$total
  equal <- all(max(equal) - min(equal) == 0)
  if (!equal){
    warning(paste0(
      "Marginal totals are not equivalent. ",
      "The percentage distribution will still match all marginals. ",
      "Final weight total will match first marginal."
    ))
  }
  
  # Create a percent column to use in the IPF
  marginals <- marginals %>%
    dplyr::group_by(marginal) %>%
    dplyr::mutate(m_pct = value / sum(value))
  
  # IPF ---
  
  iter <- 1
  while (!converged & iter <= max_iterations){
    
    # In the following loop, track the maximum gap in this vector
    gap <- vector("numeric", nrow(marginals))
    
    # For each row in the marginal table
    for (i in 1:nrow(marginals)) {
      mName <- marginals[[1, i]]
      
      supressMessages(
        seed_summary <- seed %>%
          dplyr::group_by_(mName) %>%
          dplyr::summarize(totalweight = sum(weight)) %>%
          dplyr::mutate(s_pct = totalweight / sum(totalweight)) %>%
          dplyr::left_join(marginals, setNames("category", mName)) %>%
          dplyr::filter(marginal == mName) %>%
          dplyr::mutate(factor = m_pct / s_pct) %>%
          dplyr::select_(mName, "factor")
      )
      
      suppressMessages(
        seed <- seed %>%
          dplyr::left_join(seed_summary) %>%
          # calc new weight while preventing the creation of zeros
          dplyr::mutate(
            weight = weight * factor,
            weight = ifelse(weight < min_weight, min_weight, weight)
          )
      )
      
      gap[i] <- max(abs(seed$factor - 1))
      
      seed <- seed %>% dplyr::select(-factor)
    }
    
    # Check for convergence and increment iter
    converged <- all(gap <= min_weight)
    iter = iter + 1
  }
  
  # When finished, scale up the weights to match the first marginal total
  seed$weight <- seed$weight * (sum(marginals[[1]]$value) / sum(seed$weight))
        
  # if iterations exceeded, throw a warning.
  if(iter == max_iterations){
    warning("Failed to converge after ", iter, " iterations")
  }
  
  # return the new weights
  seed$weight
}

   







