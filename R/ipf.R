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
#'    or "wrk2").  The character portion (\code{persons} or \code{wrk}) must 
#'    correspond to a column name in the \code{seed} table.
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
    seed <- dplyr::mutate(seed, weight = 1)
  } else {
    warning("weight_var not specified.  Initializing with equal weights.")
    seed <- dplyr::rename_(seed, weight = weight_var)
  }
  
  # Split the marginal column into marginal and category columns
  marginals <- marginals %>%
    dplyr::mutate(
      category = gsub("[A-z]", "", marginal),
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
  


  
  variables <- names(marginals)
  
  # IPF ---
  
  for(iter in 1:max_iterations){
    
    # For each marginal table
    for (i in 1:length(marginals)) {
      
      variable <- variables[i]
      
      marginal <- marginals[[variable]]  %>%
        
        # Normalize marginals to percents
        dplyr::mutate(m_pct = value / sum(value)) %>%
        dplyr::select_(variable, "m_pct")
      
      suppressMessages(
        seed_summary <- seed %>%
          dplyr::group_by_(variable) %>%
          dplyr::summarize(totalweight = sum(weight)) %>%
          dplyr::mutate(s_pct = totalweight / sum(totalweight)) %>%
          dplyr::left_join(marginal) %>%
          dplyr::mutate(factor = m_pct / s_pct) %>%
          dplyr::select_(variable, "factor")
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
      
      seed <- seed %>% dplyr::select(-factor)
    }
    
    # check tolerance after each pass
    converged <- FALSE
    for (i in 1:length(marginals)){
      variable <- variables[i]
      
      marginal <- marginals[[variable]]  %>%
        
        # Normalize marginals to percents
        dplyr::mutate(m_pct = value / sum(value)) %>%
        dplyr::select_(variable, "m_pct")
      
      suppressMessages(
        seed_summary <- seed %>%
          dplyr::group_by_(variable) %>%
          dplyr::summarize(totalweight = sum(weight)) %>%
          dplyr::mutate(s_pct = totalweight / sum(totalweight)) %>%
          dplyr::left_join(marginal) %>%
          dplyr::mutate(factor = m_pct / s_pct) %>%
          dplyr::select_(variable, "factor")
      )
      
      if (i == 1){
        max_factor <- max(seed_summary$factor)
      } else {
        max_factor <- max(max_factor, max(seed_summary$factor))
      }
      
      if (abs(max_factor - 1) < relative_gap){
        converged <-  TRUE
        print(paste0("Converged after ", iter, " iterations."))
        break
      }
    }
    
    if (converged){
      break
    }
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

   







