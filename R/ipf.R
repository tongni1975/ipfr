#' Reweight a Seed Table to Marginal Controls
#' 
#' @param seed a \code{data.frame} including necessary colums for matching to 
#' marginals and an optional weight field.
#'   
#' @param weight_var existing weights column in seed table. 
#' If \code{NULL}, defaults to 1.
#' 
#' @param variables a character vector of the table variable names that will
#'   match against marginal tables.
#' 
#' @param marginals a list of \code{data.frame} objects, each a marginal table corresponding 
#'   to \code{variables}. The tables should have columns corresponding to 
#'   possible values for the appropriate \code{var*} variable.
#'   
#' @param verbose Print the maximum expansion factor with each iteration? 
#'   Default \code{FALSE}. 
#'   
#' @return a vector of weights for each row in \code{seed}
#' 
#' @export
#' 
#' @import dplyr
ipf <- function(seed, weight_var = NULL, variables, marginals,
                            relative_gap = 0.01, max_iterations = 50, 
                            verbose = FALSE){
  
  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    seed <- dplyr::mutate(seed, weight = 1)
  } else {
    seed <- dplyr::rename_(seed, weight = weight_var)
  }
  
  # IPF ---
  gap <- vector("numeric", length(marginals))
  
  for(iter in 1:max_iterations){
    
    # For each marginal table
    for (i in 1:length(marginals)) {
      
      variable <- variables[i]
      
      marginal <- marginals[[i]]  %>%
        gather_(variable, "value", colnames(.), convert = TRUE) %>%
        
        # Normalize marginals to percents
        mutate(m_pct = value / sum(value)) %>%
        select_(variable, "m_pct")
      
      suppressMessages(
        seed_summary <- seed %>%
          group_by_(variable) %>%
          summarize(totalweight = sum(weight)) %>%
          mutate(s_pct = totalweight / sum(totalweight)) %>%
          left_join(marginal) %>%
          mutate(factor = m_pct / s_pct) %>%
          select(persons, factor)
      )
      
      suppressMessages(
        seed <- seed %>%
          left_join(seed_summary) %>%
          # calc new weight while preventing the creation of zeros
          mutate(
            weight = weight * factor,
            weight = ifelse(weight < .0001, .0001, weight)
            )
      )
      
      # convergence criteria for each marginal
      gap[i] <- abs(max(seed$factor) - 1)
    }
    
    # print statistic
    if(verbose){print(gap)}
    
    # check tolerance
    if(max(gap) < relative_gap){
      print(paste("Converged after", iter, "iterations"))
      break
    }
    
    
  }
      
  # if iterations exceeded, throw a warning.
  if(iter == max_iterations){
    warning("Failed to converge after ", iter, " iterations")
  }
  
  # return the new weights
  seed$weight
}

   







