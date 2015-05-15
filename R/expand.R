#' Reweight a Survey to Marginal Controls
#' 
#' @param survey a \code{data.frame} with survey responses, including
#'   necessary colums for matching to marginals.
#'   
#' @param weight_var existing weights variable. If \code{NULL}, defaults to 1.
#' 
#' @param variables a character vector of the variables \code{survey} that will
#'   match against marginal tables.
#' @param marginals a list of \code{data.frame} objects, each a marginal table corresponding 
#'   to \code{variables}. The tables should have rows for each \code{geo_var} and
#'   columns corresponding to possible values for the appropriate \code{var*}
#'   variable.
#' @param verbose Print the maximum expansion factor with each iteration? 
#'   Default \code{FALSE}. 
#'   
#' @return a vector of weights for each row in \code{survey}
#' 
#' @export
#' 
#' @import dplyr
reweight_survey <- function(survey, weight_var = NULL, variables, marginals,
                            relative_gap = 0.01, max_iterations = 50, 
                            verbose = FALSE){
  
  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    survey <- dplyr::mutate(survey, weight = 1)
  } else {
    survey <- dplyr::rename_(survey, weight = weight_var)
  }
  
  
  # IPF ---
  for(iter in 1:max_iterations){
    
    # For each marginal table
    for (i in length(marginals)) {
      
      variable <- variables[i]
      
      marginal <- marginals[[i]]  %>%
        gather_(variable, "value", names(.)[-1], convert = TRUE) %>% 
        group_by_("geoid") %>%
        
        # Calculate mariginal contribution
        mutate(m_total = sum(value)) %>%
        select_(variable, "m_total")
      
      
      suppressMessages(
      survey <- survey %>%
        # Calculate survey contribution
        group_by_("geoid", as.name(variable)) %>%
        mutate(s_total = sum(weight)) %>%
      
        # join marginal table
        left_join(marginal) %>%
        
        # calculate new expansion factor
        mutate(
          exp_factor = m_total / s_total,
          weight = weight * exp_factor
        ) %>%
        select(-m_total)
      )
       
    }
    
    gap <- abs(max(survey$exp_factor) - 1)
    
    # print statistic
    if(verbose){print(gap)}
    
    # check tolerance
    if(gap < relative_gap){
      print(paste("Converged after", iter, "iterations"))
      break
    }
    
    
  }
      
  # if iterations exceeded, throw a warning.
  if(iter == max_iterations){
    warning("Failed to converge after ", iter, " iterations")
  }
  
  # return the new weights
  survey$weight
}

   







