#' Reweight a Seed Table to Marginal Controls
#' 
#' @param seed a \code{data.frame} including necessary colums for matching to 
#' marginals and an optional weight field.
#'   
#' @param weight_var existing weights column in seed table. 
#' If \code{NULL}, defaults to 1.
#' 
#' @param marginals a named list of \code{data.frame} objects. Each name must 
#'    match a column in the seed table, with each table describing the marginal
#'    distribution. The tables should have rows corresponding to 
#'    possible values for the appropriate \code{var*} variable.
#'
#' @param relative_gap defines convergence.  If if no cells are factored by more 
#'    than this amount, the process is said to have converged.  For example, the
#'    default value of .01 imlies that convergence is reached when no cell in 
#'    seed table is changed by more than 10 percent.
#'
#' @param max_iterations maximimum number of iterations to perform, even if 
#'    convergence is not reached.
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
                max_iterations = 50, verbose = FALSE){
  
  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    seed <- dplyr::mutate(seed, weight = 1)
  } else {
    seed <- dplyr::rename_(seed, weight = weight_var)
  }
  
  # Check to see if the marginal totals match
  if (length(marginals) > 1){
    for (i in 1:length(marginals)){
      marginal <- marginals[[i]]
      
      if (i == 1){
        check <- c(sum(marginal$value))
      } else {
        check <- append(check, sum(marginal$value))
      }
    }
    if (!all(max(check) - min(check) == 0)){
      warning(paste0(
        "Marginal totals are not equivalent. ",
        "The percentage distribution will still match all marginals. ",
        "Final weight total will match first marginal."
      ))
    }
  }
  
  variables <- names(marginals)
  
  # IPF ---
  gap <- vector("numeric", length(marginals))
  
  for(iter in 1:max_iterations){
    
    # For each marginal table
    for (i in 1:length(marginals)) {
      
      variable <- variables[i]
      
      marginal <- marginals[[variable]]  %>%
        #tidyr::gather_(variable, "value", colnames(.), convert = TRUE) %>%
        
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
            weight = ifelse(weight < .0001, .0001, weight)
            )
      )
      
      # convergence criteria for each marginal
      gap[i] <- abs(max(seed$factor) - 1)
      
      seed <- seed %>% select(-factor)
    }
    
    # print statistic
    if(verbose){print(gap)}
    
    # check tolerance
    if(max(gap) < relative_gap){
      print(paste("Converged after", iter, "iterations"))
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

   







