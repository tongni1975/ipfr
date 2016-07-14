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
#'    Default \code{FALSE}. 
#'   
#' @return a vector of weights for each row in \code{seed}
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
#' 
ipf <- function(seed, weight_var = NULL, marginals,
                relative_gap = 0.01, max_iterations = 50, min_weight = .0001,
                verbose = FALSE){

  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    warning("weight_var not specified.  Initializing with equal weights.")
    utils::flush.console()
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
  vals <- marginals %>%
    dplyr::group_by(marginal) %>%
    dplyr::summarize(total = sum(value)) %>%
    .$total
  equal <- all(max(vals) - min(vals) == 0)
  if (!equal){
    warning(paste0(
      "Marginal totals are not equivalent. ",
      "The percentage distribution will still match all marginals. ",
      "Final weight total will match first marginal."
    ))
    utils::flush.console()
  }
  
  # If all the marginals add up to zero, return vector of zeros.
  if (sum(vals) == 0) {
    seed$weight <- 0
    return(seed$weight)
  }
  
  # Create a percent column to use in the IPF
  marginals <- marginals %>%
    dplyr::group_by(marginal) %>%
    dplyr::mutate(m_pct = value / sum(value))
  
  # Check that at least one row of seed information exists for each
  # marginal category.
  margs <- unique(marginals$marginal)
  temp <- list()
  for (marg in margs){
    cats <- marginals %>%
      dplyr::ungroup() %>%
      dplyr::filter(marginal == marg)
    cats <- cats$category
    
    for (cat in cats){
      ok <- any(seed[, marg] == cat)
      
      if (!ok) {
        warning(paste0(
          "The seed table has no observations of marginal: ", marg,
          " category: ", cat
        ))
        utils::flush.console()
        return()
      }
    }
  }
  
  # IPF ---
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations){
    
    # In the following loop, track the maximum gap in this vector
    gap <- vector("numeric", nrow(marginals))
    
    # For each row in the marginal table
    for (i in 1:nrow(marginals)) {
      mName <- marginals[[i, 1]]
      
      suppressMessages(
        seed_summary <- seed %>%
          dplyr::group_by_(mName) %>%
          dplyr::summarize(totalweight = sum(weight)) %>%
          dplyr::mutate(s_pct = totalweight / sum(totalweight)) %>%
          dplyr::left_join(marginals, stats::setNames("category", mName)) %>%
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
  firstMarg <- marginals$marginal[1]
  target <- marginals %>%
    dplyr::summarize(total = sum(value))
  target <- target$total[1]
  
  seed$weight <- seed$weight * (target) / sum(seed$weight)
        
  # if iterations exceeded, throw a warning.
  if(iter == max_iterations & !converged){
    warning("Failed to converge after ", iter, " iterations")
    utils::flush.console()
  }
  
  # return the new weights
  seed$weight
}

   


#' Perform the ipf procedure for multiple marginal sets and returns a 
#' \code{data.frame}.
#' 
#' @inheritParams ipf
#' 
#' @param margTbl A \code{data.frame} with a row for each set of marginals to
#'    use in an ipf procedure.  Every column other than the \code{id_field} will
#'    be treated as a marginal column.
#'    
#' @param id_field The name of the identifying field not to be used as a 
#'    marginal.
#'
#' @return a \code{data.frame} with a row for each \code{id_field}. The columns 
#'    of which contain the joint-distribution weights after fitting to marginals.
#'
#' @export
#'    
ipf_multi <- function(seed, weight_var = "weight", margTbl, id_field,
                   relative_gap = 0.01, max_iterations = 50, 
                   min_weight = .0001, verbose = FALSE){
  
  # set weights variable ----
  if(is.null(weight_var)){
    # if none given, set to 1.
    warning("weight_var not specified.  Initializing with equal weights.")
    utils::flush.console()
    seed <- dplyr::mutate(seed, weight = 1)
  } else {
    seed <- dplyr::rename_(seed, weight = weight_var)
  }
  
  # Collect the marginal names
  mNameCats <- margTbl %>%
    dplyr::select(-ID) %>%
    names()
  
  mNames <- mNameCats %>%
    gsub("[0-9]", "", .) %>%
    unique()
    
  # Expand the seed table to be repeated for each id_field, add marginal, and
  # convert weight field to a percent.
  seed_long <- merge(margTbl$ID, seed) %>%
    dplyr::rename(ID = x) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(weight = weight / sum(weight)) %>%
    dplyr::ungroup()
  
  for (i in 1:length(mNames)) {
    name <- mNames[i]
    namecats <- mNameCats[grep(name, mNameCats)]
    
    # save first marginal for later scaling
    if (i == 1) {
      firstMarg <- margTbl %>%
        dplyr::select(ID, dplyr::one_of(namecats)) %>%
        tidyr::gather(key = m, value = value, -ID) %>%
        dplyr::group_by(ID) %>%
        dplyr::summarize(total = sum(value))
    }
    
    tojoin <- margTbl %>%
      dplyr::select(ID, dplyr::one_of(namecats)) %>%
      tidyr::gather(key = m, value = value, -ID) %>%
      dplyr::group_by(ID) %>%
      dplyr::mutate(
        m = as.numeric(gsub("[A-z]", "", m)),
        total = sum(value),
        value = ifelse(total == 0, 0, value / total)
      ) %>%
      dplyr::rename_(
        .dots = stats::setNames(c("m", "value"), c(name, paste0(name,"marg")))
      ) %>%
      dplyr::ungroup()
    
    seed_long <- seed_long %>%
      dplyr::left_join(tojoin)
  }
  
  # Perform IPF
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # In the following loop, track the maximum gap in this vector
    gap <- vector("numeric", length(mNames))
    
    # For each marginal
    for (i in 1:length(mNames)) {
      name <- mNames[i]
      namecat <- paste0(name,"marg")
      
      dots <- list(
        lazyeval::interp(~x / y, x = as.name(namecat), y = as.name("total"))
      )
      
      seed_long <- seed_long %>%
        dplyr::group_by_("ID", name) %>%
        dplyr::mutate(total = sum(weight)) %>%
        dplyr::mutate_(.dots = stats::setNames(dots, "factor")) %>%
        dplyr::mutate(
          weight = weight * factor,
          # Don't allow weight to drop below min_weight
          weight = ifelse(weight < min_weight, min_weight, weight),
          # If a factor was 0 (meaning the marginal was zero), set the factor
          # to 1 to allow convergence.
          factor = ifelse(factor == 0, 1, factor)
        )
      
      gap[i] <- max(abs(seed_long$factor - 1))
    }
    
    # Check for convergence and increment iter
    converged <- all(gap <= min_weight)
    iter = iter + 1
  }
  
  # if iterations exceeded, throw a warning.
  if(iter > max_iterations){
    warning("Failed to converge after ", max_iterations, " iterations")
    utils::flush.console()
  }
  
  # When finished, scale up the weights to match the first marginal total
  seed_long <- seed_long %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, dplyr::one_of(mNames), weight) %>%
    dplyr::left_join(firstMarg) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(weight = weight * total / sum(weight)) %>%
    dplyr::select(-total) %>%
    dplyr::ungroup()
    
  # Format final data frame
  final <- seed_long
  for (name in mNames) {
    final[, name] <- paste0(name, final[[name]])
  }
  #final$category <- do.call(paste0, select(final, one_of(mNames)))
  final <- tidyr::unite(final, col = category, dplyr::one_of(mNames), sep = "_")
  final <- final %>%
    dplyr::select(ID, category, weight) %>%
    tidyr::spread(key = category, value = weight)
  
  return(final)
}





