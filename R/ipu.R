#' Iterative Proportional Updating
#' 
#' @description A special case of iterative proportional fitting. It can satisfy
#'   two, disparate sets of marginals that do not agree on a single total. A
#'   common example is balancing population data using household- and person-level
#'   marginal controls. This could be for survey expansion or synthetic
#'   population creation. The second set of marginal/seed data is optional, meaning
#'   it can also be used for more basic IPF tasks.
#'   
#'   Vignette: \url{http://pbsag.github.io/ipfr/}
#' 
#' @references \url{http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf}
#' 
#' @param primary_seed In population synthesis or household survey expansion, 
#'   this would be the household seed table (each record would represent a 
#'   household). It could also be a trip table, where each row represents an 
#'   origin-destination pair. Must contain a \code{pid} ("primary ID") field
#'   that is unique for each row. Must also contain a geography field that
#'   starts with "geo_".
#' 
#' @param primary_targets A \code{named list} of data frames.  Each name in the 
#'   list defines a marginal dimension and must match a column from the 
#'   \code{primary_seed} table. The data frame associated with each named list
#'   element must contain a geography field (starts with "geo_"). Each row in
#'   the target table defines a new geography (these could be TAZs, tracts,
#'   clusters, etc.). The other column names define the marginal categories that
#'   targets are provided for. The vignette provides more detail.
#' 
#' @param secondary_seed Most commonly, if the primary_seed describes households, the 
#'   secondary seed table would describe a unique person with each row. Must
#'   also contain the \code{pid} column that links each person to their 
#'   respective household in \code{primary_seed}. Must not contain any geography
#'   fields (starting with "geo_").
#' 
#' @param secondary_targets Same format as \code{primary_targets}, but they constrain 
#'   the \code{secondary_seed} table.
#' 
#' @param relative_gap After each iteration, the weights are compared to the
#' previous weights and the %RMSE is calculated. If the %RMSE is less than
#' the \code{relative_gap} threshold, then the process terminates.
#' 
#' @param max_iterations maximimum number of iterations to perform, even if 
#'    \code{relative_gap} is not reached.
#'    
#' @param absolute_diff Upon completion, the \code{ipu()} function will report
#'   the worst-performing marginal category and geography based on the percent
#'   difference from the target. \code{absolute_diff} is a threshold below which
#'   percent differences don't matter.
#'   
#'   For example, if if a target value was 2, and the expanded weights equaled
#'   1, that's a 100% difference, but is not important because the absolute value
#'   is only 1.
#'   
#'   Defaults to 10.
#'   
#' @param min_weight Minimum weight to allow in any cell to prevent zero weights.
#'    Set to .0001 by default.  Should be arbitrarily small compared to your 
#'    seed table weights.
#'   
#' @param verbose Print details on the maximum expansion factor with each 
#'    iteration? Default \code{FALSE}. 
#' 
#' @return the \code{primary_seed} with a revised weight column.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' hh_seed <- data.frame(
#'   hhid = c(1, 2, 3, 4),
#'   siz = c(1, 2, 2, 1),
#'   weight = c(1, 1, 1, 1),
#'   geo_cluster = c(1, 1, 2, 2)
#' )
#' 
#' hh_targets <- list()
#' hh_targets$siz <- data.frame(
#'   geo_cluster = c(1, 2),
#'   `1` = c(75, 100),
#'   `2` = c(25, 150)
#' )
#' 
#' result <- ipu(hh_seed, hh_targets, max_iterations = 10)
#' }
#' 
#' @importFrom magrittr "%>%"
ipu <- function(primary_seed, primary_targets, secondary_seed = NULL, secondary_targets = NULL,
                relative_gap = 0.01, max_iterations = 100, absolute_diff = 10,
                min_weight = .0001, verbose = FALSE){
  
  # If person data is provided, both seed and targets must be
  if (xor(!is.null(secondary_seed), !is.null(secondary_targets))) {
    stop("You provided either secondary_seed or secondary_targets, but not both.")
  }
  
  # Check hh and person tables
  if (!is.null(secondary_seed)) {
    check_tables(primary_seed, primary_targets, secondary_seed, secondary_targets)
  } else {
    check_tables(primary_seed, primary_targets)
  }
  
  # Pull off the geo information into a separate equivalency table
  # to be used as needed.
  geo_equiv <- primary_seed %>%
    dplyr::select(dplyr::starts_with("geo_"), "pid")
  primary_seed_mod <- primary_seed %>%
    dplyr::select(-dplyr::starts_with("geo_"))
  
  # Modify the household seed to the required format. Use one-hot-encoding to
  # expand standard columns into dummy columns.
  col_names <- names(primary_targets)
  primary_seed_mod <- primary_seed_mod %>%
    # Keep only the fields of interest (marginal columns and pid)
    dplyr::select(dplyr::one_of(c(col_names, "pid"))) %>%
    # Convert to factors and then to dummy columns
    dplyr::mutate_at(
      .vars = col_names,
      .funs = dplyr::funs(as.factor(.))
    ) %>%
    mlr::createDummyFeatures()
  
  if (!is.null(secondary_seed)) {
    # Modify the person seed table the same way, but sum by household ID
    col_names <- names(secondary_targets)
    secondary_seed_mod <- secondary_seed %>%
      # Keep only the fields of interest
      dplyr::select(dplyr::one_of(c(col_names, "pid"))) %>%
      dplyr::mutate_at(
        .vars = col_names,
        .funs = dplyr::funs(as.factor(.))
      ) %>%
      mlr::createDummyFeatures() %>%
      dplyr::group_by(pid) %>%
      dplyr::summarize_all(
        .funs = sum
      )
    
    # combine the hh and per seed tables into a single table
    seed <- primary_seed_mod %>%
      dplyr::left_join(secondary_seed_mod, by = "pid")
  } else {
    seed <- primary_seed_mod
  }
  
  # Add the geo information back.
  seed <- seed %>%
    dplyr::mutate(weight = 1)  %>%
    dplyr::left_join(geo_equiv, by = "pid")
  
  # store a vector of attribute column names to loop over later.
  # don't include 'pid' or 'weight' in the vector.
  geo_pos <- grep("geo_", colnames(seed))
  pid_pos <- grep("pid", colnames(seed))
  weight_pos <- grep("weight", colnames(seed))
  seed_attribute_cols <- colnames(seed)[-c(geo_pos, pid_pos, weight_pos)]
  
  # modify the targets to match the new seed column names and
  # join them to the seed table
  if (!is.null(secondary_seed)) {
    targets <- c(primary_targets, secondary_targets)
  } else {
    targets <- primary_targets
  }
  for (name in names(targets)) {
    # targets[[name]] <- targets[[name]] %>%
    temp <- targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -dplyr::starts_with("geo_")) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target")) %>%
      tidyr::spread(key = key, value = target)
    
    # Get the name of the geo column
    pos <- grep("geo_", colnames(temp))
    geo_colname <- colnames(temp)[pos]
    
    seed <- seed %>%
      dplyr::left_join(temp, by = geo_colname)
  }
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # Loop over each target and upate weights
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      
      # Get the name of the geo column
      target_tbl <- targets[[target_tbl_name]]
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      # Adjust weights
      seed <- seed %>%
        dplyr::mutate(
          geo = !!as.name(geo_colname),
          attr = !!as.name(seed_attribute),
          target = !!as.name(target_name)
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(
          factor = target / sum(attr * weight),
          weight = ifelse(attr > 0, weight * factor, weight),
          # Implement the floor on minimum weight
          weight = pmax(weight, min_weight)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-geo, -attr, -target, -factor)
    }
    
    # Determine percent differences (by geo field)
    pct_diff <- 0
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      target_tbl <- targets[[target_tbl_name]]
      
      # Get the name of the geo column
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      diff_tbl <- seed %>%
        dplyr::filter((!!as.name(seed_attribute)) > 0) %>%
        dplyr::select(
          geo = !!geo_colname, pid, attr = !!seed_attribute, weight,
          target = !!target_name
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(
          abs_diff = abs((target - sum(attr * weight))),
          pct_diff = abs_diff / (target + .0000001) # avoid dividing by zero
        ) %>%
        # Removes rows where the absolute gap is smaller than 'absolute_diff'
        dplyr::filter(abs_diff > absolute_diff) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      
      # If any records are left in the diff_tbl, record worst percent difference 
      # and save that percent difference table for reporting.
      if (nrow(diff_tbl) > 0) {
        if (max(diff_tbl$pct_diff) > pct_diff) {
          pct_diff <- max(diff_tbl$pct_diff)
          saved_diff_tbl <- diff_tbl
          saved_category <- seed_attribute
          saved_geo <- geo_colname
        }
      }
      
    }
    
    # Test for convergence
    if (iter > 1) {
      rmse <- mlr::measureRMSE(prev_weights, seed$weight)
      pct_rmse <- rmse / mean(prev_weights) * 100
      converged <- ifelse(pct_rmse <= relative_gap, TRUE, FALSE)
      if(verbose){
        cat("\r Finished iteration ", iter, ". %RMSE = ", pct_rmse)
      }
    }
    prev_weights <- seed$weight
    iter <- iter + 1
  }
  
  if (verbose) {
    message(ifelse(converged, "IPU converged", "IPU did not converge"))
    message("Worst marginal stats:")
    position <- which(saved_diff_tbl$pct_diff == pct_diff)[1]
    message("Category: ", saved_category)
    message(saved_geo, ": ", saved_diff_tbl$geo[position])
    message("Max % Diff: ", round(pct_diff * 100, 2), "%")
    message("Absolute Diff: ", round(saved_diff_tbl$abs_diff[position], 2))
    utils::flush.console()
  }
  
  primary_seed$weight <- seed$weight
  return(primary_seed)
  
}

#' Check seed and target tables for completeness
#' 
#' @description Given seed and targets, checks to make sure that at least one
#'   observation of each marginal category exists in the seed table.  Otherwise,
#'   ipf/ipu would produce wrong answers without throwing errors.
#'
#' @inheritParams ipu

check_tables <- function(primary_seed, primary_targets, secondary_seed = NULL, secondary_targets = NULL){
  
  # If person data is provided, both seed and targets must be
  if (xor(!is.null(secondary_seed), !is.null(secondary_targets))) {
    stop("You provided either secondary_seed or secondary_targets, but not both.")
  }
  
  ## Primary checks ##
  
  # Check that there are no NA values in seed or targets
  if (any(is.na(unlist(primary_seed)))) {
    stop("primary_seed table contains NAs")
  }
  if (any(is.na(unlist(primary_targets)))) {
    stop("primary_targets table contains NAs")
  }
  
  # Check that primary_seed table has a pid field and that it has a unique
  # value on each row.
  if (!"pid" %in% colnames(primary_seed)) {
    stop("The primary seed table does not have field 'pid'.")
  }
  unique_pids <- unique(primary_seed$pid)
  if (length(unique_pids) != nrow(primary_seed)) {
    stop("The primary seed's pid field has duplicate values.")
  }
  
  # check hh tables for correctness
  for (name in names(primary_targets)) {
    tbl <- primary_targets[[name]]
    
    # Check that each target table has a geo field
    check <- grepl("geo_", colnames(tbl))
    if (!any(check)) {
      stop("primary_target table '", name, "' does not have a geo column (must start with 'geo_')")
    }
    if (sum(check) > 1) {
      stop("primary_target table '", name, "' has more than one geo column (starts with 'geo_'")
    }
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(tbl))
    geo_colname <- colnames(tbl)[pos]
    
    # Get vector of other column names
    col_names <- colnames(tbl)
    col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
    
    # Check that at least one observation of the current target is in every geo
    for (geo in unique(unlist(primary_seed[, geo_colname]))){
      test <- match(col_names, primary_seed[[name]][primary_seed[, geo_colname] == geo])
      if (any(is.na(test))) {
        prob_cat <- col_names[which(is.na(test))]
        stop(
          "Marginal ", name, "; category ", prob_cat[1], " is missing from ",
          geo_colname, " ", geo, " in the primary_seed table."
        )
      }   
    }
  }
  
  
  ## Secondary checks (if provided) ##
  
  if (!is.null(secondary_seed)) {
    # Check for NAs
    if (any(is.na(unlist(secondary_seed)))) {
      stop("secondary_seed table contains NAs")
    }
    if (any(is.na(unlist(secondary_targets)))) {
      stop("secondary_targets table contains NAs")
    }
    
    # Check that secondary seed table has a pid field
    if (!"pid" %in% colnames(secondary_seed)) {
      stop("The primary seed table does not have field 'pid'.")
    }
    
    # Check that the secondary seed table does not have any geo columns
    check <- grepl("geo_", colnames(secondary_seed))
    if (any(check)) {
      stop("Do not include geo fields in the secondary_seed table (primary_seed only).")
    }
    
    # check the per tables for correctness
    for (name in names(secondary_targets)) {
      tbl <- secondary_targets[[name]]
      
      # Check that each target table has a geo field
      check <- grepl("geo_", colnames(tbl))
      if (!any(check)) {
        stop("secondary_target table '", name, "' does not have a geo column (must start with 'geo_')")
      }
      if (sum(check) > 1) {
        stop("secondary_target table '", name, "' has more than one geo column (starts with 'geo_'")
      }
      
      # Get the name of the geo field
      pos <- grep("geo_", colnames(tbl))
      geo_colname <- colnames(tbl)[pos]
      
      # Add the geo field from the primary_seed before checking
      secondary_seed <- secondary_seed %>%
        dplyr::left_join(
          primary_seed %>% dplyr::select(pid, geo_colname),
          by = "pid"
        )
      
      # Get vector of other column names
      col_names <- colnames(tbl)
      col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
      
      # Check that at least one observation of the current target is in every geo
      for (geo in unique(unlist(secondary_seed[, geo_colname]))){  
        test <- match(col_names, secondary_seed[[name]][secondary_seed[, geo_colname] == geo])
        if (any(is.na(test))) {
          prob_cat <- col_names[which(is.na(test))]
          stop(
            "Marginal ", name, "; category ", prob_cat[1], " is missing from ",
            geo_colname, " ", geo, " in the secondary_seed table."
          )
        }   
      }
    }
  }
}


