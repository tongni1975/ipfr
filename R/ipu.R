#' Iterative Proportional Updating
#' 
#' @description A special case of iterative proportional fitting when two,
#' disparate sets of marginals need to be satisfied that do not agree on a
#' single total. A common example is balance population data using household-
#' and person-level marginal controls. This could be for survey expansion or
#' synthetic population creation. The original Arizona IPU method is improved
#' by dampening the adjustment factor.
#' 
#' @references \url{http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf}
#' 
#' @param hh_seed Same as \code{seed} in \link[ipfr]{ipf}. Each row is a
#' household. Must contain an \code{hhid} column.
#' 
#' @param hh_targets Same as \code{targets} in \link[ipfr]{ipf}.
#' 
#' @param per_seed Seed table of person attributes. Each row is a person.
#' Must also contain an \code{hhid} column that links each person to their
#' respective household in \code{hh_seed}.
#' 
#' @param per_targets Similar to \code{targets} from \link[ipfr]{ipf}, but for
#' person-level targets.
#' 
#' @param damp_factor Should be a number between 0 and 1. Reduces the adjustment
#'   factor calculated in each iteration before updating weights. Lower values
#'   require more iterations but are less likely to overfit a single marginal
#'   distribution at the expense of the others. Defaults to 0.75.
#' 
#' @param relative_gap After each iteration, the weights are compared to the
#' previous weights and an RMSE metric is calculated. If the RMSE is less than
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
#' 
#' @return the \code{hh_seed} with a revised weight column.
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
ipu <- function(hh_seed, hh_targets, per_seed, per_targets, damp_factor = .75,
                relative_gap = 0.01,max_iterations = 100, absolute_diff = 10,
                min_weight = .0001, verbose = FALSE){
  
  # Check hh and person tables
  check_tables(hh_seed, hh_targets, per_seed, per_targets)
  
  # Pull off the geo information into a separate equivalency table
  # to be used as needed.
  geo_equiv <- hh_seed %>%
    dplyr::select(dplyr::starts_with("geo_"), "hhid")
  hh_seed_mod <- hh_seed %>%
    dplyr::select(-dplyr::starts_with("geo_"))
  
  # Modify the household seed to the required format. Use one-hot-encoding to
  # expand standard columns into dummy columns.
  col_names <- names(hh_targets)
  hh_seed_mod <- hh_seed_mod %>%
    # Keep only the fields of interest (marginal columns and hhid)
    dplyr::select(dplyr::one_of(c(col_names, "hhid"))) %>%
    # Convert to factors and then to dummy columns
    dplyr::mutate_at(
      .vars = col_names,
      .funs = dplyr::funs(as.factor(.))
    ) %>%
    mlr::createDummyFeatures()
  
  # Modify the person seed table the same way, but sum by household ID
  col_names <- names(per_targets)
  per_seed_mod <- per_seed %>%
    # Keep only the fields of interest
    dplyr::select(dplyr::one_of(c(col_names, "hhid"))) %>%
    dplyr::mutate_at(
      .vars = col_names,
      .funs = dplyr::funs(as.factor(.))
    ) %>%
    mlr::createDummyFeatures() %>%
    dplyr::group_by(hhid) %>%
    dplyr::summarize_all(
      .funs = sum
    )
  
  # combine the hh and per seed tables into a single table
  # and add the geo information back.
  seed <- hh_seed_mod %>%
    dplyr::left_join(per_seed_mod, by = "hhid") %>%
    dplyr::mutate(weight = 1)  %>%
    dplyr::left_join(geo_equiv, by = "hhid")
  
  # store a vector of attribute column names to loop over later.
  # don't include 'hhid' or 'weight' in the vector.
  geo_pos <- grep("geo_", colnames(seed))
  hhid_pos <- grep("hhid", colnames(seed))
  weight_pos <- grep("weight", colnames(seed))
  seed_attribute_cols <- colnames(seed)[-c(geo_pos, hhid_pos, weight_pos)]
  
  # modify the targets to match the new seed column names and
  # join them to the seed table
  targets <- c(hh_targets, per_targets)
  for (name in names(targets)) {
    # targets[[name]] <- targets[[name]] %>%
    temp <- targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -dplyr::starts_with("geo_")) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target")) %>%
      spread(key = key, value = target)
    
    # Get the name of the geo column
    pos <- grep("geo_", colnames(temp))
    geo_colname <- colnames(temp)[pos]
    
    seed <- seed %>%
      left_join(temp, by = geo_colname)
  }
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # Loop over each target and upate weights
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      target_tbl <- targets[[target_tbl_name]]
      
      # Get the name of the geo column
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      # Calculate adjustment factor
      fac_tbl <- seed %>%
        dplyr::filter((!!as.name(seed_attribute)) > 0) %>%
        dplyr::select(
          geo = !!geo_colname, hhid, attr = !!seed_attribute, weight, 
          target = !!target_name
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(factor = (target) / sum(attr * weight) * damp_factor) %>%
        dplyr::ungroup() %>%
        dplyr::select(hhid, factor)
      
      seed <- seed %>%
        dplyr::left_join(fac_tbl, by = "hhid") %>%
        dplyr::mutate(weight = ifelse(!is.na(factor), weight * factor, weight)) %>%
        # Implement the floor on minimum weight
        dplyr::mutate(weight = pmax(weight, min_weight)) %>%
        dplyr::select(-factor)
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
          geo = !!geo_colname, hhid, attr = !!seed_attribute, weight,
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
      converged <- ifelse(rmse <= relative_gap, TRUE, FALSE)
    }
    prev_weights <- seed$weight
    iter <- iter + 1
  }
  
  if (verbose) {
    message("Relative Gap: ", rmse)
    message(ifelse(converged, "IPU converged", "IPU did not converge"))
    message("Worst marginal stats:")
    position <- which(saved_diff_tbl$pct_diff == pct_diff)[1]
    message("Category: ", saved_category)
    message(saved_geo, ": ", saved_diff_tbl$geo[position])
    message("Max % Diff: ", round(pct_diff * 100, 2), "%")
    message("Absolute Diff: ", round(saved_diff_tbl$abs_diff[position], 2))
    utils::flush.console()
  }
  
  hh_seed$weight <- seed$weight
  return(hh_seed)
  
}

#' Check seed and target tables for completeness
#' 
#' @description Given seed and targets, checks to make sure that at least one
#'   observation of each marginal category exists in the seed table.  Otherwise,
#'   ipf/ipu would produce wrong answers without throwing errors.
#'
#' @inheritParams ipu

check_tables <- function(hh_seed, hh_targets, per_seed, per_targets){
  
  # Check check that seed and target are provided
  if (is.null(hh_seed)) {
    stop("hh_seed table not provided")
  }
  if (is.null(hh_targets)) {
    stop("hh_targets not provided")
  }
  if (is.null(per_seed)) {
    stop("per_seed table not provided")
  }
  if (is.null(per_targets)) {
    stop("per_targets not provided")
  }
  
  # Check that there are no NA values in seed or targets
  if (any(is.na(unlist(hh_seed)))) {
    stop("hh_seed table contains NAs")
  }
  if (any(is.na(unlist(hh_targets)))) {
    stop("hh_targets table contains NAs")
  }
  if (any(is.na(unlist(per_seed)))) {
    stop("per_seed table contains NAs")
  }
  if (any(is.na(unlist(per_targets)))) {
    stop("per_targets table contains NAs")
  }
  
  # Check that the person seed table does not have any geo columns
  check <- grepl("geo_", colnames(per_seed))
  if (any(check)) {
    stop("Do not include geo fields in the per_seed table (hh_seed only).")
  }
  
  # check hh tables for completeness
  for (name in names(hh_targets)) {
    tbl <- hh_targets[[name]]
    
    # Check that each target table has a geo field
    check <- grepl("geo_", colnames(tbl))
    if (!any(check)) {
      stop("hh_target table '", name, "' does not have a geo column (must start with 'geo_')")
    }
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(tbl))
    geo_colname <- colnames(tbl)[pos]
    
    # Get vector of other column names
    col_names <- colnames(tbl)
    col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
    
    for (geo in hh_seed[, geo_colname]){
      test <- match(col_names, hh_seed[[name]][hh_seed[, geo_colname] == geo])
      if (any(is.na(test))) {
        prob_cat <- col_names[which(is.na(test))]
        stop(
          "Marginal ", name, "; category ", prob_cat[1], " is missing from ",
          geo_colname, " ", geo, " in the hh_seed table."
        )
      }   
    }
  }
  
  # check the per tables for completeness
  for (name in names(per_targets)) {
    tbl <- per_targets[[name]]
    
    # Check that each target table has a geo field
    check <- grepl("geo_", colnames(tbl))
    if (!any(check)) {
      stop("hh_target table '", name, "' does not have a geo column (must start with 'geo_')")
    }
    
    # Get the name of the geo field
    pos <- grep("geo_", colnames(tbl))
    geo_colname <- colnames(tbl)[pos]
    
    # Add the geo field from the hh_seed before checking
    per_seed <- per_seed %>%
      dplyr::left_join(
        hh_seed %>% dplyr::select(hhid, geo_colname),
        by = "hhid"
      )
    
    # Get vector of other column names
    col_names <- colnames(tbl)
    col_names <- type.convert(col_names[!col_names == geo_colname], as.is = TRUE)
    
    for (geo in per_seed[, geo_colname]){
      test <- match(col_names, per_seed[[name]][per_seed[, geo_colname] == geo])
      if (any(is.na(test))) {
        prob_cat <- col_names[which(is.na(test))]
        stop(
          "Marginal ", name, "; category ", prob_cat[1], " is missing from ",
          geo_colname, " ", geo, " in the per_seed table."
        )
      }   
    }
  }
}


