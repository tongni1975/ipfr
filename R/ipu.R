#' Iterative Proportional Updating
#' 
#' @description A special case of iterative proportional fitting when two,
#' disparate sets of marginals need to be satisfied that do not agree on a
#' single total. A common example is balance population data using household-
#' and person-level marginal controls. This could be for survey expansion or
#' synthetic population creation.
#' 
#' @references \url{http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf}
#' 
#' @inheritParams ipf
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
#' @return the \code{hh_seed} with a revised weight column.
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
ipu <- function(hh_seed, hh_targets, per_seed, per_targets,
                relative_gap = 0.01, absolute_gap = 1, max_iterations = 50,
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
  seed <- hh_seed_mod %>%
    dplyr::left_join(per_seed_mod, by = "hhid") %>%
    dplyr::mutate(weight = 1)
  
  # store a vector of attribute column names to loop over later.
  # don't include 'hhid' or 'weight' in the vector.
  hhid_pos <- grep("hhid", colnames(seed))
  weight_pos <- grep("weight", colnames(seed))
  seed_attribute_cols <- colnames(seed)[-c(hhid_pos, weight_pos)]
  
  # modify the targets to match the new seed column names
  targets <- c(hh_targets, per_targets)
  for (name in names(targets)) {
    targets[[name]] <- targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -dplyr::starts_with("geo_")) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target"))
  }
  
  # Create a copy of the seed table that will have the final results.
  # Add back in the geo information which will be needed for joining.
  final <- seed %>%
    dplyr::left_join(geo_equiv, by = "hhid")
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # Loop over each target and upate weights
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      target_tbl <- targets[[target_tbl_name]] %>%
        dplyr::filter(key == target_name)
      
      # Get the name of the geo column
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      # Join the target table to the seed table
      fac_tbl <- final %>%
        dplyr::filter((!!as.name(seed_attribute)) > 0) %>%
        dplyr::left_join(target_tbl, by = geo_colname) %>%
        dplyr::select(
          geo = !!geo_colname, hhid, attr = !!seed_attribute, weight, target
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(factor = (target) / sum(attr * weight)) %>%
        dplyr::ungroup() %>%
        dplyr::select(hhid, factor)
      
      final <- final %>%
        dplyr::left_join(fac_tbl, by = "hhid") %>%
        dplyr::mutate(weight = ifelse(!is.na(factor), weight * factor, weight)) %>%
        # Implement the floor on minimum weight
        dplyr::mutate(weight = pmax(weight, min_weight)) %>%
        dplyr::select(-factor)
    }
    
    # Determine relative gaps (by geo field)
    rel_gap <- 0
    for (seed_attribute in seed_attribute_cols) {
      # create lookups for targets list
      target_tbl_name <- strsplit(seed_attribute, ".", fixed = TRUE)[[1]][1]
      target_name <- paste0(seed_attribute, ".", "target")
      target_tbl <- targets[[target_tbl_name]] %>%
        dplyr::filter(key == target_name)
      
      # Get the name of the geo column
      pos <- grep("geo_", colnames(target_tbl))
      geo_colname <- colnames(target_tbl)[pos]
      
      gap_tbl <- final %>%
        dplyr::filter((!!as.name(seed_attribute)) > 0) %>%
        dplyr::left_join(target_tbl, by = geo_colname) %>%
        dplyr::select(
          geo = !!geo_colname, hhid, attr = !!seed_attribute, weight, target
        ) %>%
        dplyr::group_by(geo) %>%
        dplyr::mutate(
          abs_gap = abs((target - sum(attr * weight))),
          rel_gap = abs_gap / (target + .0000001) # avoid dividing by zero
        ) %>%
        # Removes rows where the absolute gap is smaller than 'absolute_gap'
        dplyr::filter(abs_gap > absolute_gap) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      
      # If any records are left in the gap_tbl, record worst relative gap and
      # save that gap table for reporting out.
      if (nrow(gap_tbl) > 0) {
        if (max(gap_tbl$rel_gap) > rel_gap) {
          rel_gap <- max(gap_tbl$rel_gap)
          saved_gap_tbl <- gap_tbl
          saved_category <- seed_attribute
          saved_geo <- geo_colname
        }
      }
      
    }
    
    # Test for convergence
    converged <- ifelse(rel_gap <= relative_gap, TRUE, FALSE)
    iter <- iter + 1
  }
  
  if (verbose) {
    position <- which(saved_gap_tbl$rel_gap == rel_gap)[1]
    message("Max Rel Gap:", rel_gap)
    message("Absolute Gap:", saved_gap_tbl$abs_gap[position])
    message(saved_geo, ": ", saved_gap_tbl$geo[position])
    message("Category:", saved_category)
    utils::flush.console()
  }
  
  hh_seed$weight <- final$weight
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