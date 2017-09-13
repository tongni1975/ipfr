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
  
  # Check check that seed and target are provided
  if (is.null(seed)) {
    stop("Seed table not provided")
  }
  if (is.null(targets)) {
    stop("Targets not provided")
  }
  
  # Check that there are no NA values in seed or targets
  if (any(is.na(unlist(seed)))) {
    stop("'seed' contains NAs")
  }
  if (any(is.na(unlist(targets)))) {
    stop("'targets' contains NAs")
  }
  
  # If a 'cluster' column wasn't provided on hh_seed, then repeat hh_seed
  # for every cluster found in the hh_targets.
  if (!("cluster" %in% colnames(hh_seed))) {
    hh_seed_mod <- merge(hh_targets[[1]]$cluster, hh_seed) %>%
      dplyr::rename(cluster = x) %>%
      dplyr::arrange(cluster, hhid)
  } else {
    hh_seed_mod <- hh_seed
  }
  # Throw an error if the cluster is defined on the person seed
  if ("cluster" %in% colnames(per_seed)) {
    stop(
      "A 'cluster' column exists on the 'per_seed' table. Only specify on 'hh_seed'."
    )
  }
  
  # Modify the household seed to the required format. Use one-hot-encoding to
  # expand standard columns into dummy columns.
  col_names <- names(hh_targets)
  hh_seed_mod <- hh_seed_mod %>%
    # Keep only the fields of interest
    select(one_of(c(col_names, "cluster", "hhid"))) %>%
    dplyr::mutate_at(
      .vars = col_names,
      .funs = funs(as.factor(.))
    ) %>%
    mlr::createDummyFeatures()
  
  # Modify the person seed table the same way, but sum by household ID
  col_names <- names(per_targets)
  per_seed_mod <- per_seed %>%
    # Keep only the fields of interest
    select(one_of(c(col_names, "hhid"))) %>%
    dplyr::mutate_at(
      .vars = col_names,
      .funs = funs(as.factor(.))
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
  attribute_cols <- colnames(seed)[-1]
  attribute_cols <- attribute_cols[-length(attribute_cols)]
  
  # modify the targets to match the new seed column names
  targets <- NULL
  for (name in names(hh_targets)) {
    temp <- hh_targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -cluster) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target"))
    targets <- dplyr::bind_rows(targets, temp)
  }
  for (name in names(per_targets)) {
    temp <- per_targets[[name]] %>%
      tidyr::gather(key = "key", value = "target", -cluster) %>%
      dplyr::mutate(key = paste0(!!name, ".", key, ".target"))
    targets <- dplyr::bind_rows(targets, temp)
  }
  targets <- targets %>%
    tidyr::spread(key = key, value = target)
  
  # join targets to the seed table
  # final <- seed %>%
  #   dplyr::left_join(cluster_equiv, by = "hhid") %>%
  #   dplyr::left_join(targets, by = "cluster") %>%
  #   dplyr::group_by(cluster)
  
  iter <- 1
  converged <- FALSE
  while (!converged & iter <= max_iterations) {
    # Loop over each target and upate weights
    for (attribute in attribute_cols) {
      target_col <- paste0(attribute, ".", "target")
      
      to_join <- final %>%
        dplyr::filter((!!as.name(attribute)) > 0) %>%
        dplyr::select(
          cluster, hhid, attr = !!attribute, target = !!target_col, weight
        ) %>%
        dplyr::mutate(factor = (target) / sum(attr * weight)) %>%
        ungroup() %>%
        dplyr::select(hhid, factor)
      
      final <- final %>%
        dplyr::left_join(to_join, by = "hhid") %>%
        dplyr::mutate(weight = ifelse(!is.na(factor), weight * factor, weight)) %>%
        # Implement the floor on minimum weight
        dplyr::mutate(weight = pmax(weight, min_weight)) %>%
        dplyr::select(-factor)
    }
    
    # Determine relative gaps (by cluster)
    gap <- 0
    for (attribute in attribute_cols) {
      target_col <- paste0(attribute, ".", "target")
      
      gap_tbl <- final %>%
        dplyr::filter((!!as.name(attribute)) > 0) %>%
        dplyr::select(cluster, hhid, attr = !!attribute, target = !!target_col, weight) %>%
        dplyr::mutate(gap = (target - sum(attr * weight)) / target) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      
      gap <- ifelse(max(gap_tbl$gap) > gap, max(gap_tbl$gap), gap)
    }
    
    # Test for convergence
    converged <- ifelse(gap <= relative_gap, TRUE, FALSE)
    iter <- iter + 1
  }
  
  hh_seed$weight <- final$weight
  
  return(hh_seed)
  
}