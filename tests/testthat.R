library(testthat)
library(tidyverse)
library(ipfr)

# This code sets up the Arizona example, which avoids repeating this code
# in multiple tests.
setup_arizona <- function() {
  hh_seed <- tibble(
    id = c(1:8),
    hhtype = c(1, 1, 1, 2, 2, 2, 2, 2)
  )
  per_seed <- tibble(
    id = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8),
    pertype = c(1, 2, 3, 1, 3, 1, 1, 2, 1, 3, 3, 2, 2, 3, 1, 2, 1, 1, 2, 3, 3, 1, 2)
  )
  hh_targets <- list()
  hh_targets$hhtype <- tibble(
    `1` = 35,
    `2` = 65
  )
  per_targets <- list()
  per_targets$pertype <- tibble(
    `1` = 91,
    `2` = 65,
    `3` = 104
  )
  return(list(hh_seed, hh_targets, per_seed, per_targets))
}

test_check("ipfr")

# To run the tests manually
# devtools::test()