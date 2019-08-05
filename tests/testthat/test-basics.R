test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("ipu_matrix works", {
  set.seed(42)
  mtx <- matrix(data = runif(9), nrow = 3, ncol = 3)
  row_targets <- c(3, 4, 5)
  column_targets <- c(5, 4, 3)
  result <- ipu_matrix(mtx, row_targets, column_targets)
  expect_equal(class(result), "matrix")
  expect_equal(round(rowSums(result)[1], 4), 3)
  expect_equal(colSums(result)[3], 3)
})

test_that("basic ipu works", {
  result <- setup_arizona()
  hh_seed <- result$hh_seed
  hh_targets <- result$hh_targets
  per_seed <- result$per_seed
  per_targets <- result$per_targets
  
  result <- ipu(hh_seed, hh_targets, per_seed, per_targets, max_iterations = 30)
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_equal(
    names(result),
    c("weight_tbl", "weight_dist", "primary_comp", "secondary_comp")
  )
  expect_equal(round(result$weight_tbl$weight[1], 3), 3.769)
})

test_that("weight constraint works", {
  result <- setup_arizona()
  hh_seed <- result$hh_seed
  hh_targets <- result$hh_targets
  per_seed <- result$per_seed
  per_targets <- result$per_targets
  
  min_ratio <- .2
  max_ratio <- 1.2
  result <- ipu(hh_seed, hh_targets, per_seed, per_targets, max_iterations = 30,
                min_ratio = min_ratio, max_ratio = max_ratio)
  expect_true(max(result$weight_tbl$weight_factor) == max_ratio)
  expect_true(min(result$weight_tbl$weight_factor) == min_ratio)
})

test_that("secondary_importance works", {
  result <- setup_arizona()
  hh_seed <- result$hh_seed
  hh_targets <- result$hh_targets
  per_seed <- result$per_seed
  per_targets <- result$per_targets
  
  result <- ipu(hh_seed, hh_targets, per_seed, per_targets,
      secondary_importance = .5, max_iterations = 10,
      verbose = TRUE)
  expect_equal(round(result$secondary_comp$pct_diff[1], 2), -.41)
})

test_that("single value marginals work", {
  result <- setup_arizona()
  hh_seed <- result$hh_seed
  hh_targets <- result$hh_targets
  # per_seed <- result$per_seed
  # per_targets <- result$per_targets
  
  hh_seed$hhtype <- 1
  hh_targets$hhtype$`2` <- NULL
  result <- ipu(hh_seed, hh_targets)
  expect_type(ipu(hh_seed, hh_targets), "list")
})