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
  args <- setup_arizona()
  result <- ipu(args[[1]], args[[2]], args[[3]], args[[4]], max_iterations = 30)
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_equal(
    names(result),
    c("weight_tbl", "weight_dist", "primary_comp", "secondary_comp")
  )
  expect_equal(round(result$weight_tbl$weight[1], 3), 3.769)
})

test_that("weight constraint works", {
  args <- setup_arizona()
  min_ratio <- .2
  max_ratio <- 1.2
  result <- ipu(args[[1]], args[[2]], args[[3]], args[[4]], max_iterations = 30,
                min_ratio = min_ratio, max_ratio = max_ratio)
  expect_true(max(result$weight_tbl$weight_factor) == max_ratio)
  expect_true(min(result$weight_tbl$weight_factor) == min_ratio)
})