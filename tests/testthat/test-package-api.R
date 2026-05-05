test_that("comparison methods use counts and grid inputs", {
  counts <- c(1, 3, 5, 3, 1)
  grid <- seq(-2.5, 2.5, length.out = length(counts) + 1)
  x_eval <- seq(min(grid), max(grid), length.out = 25)

  do <- DensOLog(counts, grid, B = 200)
  expect_equal(mean_densolog(do), do$mu_hat)
  expect_true(is.finite(mean_densolog(do)))

  bk <- BK2002(counts, grid, bw_B = 2, max_iter = 2)
  expect_length(evaluate_BK2002(bk, x_eval), length(x_eval))

  ks <- KernSmooth(counts, grid)
  expect_length(evaluate_KernSmooth(ks, x_eval), length(x_eval))
})
