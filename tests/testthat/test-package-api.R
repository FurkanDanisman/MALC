test_that("comparison methods use counts and grid inputs", {
  counts <- c(1, 3, 5, 3, 1)
  grid <- seq(-2.5, 2.5, length.out = length(counts) + 1)
  x_eval <- seq(min(grid), max(grid), length.out = 25)

  do <- DensOLog(counts, grid, B = 200)
  expect_equal(mean_densolog(do), do$mu_hat)
  expect_true(is.finite(mean_densolog(do)))

  do_s <- DensOLog(counts, grid, B = 200, smooth = TRUE)
  expect_length(ddensolog(do, x_eval), length(x_eval))
  expect_length(ddensolog(do_s, x_eval), length(x_eval))
  expect_false(isTRUE(all.equal(ddensolog(do, x_eval), ddensolog(do_s, x_eval))))
  expect_true(is.finite(Eval(do, dnorm)))
  expect_true(is.finite(Eval(do_s, dnorm)))

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  bk <- BK2002(counts, grid, bw_B = 2, max_iter = 2)
  expect_true(is.finite(Eval(bk, dnorm)))
  expect_s3_class(Plot_Kernel(bk, eval_grid = x_eval), "data.frame")

  bn <- BinnedNP(counts, grid, bandwidth = "plugin", max_tries = 20)
  expect_true(is.finite(Eval(bn, dnorm)))
  expect_s3_class(Plot_Kernel(bn, eval_grid = x_eval), "data.frame")

  ks <- KernSmooth(counts, grid)
  expect_true(is.finite(Eval(ks, dnorm)))
  expect_s3_class(Plot_Kernel(ks, eval_grid = x_eval), "data.frame")
})
