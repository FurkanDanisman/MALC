simulate_pseudodata_uniform_in_bins <- function(counts, breaks, seed = NULL) {
  .check_counts_grid(counts, breaks)
  counts <- as.integer(round(counts))
  old_seed <- .restore_seed(seed)
  on.exit(.on_exit_restore_seed(old_seed), add = TRUE)

  x <- numeric(sum(counts))
  idx <- 1L
  for (j in seq_along(counts)) {
    nj <- counts[j]
    if (nj == 0L) {
      next
    }
    x[idx:(idx + nj - 1L)] <- stats::runif(nj, min = breaks[j], max = breaks[j + 1L])
    idx <- idx + nj
  }
  x
}

bw_SJ_from_binned <- function(counts, breaks, B = 50, seed = 1,
                              bw_method = c("ste", "dpi")) {
  .check_counts_grid(counts, breaks)
  bw_method <- match.arg(bw_method)
  hs <- numeric(B)
  for (b in seq_len(B)) {
    x_pseudo <- simulate_pseudodata_uniform_in_bins(counts, breaks, seed = seed + b - 1L)
    hs[b] <- stats::bw.SJ(x_pseudo, method = bw_method)
  }
  stats::median(hs, na.rm = TRUE)
}

fft_convolve_real_open <- function(a, b) {
  n <- length(a) + length(b) - 1L
  nfft <- 2L^ceiling(log2(n))
  A <- stats::fft(c(a, rep(0, nfft - length(a))))
  B <- stats::fft(c(b, rep(0, nfft - length(b))))
  out <- Re(stats::fft(A * B, inverse = TRUE)) / nfft
  out[seq_len(n)]
}

BK2002 <- function(counts, grid, h = NULL, bw_B = 50, bw_seed = 1,
                   bw_method = c("ste", "dpi"), max_iter = 200,
                   tol_L1 = 1e-6, verbose = FALSE) {
  .check_counts_grid(counts, grid)
  .check_uniform_grid(grid)
  bw_method <- match.arg(bw_method)

  x_grid <- .left_centers(grid)
  dx <- x_grid[2L] - x_grid[1L]
  grid_n <- length(x_grid)
  n <- sum(counts)
  p <- counts / n
  widths <- diff(grid)
  m <- length(counts)

  if (is.null(h)) {
    h <- bw_SJ_from_binned(counts, grid, B = bw_B, seed = bw_seed, bw_method = bw_method)
  }
  if (!is.finite(h) || h <= 0) {
    stop("Bandwidth h must be positive and finite.", call. = FALSE)
  }

  f <- counts / (n * widths)
  mass0 <- sum(f) * dx
  if (!is.finite(mass0) || mass0 <= 0) {
    stop("Initialization produced invalid density mass.", call. = FALSE)
  }
  f <- f / mass0

  lags <- (-(grid_n - 1L)):(grid_n - 1L)
  k_vec <- stats::dnorm(lags * dx, mean = 0, sd = h) * dx

  bin_masks <- vector("list", m)
  for (j in seq_len(m)) {
    mask <- (x_grid >= grid[j]) & (x_grid < grid[j + 1L])
    if (j == m) {
      mask <- (x_grid >= grid[j]) & (x_grid <= grid[j + 1L])
    }
    bin_masks[[j]] <- mask
  }

  for (iter in seq_len(max_iter)) {
    f_next <- numeric(grid_n)
    for (j in seq_len(m)) {
      if (p[j] == 0) {
        next
      }
      g_j <- f
      g_j[!bin_masks[[j]]] <- 0
      denom_j <- sum(g_j) * dx
      if (!is.finite(denom_j) || denom_j <= 0) {
        stop(sprintf("Denominator integral over bin %d is non-positive.", j), call. = FALSE)
      }
      conv <- fft_convolve_real_open(g_j, k_vec)
      num_j <- conv[grid_n:(2L * grid_n - 1L)]
      f_next <- f_next + p[j] * (num_j / denom_j)
    }

    total_mass <- sum(f_next) * dx
    if (!is.finite(total_mass) || total_mass <= 0) {
      stop("Iteration produced invalid density mass.", call. = FALSE)
    }
    f_next <- f_next / total_mass

    diff_L1 <- sum(abs(f_next - f)) * dx
    if (verbose) {
      message(sprintf("iter=%d  L1=%.3e  h=%.6g", iter, diff_L1, h))
    }
    f <- f_next
    if (diff_L1 < tol_L1) {
      break
    }
  }

  structure(
    list(x = x_grid, fhat = f, h = h, counts = counts, grid = grid, dx = dx),
    class = "BK2002"
  )
}

.density_BK2002 <- function(object, x) {
  if (!inherits(object, "BK2002")) {
    stop("object must be a BK2002 fit.", call. = FALSE)
  }
  stats::approx(object$x, object$fhat, xout = x, rule = 2)$y
}

BinnedNP <- function(counts, grid, bandwidth = c("boot", "plugin"),
                     seed_start = 1L, max_tries = 200L) {
  .check_counts_grid(counts, grid)
  bandwidth <- match.arg(bandwidth)
  counts <- as.integer(round(counts))
  n <- sum(counts)
  w <- counts / n

  bw <- NULL
  for (k in 0:(max_tries - 1L)) {
    old_seed <- .restore_seed(seed_start + k)
    bw_try <- tryCatch(
      binnednp::bw.dens.binned(
        n = n, y = grid, w = w, ni = counts,
        plot = FALSE, print = FALSE
      ),
      error = function(e) NULL
    )
    .on_exit_restore_seed(old_seed)
    if (!is.null(bw_try) &&
        is.finite(bw_try$h_boot) && bw_try$h_boot > 0 &&
        is.finite(bw_try$h_plugin) && bw_try$h_plugin > 0) {
      bw <- bw_try
      break
    }
  }
  if (is.null(bw)) {
    stop("bw.dens.binned failed after max_tries.", call. = FALSE)
  }

  structure(
    list(
      counts = counts,
      grid = grid,
      weights = w,
      mids = .left_centers(grid),
      bw = bw,
      bandwidth = bandwidth
    ),
    class = "BinnedNP"
  )
}

.density_BinnedNP <- function(object, x) {
  if (!inherits(object, "BinnedNP")) {
    stop("object must be a BinnedNP fit.", call. = FALSE)
  }
  h <- if (object$bandwidth == "boot") object$bw$h_boot else object$bw$h_plugin
  vapply(
    x,
    function(xx) sum(object$weights * stats::dnorm((xx - object$mids) / h)) / h,
    numeric(1)
  )
}

KernSmooth <- function(counts, grid, seed = 1, scalest = "minim",
                       level = 2L, kernel = "normal", canonical = FALSE,
                       truncate = TRUE) {
  .check_counts_grid(counts, grid)
  .check_uniform_grid(grid)
  counts <- as.integer(round(counts))
  centers <- .left_centers(grid)
  y <- rep(centers, times = counts)

  try_dpik <- function(y_in) {
    KernSmooth::dpik(
      y_in,
      scalest = scalest,
      level = level,
      kernel = kernel,
      canonical = canonical,
      gridsize = length(grid),
      range.x = range(grid),
      truncate = truncate
    )
  }

  h <- tryCatch(
    try_dpik(y),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("scale estimate is zero", msg, fixed = TRUE)) {
        y2 <- simulate_pseudodata_uniform_in_bins(counts, grid, seed = seed)
        try_dpik(y2)
      } else {
        stop(e)
      }
    }
  )

  if (!is.finite(h) || h <= 0) {
    y2 <- simulate_pseudodata_uniform_in_bins(counts, grid, seed = seed)
    h <- try_dpik(y2)
  }

  fit <- KernSmooth::bkde(
    y,
    bandwidth = h,
    kernel = kernel,
    canonical = canonical,
    gridsize = length(grid),
    range.x = range(grid),
    truncate = truncate
  )

  structure(
    list(x = fit$x, fhat = fit$y, h = h, counts = counts, grid = grid),
    class = "DensOLogKernSmooth"
  )
}

.density_KernSmooth <- function(object, x) {
  if (!inherits(object, "DensOLogKernSmooth")) {
    stop("object must be a KernSmooth fit returned by DensOLog::KernSmooth().", call. = FALSE)
  }
  stats::approx(object$x, object$fhat, xout = x, rule = 2)$y
}

.is_kernel_fit <- function(object) {
  inherits(object, "BK2002") ||
    inherits(object, "BinnedNP") ||
    inherits(object, "DensOLogKernSmooth")
}

.density_kernel <- function(object, x) {
  if (inherits(object, "BK2002")) {
    .density_BK2002(object, x)
  } else if (inherits(object, "BinnedNP")) {
    .density_BinnedNP(object, x)
  } else if (inherits(object, "DensOLogKernSmooth")) {
    .density_KernSmooth(object, x)
  } else {
    stop("object must be a fit from BK2002(), BinnedNP(), or KernSmooth().", call. = FALSE)
  }
}

Plot_Kernel <- function(object, eval_grid = NULL, add = FALSE,
                        xlab = "x", ylab = "Density", col = "darkgreen",
                        lwd = 2, ...) {
  if (!.is_kernel_fit(object)) {
    stop("object must be a fit from BK2002(), BinnedNP(), or KernSmooth().", call. = FALSE)
  }
  if (is.null(eval_grid)) {
    eval_grid <- .eval_grid_from_breaks(object$grid)
  }
  dens <- .density_kernel(object, eval_grid)
  if (add) {
    lines(eval_grid, dens, col = col, lwd = lwd, ...)
  } else {
    plot(eval_grid, dens, type = "l", xlab = xlab, ylab = ylab,
         col = col, lwd = lwd, ...)
  }
  invisible(data.frame(x = eval_grid, density = dens))
}
