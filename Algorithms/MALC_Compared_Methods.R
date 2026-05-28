library(KernSmooth)
library(binnednp)

# Helper: validate bins 
.validate_binned_input <- function(counts, breaks) {
  if (!is.numeric(counts) || !is.numeric(breaks)) stop("counts and breaks must be numeric.")
  if (length(breaks) != length(counts) + 1) stop("Need length(breaks) = length(counts)+1.")
  if (any(diff(breaks) <= 0)) stop("breaks must be strictly increasing.")
  if (any(counts < 0) || any(!is.finite(counts))) stop("counts must be finite and >= 0.")
  if (sum(counts) <= 0) stop("sum(counts) must be > 0.")
  invisible(TRUE)
}

# Helper: validate a uniform, strictly increasing grid 
.validate_uniform_grid <- function(grid, name = "grid", rtol = 1e-8, atol = 1e-12) {
  if (!is.numeric(grid) || length(grid) < 2) stop(sprintf("%s must be numeric with length >= 2.", name))
  if (any(!is.finite(grid))) stop(sprintf("%s must be finite.", name))
  d <- diff(grid)
  if (any(d <= 0)) stop(sprintf("%s must be strictly increasing.", name))
  
  # robust uniformity check: compare to median step
  dx <- stats::median(d)
  dev <- max(abs(d - dx))
  
  # allow both absolute and relative wiggle
  if (dev > atol + rtol * max(1, abs(dx))) {
    stop(sprintf(
      "%s must be (approximately) uniformly spaced. max|diff-Δ| = %.3e, Δ = %.6g (rtol=%g, atol=%g)",
      name, dev, dx, rtol, atol
    ))
  }
  invisible(TRUE)
}

# Paper-suggested bandwidth
simulate_pseudodata_uniform_in_bins <- function(counts, breaks, seed = NULL) {
  .validate_binned_input(counts, breaks)
  # if (!is.null(seed)) set.seed(seed)
  
  m <- length(counts)
  n <- sum(counts)
  
  x <- numeric(n)
  idx <- 1L
  
  for (j in seq_len(m)) {
    nj <- counts[j]
    if (nj == 0) next
    a <- breaks[j]
    b <- breaks[j + 1]
    x[idx:(idx + nj - 1L)] <- runif(nj, min = a, max = b)
    idx <- idx + nj
  }
  x
}

# SJ bandwidth for binned data (median over repetitions) 
bw_SJ_from_binned <- function(counts, breaks, B = 50, seed = 1, bw_method = c("ste", "dpi")) {
  .validate_binned_input(counts, breaks)
  bw_method <- match.arg(bw_method)
  
  hs <- numeric(B)
  for (b in seq_len(B)) {
    x_pseudo <- simulate_pseudodata_uniform_in_bins(counts, breaks, seed = seed + b - 1L)
    hs[b] <- stats::bw.SJ(x_pseudo, method = bw_method)
  }
  stats::median(hs, na.rm = TRUE)
}

# FFT convolution
fft_convolve_real_open <- function(a, b) {
  n <- length(a) + length(b) - 1L
  nfft <- 2L^ceiling(log2(n))
  A <- fft(c(a, rep(0, nfft - length(a))))
  B <- fft(c(b, rep(0, nfft - length(b))))
  out <- Re(fft(A * B, inverse = TRUE)) / nfft
  out[seq_len(n)]
}

# Main estimator: fixed point of T 

# The algorithm evaluates the density on bin CENTERS derived from grid.
nonlinear_kde_binned_BK2002 <- function(counts,
                                        breaks,
                                        h = NULL,
                                        bw_B = 50,
                                        bw_seed = 1,
                                        bw_method = c("ste", "dpi"),
                                        grid,               # interpreted as breaks you want to use
                                        extend = 0,         # keep 0 if you want to *stick* to your grid
                                        max_iter = 200,
                                        tol_L1 = 1e-6,
                                        verbose = TRUE) {
  .validate_binned_input(counts, breaks)
  bw_method <- match.arg(bw_method)
  
  # Interpret the grid as breaks (edges) and build centers 
  .validate_uniform_grid(grid = grid, name = "grid (breaks)")
  grid_breaks <- grid
  x_grid <- (grid_breaks[-1] + grid_breaks[-length(grid_breaks)]) / 2  # centers
  dx <- x_grid[2] - x_grid[1]
  grid_n <- length(x_grid)
  
  N <- sum(counts)
  p <- counts / N
  widths <- diff(breaks)
  m <- length(counts)
  
  # Bandwidth selection
  if (is.null(h)) {
    h <- bw_SJ_from_binned(counts, breaks, B = bw_B, seed = bw_seed, bw_method = bw_method)
  }
  if (!is.finite(h) || h <= 0) stop("Bandwidth h must be positive and finite.")
  
  # Initial f0: histogram density on A, 0 outside
  # NOTE: Bin densities assigned to CENTER points that fall in each original break-interval.
  f <- numeric(grid_n)
  for (j in seq_len(m)) {
    in_bin <- (x_grid >= breaks[j]) & (x_grid < breaks[j + 1])
    if (j == m) in_bin <- (x_grid >= breaks[j]) & (x_grid <= breaks[j + 1])
    if (any(in_bin)) f[in_bin] <- counts[j] / (N * widths[j])
  }
  
  # Normalize so sum(f)*dx ≈ 1
  mass0 <- sum(f) * dx
  if (!is.finite(mass0) || mass0 <= 0) stop("Initialization produced zero/invalid mass. Check that x_grid overlaps breaks.")
  f <- f / mass0
  
  # Precompute discrete Gaussian kernel weights for convolution on x_grid (uniform spacing dx)
  lags <- (-(grid_n - 1L)):(grid_n - 1L)
  k_vec <- stats::dnorm(lags * dx, mean = 0, sd = h) * dx
  
  # Bin masks on x_grid for Aj restriction
  bin_masks <- vector("list", m)
  for (j in seq_len(m)) {
    mask <- (x_grid >= breaks[j]) & (x_grid < breaks[j + 1])
    if (j == m) mask <- (x_grid >= breaks[j]) & (x_grid <= breaks[j + 1])
    bin_masks[[j]] <- mask
  }
  
  # Iterate f_{k+1} = T f_k
  for (iter in seq_len(max_iter)) {
    f_next <- numeric(grid_n)
    
    for (j in seq_len(m)) {
      if (p[j] == 0) next
      
      mask <- bin_masks[[j]]
      g_j <- f
      g_j[!mask] <- 0
      
      denom_j <- sum(g_j) * dx
      if (!is.finite(denom_j) || denom_j <= 0) {
        stop(sprintf("Denominator integral over bin %d is non-positive. This usually means x_grid has no points in that bin. Use centers (this code does), or widen/extend grid.", j))
      }
      
      # numerator via convolution
      conv <- fft_convolve_real_open(g_j, k_vec)
      num_j <- conv[grid_n:(2L * grid_n - 1L)]
      
      f_next <- f_next + p[j] * (num_j / denom_j)
    }
    
    # Normalize
    total_mass <- sum(f_next) * dx
    if (!is.finite(total_mass) || total_mass <= 0) stop("Iteration produced invalid density mass.")
    f_next <- f_next / total_mass
    
    # L1 convergence
    diff_L1 <- sum(abs(f_next - f)) * dx
    if (verbose) message(sprintf("iter=%d  L1=%.3e  h=%.6g", iter, diff_L1, h))
    
    f <- f_next
    if (diff_L1 < tol_L1) break
  }
  
  list(
    x = x_grid,           # centers
    fhat = f,             # density at centers
    h = h,
    counts = counts,
    breaks = breaks,      # original histogram breaks used for binning counts
    eval_breaks = grid_breaks, # breaks that define the evaluation grid 
    dx = dx
  )
}


# Convenience plotting 

plot_binned_kde <- function(res, add_hist = TRUE, main = "Nonlinear binned KDE (BK2002)", ...) {
  x <- res$x
  fhat <- res$fhat
  plot(x, fhat, type = "l", main = main, xlab = "x", ylab = "density", ...)
  
  if (add_hist) {
    # overlay histogram density as step function
    counts <- res$counts
    breaks <- res$breaks
    N <- sum(counts)
    widths <- diff(breaks)
    hden <- counts / (N * widths)
    
    for (j in seq_along(counts)) {
      lines(c(breaks[j], breaks[j+1]), c(hden[j], hden[j]))
      if (j < length(counts)) {
        lines(c(breaks[j+1], breaks[j+1]), c(hden[j], hden[j+1]))
      }
    }
  }
}

binned_np.func <- function(x, grid, length_grid = 1001, x_eval = NULL,
                           seed_start = 1L,
                           max_tries = 200L){
  
  counts <- hist(x, breaks = grid, plot = FALSE)$counts
  n      <- sum(counts)
  w      <- counts / n
  
  t_mid  <- (grid[-length(grid)] + grid[-1]) / 2
  if (is.null(x_eval)) {
    x_eval <- seq(min(grid),max(grid), length.out=length_grid)
  }
  
  fhat_g <- function(x, h) {
    sapply(x, function(xx) sum(w * dnorm((xx - t_mid) / h)) / h)
  }
  
  # retry bw.dens.binned if fails
  old_seed <- .Random.seed
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  
  bw <- NULL
  for (k in 0:(max_tries - 1L)) {
    set.seed(seed_start + k)
    
    bw_try <- tryCatch(
      binnednp::bw.dens.binned(
        n = n, y = grid, w = w, ni = counts,
        plot = FALSE, print = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(bw_try) &&
        is.finite(bw_try$h_boot) && bw_try$h_boot > 0 &&
        is.finite(bw_try$h_plugin) && bw_try$h_plugin > 0) {
      bw <- bw_try
      break
    }
  }
  
  if (is.null(bw)) {
    stop("bw.dens.binned failed after max_tries. Increase max_tries or inspect counts/grid.")
  }
  
  f_boot <- fhat_g(x_eval, bw$h_boot)
  return(f_boot)
}

binned_np.func_for_plot = function(x,grid,grid2){
  
  counts <- hist(x, breaks = grid, plot = FALSE)$counts
  n      <- sum(counts)
  w      <- counts / n
  
  t_mid  <- (grid[-length(grid)] + grid[-1]) / 2  # t_i in the paper
  x_eval <- grid2  # evaluate density on *bin values* if you want
  
  fhat_g <- function(x, h) {
    # Gaussian kernel K(u)=dnorm(u)
    sapply(x, function(xx) sum(w * dnorm((xx - t_mid) / h)) / h)
  }
  
  bw <- bw.dens.binned(n = n, y = grid, w = w, ni = counts, 
                       plot = FALSE, print = FALSE)
  
  # f_plugin <- fhat_g(x_eval, bw$h_plugin)
  f_boot   <- fhat_g(x_eval, bw$h_boot)
  
  return(f_boot)
  
}

library(KernSmooth)

dpik_with_fallback_jitter <- function(y, counts, breaks, grid,
                                      seed = 1,
                                      scalest = "minim",
                                      level = 2L,
                                      kernel = "normal",
                                      canonical = FALSE,
                                      truncate = TRUE) {
  
  # gridsize
  gs <- length(grid)
  
  try_dpik <- function(y_in) {
    dpik(y_in,
         scalest   = scalest,
         level     = level,
         kernel    = kernel,
         canonical = canonical,
         gridsize  = gs,
         range.x   = range(grid),
         truncate  = truncate)
  }
  
  h <- tryCatch(
    try_dpik(y),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("scale estimate is zero", msg, fixed = TRUE)) {
        # set.seed(seed)
        y2 <- simulate_pseudodata_uniform_in_bins(counts, breaks)
        return(try_dpik(y2))
      } else {
        stop(e)  # other errors should surface
      }
    }
  )
  
  if (!is.finite(h) || h <= 0) {
    # set.seed(seed)
    y2 <- simulate_pseudodata_uniform_in_bins(counts, breaks)
    h <- try_dpik(y2)
  }
  
  h
}
