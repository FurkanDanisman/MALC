G1 <- function(y) stats::pnorm(y)
G2 <- function(y) stats::dnorm(y)

EMemp_counts <- function(counts, grid, sigma = 1, start = 0,
                         max_step = 1000, eps2 = 1e-10, eps1 = 1e-5) {
  .check_counts_grid(counts, grid)
  n <- sum(counts)
  pn <- counts / n

  mu_new <- start
  mu_vec <- mu_new

  for (i in seq_len(max_step)) {
    mu_old <- mu_new
    alpha <- (grid - mu_old) / sigma
    temp <- (diff(G2(alpha)) + eps2) / (diff(G1(alpha)) + eps2)
    mu_new <- mu_old - sigma * sum(pn * temp)
    mu_vec <- c(mu_new, mu_vec)

    if (abs(mu_new - mu_old) < eps1) {
      break
    }
  }

  list(mu_hat = mu_new, mu_vec = mu_vec)
}

MALC <- function(counts, grid, smooth = FALSE, B = 10000, alpha = 2,
                 seed = 20180621, log_conc = TRUE) {
  .check_counts_grid(counts, grid)
  .check_uniform_grid(grid)

  counts <- as.integer(round(counts))
  n <- sum(counts)
  pn <- counts / n
  delta <- diff(grid)[1L]
  g0 <- grid[1L]
  grid_left <- grid[-length(grid)]

  y <- rep(grid_left, times = counts)
  y_int <- rep(seq_along(counts) - 1L, times = counts)

  mle1 <- logcondiscr::logConDiscrMLE(y_int, output = FALSE)
  phat <- exp(mle1$psi) / sum(exp(mle1$psi))

  mu_low <- sum(pn * grid_left)
  mu_high <- sum(pn * grid[-1L])
  mu_mid <- (mu_low + mu_high) / 2
  sigma <- stats::sd(y)
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- sqrt(sum(pn * (.left_centers(grid) - mu_mid)^2) + delta^2 / 12)
  }

  em_fit <- EMemp_counts(counts, grid, start = mu_mid, sigma = sigma)
  mu_n <- em_fit$mu_hat
  Delta <- mu_n - mu_low
  beta <- 2 * alpha * (Delta / delta - 0.5)
  check_z <- min(alpha + beta, alpha - beta)
  if (!is.finite(check_z) || check_z <= 0) {
    stop("MALC beta perturbation parameters are not positive for this input.", call. = FALSE)
  }

  old_seed <- .restore_seed(seed)
  on.exit(.on_exit_restore_seed(old_seed), add = TRUE)

  ystar <- sample(mle1$x, B, prob = phat, replace = TRUE)
  ystar <- g0 + ystar * delta
  zstar <- delta * stats::rbeta(B, alpha + beta, alpha - beta)
  xstar <- ystar + zstar

  mle2 <- logcondens::logConDens(xstar, smoothed = log_conc, print = FALSE)

  structure(
    list(
      fhatn = mle2,
      counts = counts,
      grid = grid,
      smooth = smooth,
      B = B,
      alpha = alpha,
      beta = beta,
      mu_hat = mu_n,
      sumphat = sum(phat),
      checkZ = check_z
    ),
    class = "MALC"
  )
}

.density_MALC <- function(object, x) {
  if (!inherits(object, "MALC")) {
    stop("object must be a MALC fit.", call. = FALSE)
  }
  if (isTRUE(object$smooth)) {
    out <- logcondens::evaluateLogConDens(x, object$fhatn, which = 4)
    ans <- out[, 5]
  } else {
    out <- logcondens::evaluateLogConDens(x, object$fhatn, which = 2)
    ans <- out[, 3]
  }
  ans[is.na(ans)] <- 0
  ans
}

dmalc <- function(object, x) {
  .density_MALC(object, x)
}

pmalc <- function(object, q, length_grid = 1001) {
  eval_grid <- .eval_grid_from_breaks(object$grid, length_grid = length_grid)
  dens <- dmalc(object, eval_grid)
  dx <- diff(eval_grid)
  cdf <- c(0, cumsum((dens[-length(dens)] + dens[-1L]) * dx / 2))
  if (max(cdf) > 0) {
    cdf <- cdf / max(cdf)
  }
  stats::approx(eval_grid, cdf, xout = q, yleft = 0, yright = 1, rule = 2)$y
}

qmalc <- function(object, p, length_grid = 1001) {
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("p must be finite probabilities in [0, 1].", call. = FALSE)
  }
  eval_grid <- .eval_grid_from_breaks(object$grid, length_grid = length_grid)
  dens <- dmalc(object, eval_grid)
  dx <- diff(eval_grid)
  cdf <- c(0, cumsum((dens[-length(dens)] + dens[-1L]) * dx / 2))
  if (max(cdf) > 0) {
    cdf <- cdf / max(cdf)
  }
  stats::approx(cdf, eval_grid, xout = p, ties = "ordered", rule = 2)$y
}

mean_malc <- function(object) {
  if (!inherits(object, "MALC")) {
    stop("object must be a MALC fit.", call. = FALSE)
  }
  object$mu_hat
}

plot.MALC <- function(x, eval_grid = NULL, add = FALSE,
                      xlab = "x", ylab = "Density", col = "navy",
                      lwd = 2, ...) {
  if (is.null(eval_grid)) {
    eval_grid <- .eval_grid_from_breaks(x$grid)
  }
  dens <- dmalc(x, eval_grid)
  if (add) {
    lines(eval_grid, dens, col = col, lwd = lwd, ...)
  } else {
    plot(eval_grid, dens, type = "l", xlab = xlab, ylab = ylab,
         col = col, lwd = lwd, ...)
  }
  invisible(data.frame(x = eval_grid, density = dens))
}

Plot <- function(object, ...) {
  plot(object, ...)
}
