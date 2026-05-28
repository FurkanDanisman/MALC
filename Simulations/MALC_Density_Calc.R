L2_Distance_calc_both = function(x,grid,pdf, length_grid = 1001){
  
  # Running the LCD Algorithm
  fhat2 <- get_fhatn(x, grid, alpha=2,log_conc = FALSE)
  
  est_grid <- seq(min(grid),max(grid), length.out=length_grid) # this is the evaluation grid
  del <- min(diff(est_grid)) # this is the difference going into the Reimann sum
  
  f0 <- pdf(est_grid)
  
  fhat <- evaluateLogConDens(est_grid, fhat2$fhatn, which=2)
  fhat_smooth <- evaluateLogConDens(est_grid, fhat2$fhatn, which=4)
  fhat_smooth[,5][is.na(fhat_smooth[,5])] = 0 
  
  L2        <- sum((f0-fhat[,3])^2)*del
  L2_smooth <- sum((f0-fhat_smooth[,5])^2)*del
  
  L2 = sqrt(L2)
  L2_smooth = sqrt(L2_smooth)
  
  c(L2,L2_smooth)
  
  return(c(L2,L2_smooth))
  
}

L2_baseline = function(x,grid,pdf, length_grid = 1001){
  
  est_grid <- seq(min(grid),max(grid), length.out=length_grid) # this is the evaluation grid
  del <- min(diff(est_grid)) # this is the difference going into the Reimann sum
  
  f0 <- pdf(est_grid)
  
  fhat <- 0
  
  L2 <- sum((f0-fhat)^2)*del
  L2 = sqrt(L2)
  
  return(L2)
  
}


L2_binnednp = function(x, grid, pdf, length_grid = 1001){
  
  # Running the LCD Algorithm
  
  est_grid <- seq(min(grid),max(grid), length.out=length_grid) # this is the evaluation grid
  estimated_grid_density <- binned_np.func(x, grid, x_eval = est_grid)
  true_density_values <- pdf(est_grid)
  delta = max(diff(est_grid))
  
  # Compute the absolute difference for L2 distance
  sqr_differance <- (true_density_values - estimated_grid_density)^2
  L2 = sum(sqr_differance)*delta
  L2 = sqrt(L2)
  
  return(L2)
  
}

L2_from_res <- function(x,grid, pdf, length_grid = 1001) {
  
  hobj <- hist(x, breaks = grid, plot = FALSE)
  est_grid <- seq(min(grid),max(grid), length.out=length_grid) # this is the evaluation grid
  
  res <- nonlinear_kde_binned_BK2002(
    counts = hobj$counts,
    breaks = hobj$breaks,
    grid = grid,
    h = NULL,        # <- forces SJ-from-binned as suggested in the paper
    bw_B = 50,
    bw_method = "ste",
    verbose = FALSE
  )
  
  fhat <- approx(res$x, res$fhat, xout = est_grid, rule = 2)$y
  dx <- est_grid[2] - est_grid[1]
  
  true <- pdf(est_grid)
  L2_mid <- sum((true - fhat)^2) * dx
  sqrt(L2_mid)
  
}

L2_from_kernsmooth <- function(x,grid, pdf, length_grid = 1001) {
  
  hobj <- hist(x, breaks=grid, plot=FALSE)
  counts <- hobj$counts
  centers <- (grid[-1] + grid[-length(grid)]) / 2
  est_grid <- seq(min(grid),max(grid), length.out=length_grid) # this is the evaluation grid
  
  y <- rep(centers, times=counts)  # expand at centers
  
  # 1) Bandwidth via dpik
  
  h <- dpik_with_fallback_jitter(
    y = y,
    counts = hobj$counts,
    breaks = hobj$breaks,
    grid = grid,
    seed = 1
  )
  
  # 2) KDE via bkde on the observed grid, then evaluate its output on est_grid
  fit <- bkde(y, bandwidth = h,
              kernel = "normal", canonical = FALSE,
              gridsize = length(grid), range.x = range(grid),
              truncate = TRUE)
  
  fhat <- approx(fit$x, fit$y, xout = est_grid, rule = 2)$y
  dx  <- est_grid[2] - est_grid[1]
  
  # 3) L2 on evaluation grid
  true <- pdf(est_grid)
  L2_mid <- sum((true - fhat)^2) * dx
  
  sqrt(L2_mid)
  
}
