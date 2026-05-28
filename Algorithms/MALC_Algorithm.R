# Jan 17 2025

library(fdrtool)
library(truncnorm)
library(pracma)
library(logcondens)
library(logcondiscr)

#----------------------------

G1 <- function(y) pnorm(y)
G2 <- function(y) dnorm(y)

EMemp <- function(x, grid, sigma=1, start=0, max_step=1000, eps2=1e-10, eps1=1e-5){
  
  n <- length(x)
  pn <- (hist(x, breaks=grid,plot = FALSE)$counts)/n
  
  mu_new <- start
  mu_vec <- mu_new
  delta  <- 10
  
  for(i in 1:max_step){
    
    mu_old <- mu_new
    
    alpha <- (grid-mu_old)/sigma
    temp  <- (diff(G2(alpha))+eps2)/(diff(G1(alpha))+eps2)
    
    mu_new <- mu_old - sigma*sum(pn*temp)
    mu_vec <- c(mu_new, mu_vec)
    
    delta <- abs(mu_new-mu_old)
    
    if(delta<eps1) break
    
  }
  
  return(list(mu_hat = mu_new, mu_vec=mu_vec))
  
}

EMemp_with_n = function(n_i,grid){
  
  delta = max(diff(grid)) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
  # (Assuming a uniform grid, delta remains the same for any i)
  
  n = sum(n_i)
  pn <- n_i/n
  a_i <- grid[-length(grid)]
  m_bar <- sum(a_i * pn)  
  y = rep(a_i,n_i)
  
  interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
  mid_point <- sum(interval_centers * pn) 
  weighted_variance <- sum(pn * (interval_centers - mid_point)^2)
  
  mu_hat = EMemp(n_i,grid,start = mid_point,sigma = sqrt(weighted_variance))$mu_hat
  
  return(mu_hat)
  
}

#-----------
# get_fhatn 


get_fhatn_prior <- function(x, grid, B=10000, alpha=2, seednum=20180621){
  
  n    <- length(x)
  pn   <- (hist(x, breaks=grid, plot=FALSE)$counts)/n
  y    <- rep(grid, times=c(round(n*pn),0))
  
  delta <- max(diff(grid))
  
  y = round(y/delta) # ADDED BY FURKAN - Turning y variables into integers (Way around for logConDiscMLE function to work)
  
  mle1 <- logConDiscrMLE(y,output = FALSE)
  phat <- exp(mle1$psi)/sum(exp(mle1$psi))
  
  mu_low <- sum(pn*grid[-length(grid)])
  mu_high <- sum(pn*grid[-1])
  mu_mid <-(mu_low+mu_high)/2
  
  sigma = sd(y)
  
  mu_n <- EMemp(x, grid, start=2,sigma=sigma)$mu_hat
  
  Delta <- mu_n-mu_low
  alpha <- alpha
  beta <- 2*alpha*(Delta/delta-0.5)
  
  #  TODO   
  #  alpha+beta
  #  alpha-beta
  
  nn<-B
  
  old_seed <- .Random.seed
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  set.seed(seednum)
  
  ystar <- sample(mle1$x * delta, nn, phat, replace=TRUE) 
  zstar <- delta*rbeta(nn, alpha+beta, alpha-beta)
  
  xstar <- ystar+zstar
  
  mle2 <- logConDens(xstar,smoothed=TRUE,print = FALSE)
  
  return(list(fhatn = mle2, sumphat=sum(phat), checkZ = min(alpha+beta,alpha-beta)))
  
}


get_fhatn <- function(x, grid, B=10000, alpha=2, seednum=20180621,log_conc = TRUE){
  
  n  <- length(x)
  pn <- hist(x, breaks=grid, plot=FALSE)$counts / n   # length = length(grid)-1
  
  delta <- diff(grid)[1]   # uniform
  g0    <- grid[1]
  
  # Use left edges as discrete support
  grid_left <- grid[-length(grid)]
  y <- rep(grid_left, times = round(n * pn))
  
  # Map to integers on the grid lattice
  y_int <- as.integer(round((y - g0) / delta))
  mle1 <- logConDiscrMLE(y_int, output = FALSE)
  phat <- exp(mle1$psi) / sum(exp(mle1$psi))
  
  # Moments computed in ORIGINAL scale
  mu_low  <- sum(pn * grid_left)
  mu_high <- sum(pn * grid[-1])
  mu_mid  <- (mu_low + mu_high) / 2
  
  sigma <- sd(y)
  
  #mids  <- (grid[-length(grid)] + grid[-1]) / 2
  #mu    <- sum(pn * mids)
  #sigma <- sqrt(sum(pn * (mids - mu)^2) + delta^2/12)
  
  mu_n <- EMemp(x, grid, start=2, sigma=sigma)$mu_hat
  
  Delta <- mu_n - mu_low
  beta  <- 2 * alpha * (Delta/delta - 0.5)
  
  nn <- B
  old_seed <- .Random.seed
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  set.seed(seednum)
  
  # Sample discrete part on integer support, then map BACK with g0 + k*delta
  ystar <- sample(mle1$x, nn, prob = phat, replace = TRUE)
  ystar <- g0 + ystar * delta
  
  zstar <- delta * rbeta(nn, alpha + beta, alpha - beta)
  xstar <- ystar + zstar
  
  if (log_conc == TRUE) {
    mle2 <- logConDens(xstar, smoothed=TRUE, print=FALSE)
  }else{
    mle2 <- logConDens(xstar, smoothed=FALSE, print=FALSE) 
  }
  
  list(fhatn = mle2, sumphat = sum(phat), checkZ = min(alpha+beta, alpha-beta))
  
}

get_fhatn_with_n <- function(n_i, grid, B=10000, alpha=2, seednum=20180621){
  
  delta = max(diff(grid)) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
  # (Assuming a uniform grid, delta remains the same for any i)
  
  n = sum(n_i)
  pn <- n_i/n
  
  delta <- diff(grid)[1]   # uniform
  g0    <- grid[1]
  
  # Use left edges as discrete support
  grid_left <- grid[-length(grid)]
  y <- rep(grid_left, times = round(n * pn))
  
  # Map to integers on the grid lattice
  y_int <- as.integer(round((y - g0) / delta))
  
  mle1 <- logConDiscrMLE(y_int, output = FALSE)
  phat <- exp(mle1$psi) / sum(exp(mle1$psi))
  
  # Moments computed in ORIGINAL scale
  mu_low  <- sum(pn * grid_left)
  mu_high <- sum(pn * grid[-1])
  mu_mid  <- (mu_low + mu_high) / 2
  
  sigma <- sd(y)
  
  mu_n <- EMemp_with_n(n_i,grid)$mu_hat
  
  Delta <- mu_n - mu_low
  beta  <- 2 * alpha * (Delta/delta - 0.5)
  
  nn <- B
  old_seed <- .Random.seed
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  set.seed(seednum)
  
  # Sample discrete part on integer support, then map BACK with g0 + k*delta
  ystar <- sample(mle1$x, nn, prob = phat, replace = TRUE)
  ystar <- g0 + ystar * delta
  
  zstar <- delta * rbeta(nn, alpha + beta, alpha - beta)
  xstar <- ystar + zstar
  
  mle2 <- logConDens(xstar, smoothed=TRUE, print=FALSE)
  
  list(fhatn = mle2, sumphat = sum(phat), checkZ = min(alpha+beta, alpha-beta))
  
}


