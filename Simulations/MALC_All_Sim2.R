library(jmuOutlier)
library(EnvStats)

# Simulation Study 

# Range 

range_laplace = c(-Inf,Inf)
range_chisq = c(0,Inf)
range_lnorm = c(0,Inf)
range_wei = c(0,Inf)
range_pareto = c(0,Inf)

# Pdfs 

pdf_laplace <- function(x) {
  dlaplace(x,mu,sd = sigma)
}

pdf_chisq <- function(x) {
  dchisq(x,df)
}

pdf_lnorm <- function(x) {
  dlnorm(x,meanlog = mul,sdlog = sigma)
}

pdf_weibull <- function(x) {
  dweibull(x,shape)
}

pdf_pareto <- function(x) {
  dpareto(x,location = location,shape = shape_p)
}

n <- 1000
B = 100

res_pdf_laplace_L1 = rep(0,B)
res_pdf_chisq_L1 = rep(0,B)
res_pdf_lnorm_L1 = rep(0,B)
res_pdf_wei_L1 = rep(0,B)
res_pdf_pareto_L1 = rep(0,B)

res_pdf_laplace_L2 = rep(0,B)
res_pdf_chisq_L2 = rep(0,B)
res_pdf_lnorm_L2 = rep(0,B)
res_pdf_wei_L2 = rep(0,B)
res_pdf_pareto_L2 = rep(0,B)

# Parameters

alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
mul = 1;
df=3;
shape_p = 4;

# Standardized Delta Values # 

laplace_deltas  = c(0.7071, 1.4142, 2.1213, 2.8284, 3.5355)
chisqr_deltas  = c(1.2247, 2.4495, 3.6742, 4.8990, 6.1237)
lnorm_deltas = c(2.9374, 5.8747, 8.8121, 11.7495,14.6869)
weibull_deltas = c(0.2316, 0.4633, 0.6949, 0.9265, 1.1581)
pareto_deltas = c(0.4714, 0.9428, 1.4142, 1.8856, 2.3570)

# Corresponding k values to Standardized Delta Values # 

laplace_ks  = c(36, 18, 12, 9, 7)
chisqr_ks  = c(22, 11, 7, 5, 4)
lnorm_ks = c(40, 20, 13, 10, 8)
weibull_ks = c(38, 19, 12, 9,  7)
pareto_ks = c(50, 25, 16, 12, 10)

# Simulations Within the Distributions 

colors1 = rep("#FDB462",5)

# Colors for Distribution

colors2 = rep("#B3DE69",5)

# Colors for Beta Distribution

colors3 = rep("#BC80BD",5)

# Colors for Bernoulli Distribution

colors4 = rep("cyan", 5)

# Colors for Poisson Distribution

colors5 = rep("#FB8072",5)

# Laplace Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

MM = qlaplace(0.9999999,mu,sigma) + 3
laplace_grid <- seq(-MM,MM, by = laplace_deltas[1])

MM1 <- ceiling(qlaplace(0.9999, mu,sigma)) * 1
MM2 <- ceiling(qlaplace(0.9999, mu,sigma)) * 1.5
MM3 <- ceiling(qlaplace(0.9999, mu,sigma)) * 2
MM4 <- ceiling(qlaplace(0.9999, mu,sigma)) * 2.5
MM5 <- ceiling(qlaplace(0.9999, mu,sigma)) * 3 

laplace_grid1 <- seq(-MM1,MM1, by = laplace_deltas[1])
laplace_grid2 <- seq(-MM2,MM2, by = laplace_deltas[1])
laplace_grid3 <- seq(-MM3,MM3, by = laplace_deltas[1])
laplace_grid4 <- seq(-MM4,MM4, by = laplace_deltas[1])
laplace_grid5 <- seq(-MM5,MM5, by = laplace_deltas[1])

make_breaks_cover <- function(x, grid) {
  # grid is your proposed breaks (sorted, equally spaced)
  step <- grid[2] - grid[1]
  lo <- min(grid[1],  floor(min(x)/step) * step)
  hi <- max(tail(grid,1), ceiling(max(x)/step) * step)
  seq(lo-step, hi+step, by = step)
}

# Method # 

res_lap_pdf_n_11_method = rep(0,B)
res_lap_pdf_n_21_method = rep(0,B)
res_lap_pdf_n_31_method = rep(0,B)
res_lap_pdf_n_41_method = rep(0,B)
res_lap_pdf_n_51_method = rep(0,B)

res_lap_pdf_n_12_method = rep(0,B)
smoothed_res_lap_pdf_n_12_method = rep(0,B)
res_lap_pdf_n_22_method = rep(0,B)
smoothed_res_lap_pdf_n_22_method = rep(0,B)
res_lap_pdf_n_32_method = rep(0,B)
smoothed_res_lap_pdf_n_32_method = rep(0,B)
res_lap_pdf_n_42_method = rep(0,B)
smoothed_res_lap_pdf_n_42_method = rep(0,B)
res_lap_pdf_n_52_method = rep(0,B)
smoothed_res_lap_pdf_n_52_method = rep(0,B)

res_lap_pdf_n_12_method2 = rep(0,B)
res_lap_pdf_n_22_method2 = rep(0,B)
res_lap_pdf_n_32_method2 = rep(0,B)
res_lap_pdf_n_42_method2 = rep(0,B)
res_lap_pdf_n_52_method2 = rep(0,B)

res_lap_pdf_n_12_method3 = rep(0,B)
res_lap_pdf_n_22_method3 = rep(0,B)
res_lap_pdf_n_32_method3 = rep(0,B)
res_lap_pdf_n_42_method3 = rep(0,B)
res_lap_pdf_n_52_method3 = rep(0,B)

res_lap_pdf_n_12_method4 = rep(0,B)
res_lap_pdf_n_22_method4 = rep(0,B)
res_lap_pdf_n_32_method4 = rep(0,B)
res_lap_pdf_n_42_method4 = rep(0,B)
res_lap_pdf_n_52_method4 = rep(0,B)

set.seed(11)

for (i in 1:B) {
  
  # Laplace Distribution - n1
  
  x <- (rlaplace(n[1], mean = mu,sd=sigma))
  laplace_grid1 <- make_breaks_cover(x, laplace_grid1)
  
  # res_lap_pdf_n_11_method[i] = L1_Distance_calc(x,laplace_grid,pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grid1,pdf_laplace)
  res_lap_pdf_n_12_method[i] = err_malc[1]
  smoothed_res_lap_pdf_n_12_method[i] = err_malc[2]
  res_lap_pdf_n_12_method2[i] = L2_from_res(x,laplace_grid1,pdf_laplace)
  res_lap_pdf_n_12_method3[i] = L2_binnednp(x,laplace_grid1,pdf_laplace)
  res_lap_pdf_n_12_method4[i] = L2_from_kernsmooth(x,laplace_grid1,pdf_laplace)
  
  # Laplace Distribution - n2
  
  x <- (rlaplace(n[2], mean = mu,sd=sigma))
  laplace_grid2 <- make_breaks_cover(x, laplace_grid2)
  
  # res_lap_pdf_n_21_method[i] = L1_Distance_calc(x,laplace_grid,pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grid2,pdf_laplace)
  res_lap_pdf_n_22_method[i] = err_malc[1]
  smoothed_res_lap_pdf_n_22_method[i] = err_malc[2]
  res_lap_pdf_n_22_method2[i] = L2_from_res(x,laplace_grid2,pdf_laplace)
  res_lap_pdf_n_22_method3[i] = L2_binnednp(x,laplace_grid2,pdf_laplace)
  res_lap_pdf_n_22_method4[i] = L2_from_kernsmooth(x,laplace_grid2,pdf_laplace)
  
  # Laplace Distribution - n3
  
  x <- (rlaplace(n[3], mean = mu,sd=sigma))
  laplace_grid3 <- make_breaks_cover(x, laplace_grid3)
  
  # res_lap_pdf_n_31_method[i] = L1_Distance_calc(x,laplace_grid,pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grid3,pdf_laplace)
  res_lap_pdf_n_32_method[i] = err_malc[1]
  smoothed_res_lap_pdf_n_32_method[i] = err_malc[2]
  res_lap_pdf_n_32_method2[i] = L2_from_res(x,laplace_grid3,pdf_laplace)
  res_lap_pdf_n_32_method3[i] = L2_binnednp(x,laplace_grid3,pdf_laplace)
  res_lap_pdf_n_32_method4[i] = L2_from_kernsmooth(x,laplace_grid3,pdf_laplace)
  
  # Laplace Distribution - n4
  
  x <- (rlaplace(n[4], mean = mu,sd=sigma))
  laplace_grid4 <- make_breaks_cover(x, laplace_grid4)
  
  # res_lap_pdf_n_41_method[i] = L1_Distance_calc(x,laplace_grid,pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grid4,pdf_laplace)
  res_lap_pdf_n_42_method[i] = err_malc[1]
  smoothed_res_lap_pdf_n_42_method[i] = err_malc[2]
  res_lap_pdf_n_42_method2[i] = L2_from_res(x,laplace_grid4,pdf_laplace)
  res_lap_pdf_n_42_method3[i] = L2_binnednp(x,laplace_grid4,pdf_laplace)
  res_lap_pdf_n_42_method4[i] = L2_from_kernsmooth(x,laplace_grid4,pdf_laplace)
  
  # Laplace Distribution - n5
  
  x <- (rlaplace(n[5], mean = mu,sd=sigma))
  laplace_grid5 <- make_breaks_cover(x, laplace_grid5)
  
  # res_lap_pdf_n_51_method[i] = L1_Distance_calc(x,laplace_grid,pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grid5,pdf_laplace)
  res_lap_pdf_n_52_method[i] = err_malc[1]
  smoothed_res_lap_pdf_n_52_method[i] = err_malc[2]
  res_lap_pdf_n_52_method2[i] = L2_from_res(x,laplace_grid5,pdf_laplace)
  res_lap_pdf_n_52_method3[i] = L2_binnednp(x,laplace_grid5,pdf_laplace)
  res_lap_pdf_n_52_method4[i] = L2_from_kernsmooth(x,laplace_grid5,pdf_laplace)
  
  print(i)
  
}

# Method # 

boxplot(res_lap_pdf_n_11_method,res_lap_pdf_n_21_method,res_lap_pdf_n_31_method,res_lap_pdf_n_41_method,res_lap_pdf_n_51_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors1,ylim=c(0,1.2))

boxplot(res_lap_pdf_n_12_method,res_lap_pdf_n_22_method,res_lap_pdf_n_32_method,res_lap_pdf_n_42_method,res_lap_pdf_n_52_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors1,ylim=c(0,0.4))

boxplot(res_lap_pdf_n_12_method2,res_lap_pdf_n_22_method2,res_lap_pdf_n_32_method2,res_lap_pdf_n_42_method2,res_lap_pdf_n_52_method2,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors1,ylim=c(0,0.4))

boxplot(res_lap_pdf_n_12_method3,res_lap_pdf_n_22_method3,res_lap_pdf_n_32_method3,res_lap_pdf_n_42_method3,res_lap_pdf_n_52_method3,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors1,ylim=c(0,0.4))

boxplot(res_lap_pdf_n_12_method4,res_lap_pdf_n_22_method4,res_lap_pdf_n_32_method4,res_lap_pdf_n_42_method4,res_lap_pdf_n_52_method4,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors1,ylim=c(0,0.4))


# Chi-Square Distribution - Sample Size Simulation 

MM = qchisq(0.999999999999,df) 
chisq_grid <- seq(0,MM, by = chisqr_deltas[1])

MM1 <- ceiling(qchisq(0.9999, df)) * 1
MM2 <- ceiling(qchisq(0.9999, df)) * 1.5
MM3 <- ceiling(qchisq(0.9999, df)) * 2
MM4 <- ceiling(qchisq(0.9999, df)) * 2.5
MM5 <- ceiling(qchisq(0.9999, df)) * 3 

chisq_grid1 <- seq(0,MM1, by = chisqr_deltas[1])
chisq_grid2 <- seq(0,MM2, by = chisqr_deltas[1])
chisq_grid3 <- seq(0,MM3, by = chisqr_deltas[1])
chisq_grid4 <- seq(0,MM4, by = chisqr_deltas[1])
chisq_grid5 <- seq(0,MM5, by = chisqr_deltas[1])

make_breaks_cover_0 <- function(x, grid,location = 0) {
  # grid is your proposed breaks (sorted, equally spaced)
  step <- grid[2] - grid[1]
  hi <- max(tail(grid,1), ceiling(max(x)/step) * step)
  seq(location, hi+step, by = step)
}

# Method # 

res_chisq_pdf_n_11_method = rep(0,B)
res_chisq_pdf_n_21_method = rep(0,B)
res_chisq_pdf_n_31_method = rep(0,B)
res_chisq_pdf_n_41_method = rep(0,B)
res_chisq_pdf_n_51_method = rep(0,B)

res_chisq_pdf_n_12_method = rep(0,B)
smoothed_res_chisq_pdf_n_12_method = rep(0,B)
res_chisq_pdf_n_22_method = rep(0,B)
smoothed_res_chisq_pdf_n_22_method = rep(0,B)
res_chisq_pdf_n_32_method = rep(0,B)
smoothed_res_chisq_pdf_n_32_method = rep(0,B)
res_chisq_pdf_n_42_method = rep(0,B)
smoothed_res_chisq_pdf_n_42_method = rep(0,B)
res_chisq_pdf_n_52_method = rep(0,B)
smoothed_res_chisq_pdf_n_52_method = rep(0,B)

res_chisq_pdf_n_12_method2 = rep(0,B)
res_chisq_pdf_n_22_method2 = rep(0,B)
res_chisq_pdf_n_32_method2 = rep(0,B)
res_chisq_pdf_n_42_method2 = rep(0,B)
res_chisq_pdf_n_52_method2 = rep(0,B)

res_chisq_pdf_n_12_method3 = rep(0,B)
res_chisq_pdf_n_22_method3 = rep(0,B)
res_chisq_pdf_n_32_method3 = rep(0,B)
res_chisq_pdf_n_42_method3 = rep(0,B)
res_chisq_pdf_n_52_method3 = rep(0,B)

res_chisq_pdf_n_12_method4 = rep(0,B)
res_chisq_pdf_n_22_method4 = rep(0,B)
res_chisq_pdf_n_32_method4 = rep(0,B)
res_chisq_pdf_n_42_method4 = rep(0,B)
res_chisq_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Chi-square Distribution - n1
  x <- sort(rchisq(n[1], df))
  chisq_grid1 <- make_breaks_cover_0(x, chisq_grid1)
  
  err_malc = L2_Distance_calc_both(x,chisq_grid1,pdf_chisq)
  res_chisq_pdf_n_12_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_n_12_method[i] = err_malc[2]
  res_chisq_pdf_n_12_method2[i] <- L2_from_res(x,chisq_grid1,pdf_chisq)
  res_chisq_pdf_n_12_method3[i] <- L2_binnednp(x,chisq_grid1,pdf_chisq)
  res_chisq_pdf_n_12_method4[i] <- L2_from_kernsmooth(x,chisq_grid1,pdf_chisq)
  
  # n2
  x <- sort(rchisq(n[2], df))
  chisq_grid2 <- make_breaks_cover_0(x, chisq_grid2)
  
  err_malc = L2_Distance_calc_both(x,chisq_grid2,pdf_chisq)
  res_chisq_pdf_n_22_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_n_22_method[i] = err_malc[2]
  res_chisq_pdf_n_22_method2[i] <- L2_from_res(x,chisq_grid2,pdf_chisq)
  res_chisq_pdf_n_22_method3[i] <- L2_binnednp(x,chisq_grid2,pdf_chisq)
  res_chisq_pdf_n_22_method4[i] <- L2_from_kernsmooth(x,chisq_grid2,pdf_chisq)
  
  # n3
  x <- sort(rchisq(n[3], df))
  chisq_grid3 <- make_breaks_cover_0(x, chisq_grid3)
  
  err_malc = L2_Distance_calc_both(x,chisq_grid3,pdf_chisq)
  res_chisq_pdf_n_32_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_n_32_method[i] = err_malc[2]
  res_chisq_pdf_n_32_method2[i] <- L2_from_res(x,chisq_grid3,pdf_chisq)
  res_chisq_pdf_n_32_method3[i] <- L2_binnednp(x,chisq_grid3,pdf_chisq)
  res_chisq_pdf_n_32_method4[i] <- L2_from_kernsmooth(x,chisq_grid3,pdf_chisq)
  
  # n4
  x <- sort(rchisq(n[4], df))
  chisq_grid4 <- make_breaks_cover_0(x, chisq_grid4)
  
  err_malc = L2_Distance_calc_both(x,chisq_grid4,pdf_chisq)
  res_chisq_pdf_n_42_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_n_42_method[i] = err_malc[2]
  res_chisq_pdf_n_42_method2[i] <- L2_from_res(x,chisq_grid4,pdf_chisq)
  res_chisq_pdf_n_42_method3[i] <- L2_binnednp(x,chisq_grid4,pdf_chisq)
  res_chisq_pdf_n_42_method4[i] <- L2_from_kernsmooth(x,chisq_grid4,pdf_chisq)
  
  # n5
  x <- sort(rchisq(n[5], df))
  chisq_grid5 <- make_breaks_cover_0(x, chisq_grid5)
  
  err_malc = L2_Distance_calc_both(x,chisq_grid5,pdf_chisq)
  res_chisq_pdf_n_52_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_n_52_method[i] = err_malc[2]
  res_chisq_pdf_n_52_method2[i] <- L2_from_res(x,chisq_grid5,pdf_chisq)
  res_chisq_pdf_n_52_method3[i] <- L2_binnednp(x,chisq_grid5,pdf_chisq)
  res_chisq_pdf_n_52_method4[i] <- L2_from_kernsmooth(x,chisq_grid5,pdf_chisq)
  
  print(i)
  
}

# Method # 

boxplot(res_chisq_pdf_n_12_method,res_chisq_pdf_n_22_method,res_chisq_pdf_n_32_method,res_chisq_pdf_n_42_method,res_chisq_pdf_n_52_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_n_12_method2,res_chisq_pdf_n_22_method2,res_chisq_pdf_n_32_method2,res_chisq_pdf_n_42_method2,res_chisq_pdf_n_52_method2,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_n_12_method3,res_chisq_pdf_n_22_method3,res_chisq_pdf_n_32_method3,res_chisq_pdf_n_42_method3,res_chisq_pdf_n_52_method3,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_n_12_method4,res_chisq_pdf_n_22_method4,res_chisq_pdf_n_32_method4,res_chisq_pdf_n_42_method4,res_chisq_pdf_n_52_method4,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2,ylim=c(0,0.3))


# Log-Normal Distribution - Sample Size Simulation 

MM = qlnorm(0.99999999999999,mu,sigma) + 30
lnrom_grid <- seq(0,MM, by = lnorm_deltas[1])

MM1 <- ceiling(qlnorm(0.9999,mu,sigma)) * 1
MM2 <- ceiling(qlnorm(0.9999,mu,sigma)) * 1.5
MM3 <- ceiling(qlnorm(0.9999,mu,sigma)) * 2
MM4 <- ceiling(qlnorm(0.9999,mu,sigma)) * 2.5
MM5 <- ceiling(qlnorm(0.9999,mu,sigma)) * 3 

lnrom_grid1 <- seq(0,MM1, by = lnorm_deltas[1])
lnrom_grid2 <- seq(0,MM2, by = lnorm_deltas[1])
lnrom_grid3 <- seq(0,MM3, by = lnorm_deltas[1])
lnrom_grid4 <- seq(0,MM4, by = lnorm_deltas[1])
lnrom_grid5 <- seq(0,MM5, by = lnorm_deltas[1])

# Method # 

res_logn_pdf_n_11_method = rep(0,B)
res_logn_pdf_n_21_method = rep(0,B)
res_logn_pdf_n_31_method = rep(0,B)
res_logn_pdf_n_41_method = rep(0,B)
res_logn_pdf_n_51_method = rep(0,B)

res_logn_pdf_n_12_method = rep(0,B)
smoothed_res_logn_pdf_n_12_method = rep(0,B)
res_logn_pdf_n_22_method = rep(0,B)
smoothed_res_logn_pdf_n_22_method = rep(0,B)
res_logn_pdf_n_32_method = rep(0,B)
smoothed_res_logn_pdf_n_32_method = rep(0,B)
res_logn_pdf_n_42_method = rep(0,B)
smoothed_res_logn_pdf_n_42_method = rep(0,B)
res_logn_pdf_n_52_method = rep(0,B)
smoothed_res_logn_pdf_n_52_method = rep(0,B)

res_logn_pdf_n_12_method2 = rep(0,B)
res_logn_pdf_n_22_method2 = rep(0,B)
res_logn_pdf_n_32_method2 = rep(0,B)
res_logn_pdf_n_42_method2 = rep(0,B)
res_logn_pdf_n_52_method2 = rep(0,B)

res_logn_pdf_n_12_method3 = rep(0,B)
res_logn_pdf_n_22_method3 = rep(0,B)
res_logn_pdf_n_32_method3 = rep(0,B)
res_logn_pdf_n_42_method3 = rep(0,B)
res_logn_pdf_n_52_method3 = rep(0,B)

res_logn_pdf_n_12_method4 = rep(0,B)
res_logn_pdf_n_22_method4 = rep(0,B)
res_logn_pdf_n_32_method4 = rep(0,B)
res_logn_pdf_n_42_method4 = rep(0,B)
res_logn_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Log-Normal Distribution - n1
  x <- sort(rlnorm(n[1], mul, sigma))
  lnrom_grid1 <- make_breaks_cover(x, lnrom_grid1)
  
  err_malc = L2_Distance_calc_both(x,lnrom_grid1,pdf_lnorm)
  res_logn_pdf_n_12_method[i] = err_malc[1]
  smoothed_res_logn_pdf_n_12_method[i] = err_malc[2]
  res_logn_pdf_n_12_method2[i] <- L2_from_res(x,lnrom_grid1,pdf_lnorm)
  res_logn_pdf_n_12_method3[i] <- L2_binnednp(x,lnrom_grid1,pdf_lnorm)
  res_logn_pdf_n_12_method4[i] <- L2_from_kernsmooth(x,lnrom_grid1,pdf_lnorm)
  
  # n2
  x <- sort(rlnorm(n[2], mul, sigma))
  lnrom_grid2 <- make_breaks_cover(x, lnrom_grid2)
  
  err_malc = L2_Distance_calc_both(x,lnrom_grid2,pdf_lnorm)
  res_logn_pdf_n_22_method[i] = err_malc[1]
  smoothed_res_logn_pdf_n_22_method[i] = err_malc[2]
  res_logn_pdf_n_22_method2[i] <- L2_from_res(x,lnrom_grid2,pdf_lnorm)
  res_logn_pdf_n_22_method3[i] <- L2_binnednp(x,lnrom_grid2,pdf_lnorm)
  res_logn_pdf_n_22_method4[i] <- L2_from_kernsmooth(x,lnrom_grid2,pdf_lnorm)
  
  # n3
  x <- sort(rlnorm(n[3], mul, sigma))
  lnrom_grid3 <- make_breaks_cover(x, lnrom_grid3)
  
  err_malc = L2_Distance_calc_both(x,lnrom_grid3,pdf_lnorm)
  res_logn_pdf_n_32_method[i] = err_malc[1]
  smoothed_res_logn_pdf_n_32_method[i] = err_malc[2]
  res_logn_pdf_n_32_method2[i] <- L2_from_res(x,lnrom_grid3,pdf_lnorm)
  res_logn_pdf_n_32_method3[i] <- L2_binnednp(x,lnrom_grid3,pdf_lnorm)
  res_logn_pdf_n_32_method4[i] <- L2_from_kernsmooth(x,lnrom_grid3,pdf_lnorm)
  
  # n4
  x <- sort(rlnorm(n[4], mul, sigma))
  lnrom_grid4 <- make_breaks_cover(x, lnrom_grid4)
  
  err_malc = L2_Distance_calc_both(x,lnrom_grid4,pdf_lnorm)
  res_logn_pdf_n_42_method[i] = err_malc[1]
  smoothed_res_logn_pdf_n_42_method[i] = err_malc[2]
  res_logn_pdf_n_42_method2[i] <- L2_from_res(x,lnrom_grid4,pdf_lnorm)
  res_logn_pdf_n_42_method3[i] <- L2_binnednp(x,lnrom_grid4,pdf_lnorm)
  res_logn_pdf_n_42_method4[i] <- L2_from_kernsmooth(x,lnrom_grid4,pdf_lnorm)
  
  # n5
  x <- sort(rlnorm(n[5], mul, sigma))
  lnrom_grid5 <- make_breaks_cover(x, lnrom_grid5)
  
  err_malc = L2_Distance_calc_both(x,lnrom_grid5,pdf_lnorm)
  res_logn_pdf_n_52_method[i] = err_malc[1]
  smoothed_res_logn_pdf_n_52_method[i] = err_malc[2]
  res_logn_pdf_n_52_method2[i] <- L2_from_res(x,lnrom_grid5,pdf_lnorm)
  res_logn_pdf_n_52_method3[i] <- L2_binnednp(x,lnrom_grid5,pdf_lnorm)
  res_logn_pdf_n_52_method4[i] <- L2_from_kernsmooth(x,lnrom_grid5,pdf_lnorm)
  
  print(i)
  
}

# Method # 

boxplot(res_logn_pdf_n_11_method,res_logn_pdf_n_21_method,res_logn_pdf_n_31_method,res_logn_pdf_n_41_method,res_logn_pdf_n_51_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3,ylim=c(0,1.6))

boxplot(res_logn_pdf_n_12_method,res_logn_pdf_n_22_method,res_logn_pdf_n_32_method,res_logn_pdf_n_42_method,res_logn_pdf_n_52_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3,ylim=c(0,0.25))

boxplot(res_logn_pdf_n_12_method2,res_logn_pdf_n_22_method2,res_logn_pdf_n_32_method2,res_logn_pdf_n_42_method2,res_logn_pdf_n_52_method2,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3,ylim=c(0,0.25))

boxplot(res_logn_pdf_n_12_method3,res_logn_pdf_n_22_method3,res_logn_pdf_n_32_method3,res_logn_pdf_n_42_method3,res_logn_pdf_n_52_method3,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3,ylim=c(0,0.25))

boxplot(res_logn_pdf_n_12_method4,res_logn_pdf_n_22_method4,res_logn_pdf_n_32_method4,res_logn_pdf_n_42_method4,res_logn_pdf_n_52_method4,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3,ylim=c(0,0.25))

# Weibull Distribution - Sample Size Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grid <- seq(0,MM, by = weibull_deltas[1])

# Method # 

res_wei_pdf_n_11_method = rep(0,B)
res_wei_pdf_n_21_method = rep(0,B)
res_wei_pdf_n_31_method = rep(0,B)
res_wei_pdf_n_41_method = rep(0,B)
res_wei_pdf_n_51_method = rep(0,B)

res_wei_pdf_n_12_method = rep(0,B)
smoothed_res_wei_pdf_n_12_method = rep(0,B)
res_wei_pdf_n_22_method = rep(0,B)
smoothed_res_wei_pdf_n_22_method = rep(0,B)
res_wei_pdf_n_32_method = rep(0,B)
smoothed_res_wei_pdf_n_32_method = rep(0,B)
res_wei_pdf_n_42_method = rep(0,B)
smoothed_res_wei_pdf_n_42_method = rep(0,B)
res_wei_pdf_n_52_method = rep(0,B)
smoothed_res_wei_pdf_n_52_method = rep(0,B)

res_wei_pdf_n_12_method2 = rep(0,B)
res_wei_pdf_n_22_method2 = rep(0,B)
res_wei_pdf_n_32_method2 = rep(0,B)
res_wei_pdf_n_42_method2 = rep(0,B)
res_wei_pdf_n_52_method2 = rep(0,B)

res_wei_pdf_n_12_method3 = rep(0,B)
res_wei_pdf_n_22_method3 = rep(0,B)
res_wei_pdf_n_32_method3 = rep(0,B)
res_wei_pdf_n_42_method3 = rep(0,B)
res_wei_pdf_n_52_method3 = rep(0,B)

res_wei_pdf_n_12_method4 = rep(0,B)
res_wei_pdf_n_22_method4 = rep(0,B)
res_wei_pdf_n_32_method4 = rep(0,B)
res_wei_pdf_n_42_method4 = rep(0,B)
res_wei_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Weibull Distribution - n1
  
  x = (rweibull(n[1],shape))
  
  # res_wei_pdf_n_11_method[i] = L1_Distance_calc(x,weibull_grid,pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_12_method[i] = err_malc[1]
  smoothed_res_wei_pdf_n_12_method[i] = err_malc[2]
  res_wei_pdf_n_12_method2[i] = L2_from_res(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_12_method3[i] = L2_binnednp(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_12_method4[i] = L2_from_kernsmooth(x,weibull_grid,pdf_weibull)
  
  # Weibull Distribution - n2
  
  x = (rweibull(n[2],shape))
  
  # res_wei_pdf_n_21_method[i] = L1_Distance_calc(x,weibull_grid,pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_22_method[i] = err_malc[1]
  smoothed_res_wei_pdf_n_22_method[i] = err_malc[2]
  res_wei_pdf_n_22_method2[i] = L2_from_res(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_22_method3[i] = L2_binnednp(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_22_method4[i] = L2_from_kernsmooth(x,weibull_grid,pdf_weibull)
  
  # Weibull Distribution - n3
  
  x = (rweibull(n[3],shape))
  
  # res_wei_pdf_n_31_method[i] = L1_Distance_calc(x,weibull_grid,pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_32_method[i] = err_malc[1]
  smoothed_res_wei_pdf_n_32_method[i] = err_malc[2]
  res_wei_pdf_n_32_method2[i] = L2_from_res(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_32_method3[i] = L2_binnednp(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_32_method4[i] = L2_from_kernsmooth(x,weibull_grid,pdf_weibull)
  
  # Weibull Distribution - n4
  
  x = (rweibull(n[4],shape))
  
  # res_wei_pdf_n_41_method[i] = L1_Distance_calc(x,weibull_grid,pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_42_method[i] = err_malc[1]
  smoothed_res_wei_pdf_n_42_method[i] = err_malc[2]
  res_wei_pdf_n_42_method2[i] = L2_from_res(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_42_method3[i] = L2_binnednp(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_42_method4[i] = L2_from_kernsmooth(x,weibull_grid,pdf_weibull)
  
  # Weibull Distribution - n5
  
  x = (rweibull(n[5],shape))
  
  # res_wei_pdf_n_51_method[i] = L1_Distance_calc(x,weibull_grid,pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_52_method[i] = err_malc[1]
  smoothed_res_wei_pdf_n_52_method[i] = err_malc[2]
  res_wei_pdf_n_52_method2[i] = L2_from_res(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_52_method3[i] = L2_binnednp(x,weibull_grid,pdf_weibull)
  res_wei_pdf_n_52_method4[i] = L2_from_kernsmooth(x,weibull_grid,pdf_weibull)
  
  print(i)
  
}

# Method # 

boxplot(res_wei_pdf_n_11_method,res_wei_pdf_n_21_method,res_wei_pdf_n_31_method,res_wei_pdf_n_41_method,res_wei_pdf_n_51_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4,ylim=c(0,1.0))

boxplot(res_wei_pdf_n_12_method,res_wei_pdf_n_22_method,res_wei_pdf_n_32_method,res_wei_pdf_n_42_method,res_wei_pdf_n_52_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_wei_pdf_n_12_method2,res_wei_pdf_n_22_method2,res_wei_pdf_n_32_method2,res_wei_pdf_n_42_method2,res_wei_pdf_n_52_method2,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_wei_pdf_n_12_method3,res_wei_pdf_n_22_method3,res_wei_pdf_n_32_method3,res_wei_pdf_n_42_method3,res_wei_pdf_n_52_method3,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_wei_pdf_n_12_method4,res_wei_pdf_n_22_method4,res_wei_pdf_n_32_method4,res_wei_pdf_n_42_method4,res_wei_pdf_n_52_method4,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4,ylim=c(0,0.2))

# Pareto Distribution - Sample Size Simulation 

MM = qpareto(0.999999999,location,shape_p) + 3
pareto_grid <- seq(location,MM, by = pareto_deltas[1])

MM1 <- ceiling(qpareto(0.9999,location,shape_p)) * 1
MM2 <- ceiling(qpareto(0.9999,location,shape_p)) * 1.5
MM3 <- ceiling(qpareto(0.9999,location,shape_p)) * 2
MM4 <- ceiling(qpareto(0.9999,location,shape_p)) * 2.5
MM5 <- ceiling(qpareto(0.9999,location,shape_p)) * 3 

pareto_grid1 <- seq(location,MM1, by = pareto_deltas[1])
pareto_grid2 <- seq(location,MM2, by = pareto_deltas[1])
pareto_grid3 <- seq(location,MM3, by = pareto_deltas[1])
pareto_grid4 <- seq(location,MM4, by = pareto_deltas[1])
pareto_grid5 <- seq(location,MM5, by = pareto_deltas[1])

# Method # 

res_pareto_pdf_n_11_method = rep(0,B)
res_pareto_pdf_n_21_method = rep(0,B)
res_pareto_pdf_n_31_method = rep(0,B)
res_pareto_pdf_n_41_method = rep(0,B)
res_pareto_pdf_n_51_method = rep(0,B)

res_pareto_pdf_n_12_method = rep(0,B)
smoothed_res_pareto_pdf_n_12_method = rep(0,B)
res_pareto_pdf_n_22_method = rep(0,B)
smoothed_res_pareto_pdf_n_22_method = rep(0,B)
res_pareto_pdf_n_32_method = rep(0,B)
smoothed_res_pareto_pdf_n_32_method = rep(0,B)
res_pareto_pdf_n_42_method = rep(0,B)
smoothed_res_pareto_pdf_n_42_method = rep(0,B)
res_pareto_pdf_n_52_method = rep(0,B)
smoothed_res_pareto_pdf_n_52_method = rep(0,B)

res_pareto_pdf_n_12_method2 = rep(0,B)
res_pareto_pdf_n_22_method2 = rep(0,B)
res_pareto_pdf_n_32_method2 = rep(0,B)
res_pareto_pdf_n_42_method2 = rep(0,B)
res_pareto_pdf_n_52_method2 = rep(0,B)

res_pareto_pdf_n_12_method3 = rep(0,B)
res_pareto_pdf_n_22_method3 = rep(0,B)
res_pareto_pdf_n_32_method3 = rep(0,B)
res_pareto_pdf_n_42_method3 = rep(0,B)
res_pareto_pdf_n_52_method3 = rep(0,B)

res_pareto_pdf_n_12_method4 = rep(0,B)
res_pareto_pdf_n_22_method4 = rep(0,B)
res_pareto_pdf_n_32_method4 = rep(0,B)
res_pareto_pdf_n_42_method4 = rep(0,B)
res_pareto_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Pareto Distribution - n1
  x <- rpareto(n[1], location, shape_p)
  pareto_grid1 <- make_breaks_cover(x, pareto_grid1)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_n_12_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_n_12_method[i] = err_malc[2]
  res_pareto_pdf_n_12_method2[i] <- L2_from_res(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_n_12_method3[i] <- L2_binnednp(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_n_12_method4[i] <- L2_from_kernsmooth(x,pareto_grid1,pdf_pareto)
  
  # n2
  x <- rpareto(n[2], location, shape_p)
  pareto_grid2 <- make_breaks_cover(x, pareto_grid2)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_n_22_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_n_22_method[i] = err_malc[2]
  res_pareto_pdf_n_22_method2[i] <- L2_from_res(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_n_22_method3[i] <- L2_binnednp(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_n_22_method4[i] <- L2_from_kernsmooth(x,pareto_grid2,pdf_pareto)
  
  # n3
  x <- rpareto(n[3], location, shape_p)
  pareto_grid3 <- make_breaks_cover(x, pareto_grid3)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_n_32_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_n_32_method[i] = err_malc[2]
  res_pareto_pdf_n_32_method2[i] <- L2_from_res(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_n_32_method3[i] <- L2_binnednp(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_n_32_method4[i] <- L2_from_kernsmooth(x,pareto_grid3,pdf_pareto)
  
  # n4
  x <- rpareto(n[4], location, shape_p)
  pareto_grid4 <- make_breaks_cover(x, pareto_grid4)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_n_42_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_n_42_method[i] = err_malc[2]
  res_pareto_pdf_n_42_method2[i] <- L2_from_res(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_n_42_method3[i] <- L2_binnednp(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_n_42_method4[i] <- L2_from_kernsmooth(x,pareto_grid4,pdf_pareto)
  
  # n5
  x <- rpareto(n[5], location, shape_p)
  pareto_grid5 <- make_breaks_cover(x, pareto_grid5)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_n_52_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_n_52_method[i] = err_malc[2]
  res_pareto_pdf_n_52_method2[i] <- L2_from_res(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_n_52_method3[i] <- L2_binnednp(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_n_52_method4[i] <- L2_from_kernsmooth(x,pareto_grid5,pdf_pareto)
  
  print(i)
  
}

# Method #  

boxplot(res_pareto_pdf_n_11_method,res_pareto_pdf_n_21_method,res_pareto_pdf_n_31_method,res_pareto_pdf_n_41_method,res_pareto_pdf_n_51_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5,ylim=c(0,1.6))

boxplot(res_pareto_pdf_n_12_method,res_pareto_pdf_n_22_method,res_pareto_pdf_n_32_method,res_pareto_pdf_n_42_method,res_pareto_pdf_n_52_method,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5,ylim=c(0,1.2))

boxplot(res_pareto_pdf_n_12_method2,res_pareto_pdf_n_22_method2,res_pareto_pdf_n_32_method2,res_pareto_pdf_n_42_method2,res_pareto_pdf_n_52_method2,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5,ylim=c(0,1.2))

boxplot(res_pareto_pdf_n_12_method3,res_pareto_pdf_n_22_method3,res_pareto_pdf_n_32_method3,res_pareto_pdf_n_42_method3,res_pareto_pdf_n_52_method3,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5,ylim=c(0,1.2))

boxplot(res_pareto_pdf_n_12_method4,res_pareto_pdf_n_22_method4,res_pareto_pdf_n_32_method4,res_pareto_pdf_n_42_method4,res_pareto_pdf_n_52_method4,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5,ylim=c(0,1.2))


n = c(10^3)
B = 100

# Laplace Distribution - Grid-Width Simulation 

MM = qlaplace(0.9999999,mu,sigma) + 200
laplace_grids = list(seq(-MM,MM, by = laplace_deltas[5]),
                     seq(-MM,MM, by = laplace_deltas[4]),
                     seq(-MM,MM, by = laplace_deltas[3]),
                     seq(-MM,MM, by = laplace_deltas[2]),
                     seq(-MM,MM, by = laplace_deltas[1]))

# Method # 

res_lap_pdf_gw_11_method = rep(0,B)
res_lap_pdf_gw_21_method = rep(0,B)
res_lap_pdf_gw_31_method = rep(0,B)
res_lap_pdf_gw_41_method = rep(0,B)
res_lap_pdf_gw_51_method = rep(0,B)

res_lap_pdf_gw_12_method = rep(0,B)
smoothed_res_lap_pdf_gw_12_method = rep(0,B)
res_lap_pdf_gw_22_method = rep(0,B)
smoothed_res_lap_pdf_gw_22_method = rep(0,B)
res_lap_pdf_gw_32_method = rep(0,B)
smoothed_res_lap_pdf_gw_32_method = rep(0,B)
res_lap_pdf_gw_42_method = rep(0,B)
smoothed_res_lap_pdf_gw_42_method = rep(0,B)
res_lap_pdf_gw_52_method = rep(0,B)
smoothed_res_lap_pdf_gw_52_method = rep(0,B)

res_lap_pdf_gw_12_method2 = rep(0,B)
res_lap_pdf_gw_22_method2 = rep(0,B)
res_lap_pdf_gw_32_method2 = rep(0,B)
res_lap_pdf_gw_42_method2 = rep(0,B)
res_lap_pdf_gw_52_method2 = rep(0,B)

res_lap_pdf_gw_12_method3 = rep(0,B)
res_lap_pdf_gw_22_method3 = rep(0,B)
res_lap_pdf_gw_32_method3 = rep(0,B)
res_lap_pdf_gw_42_method3 = rep(0,B)
res_lap_pdf_gw_52_method3 = rep(0,B)

res_lap_pdf_gw_12_method4 = rep(0,B)
res_lap_pdf_gw_22_method4 = rep(0,B)
res_lap_pdf_gw_32_method4 = rep(0,B)
res_lap_pdf_gw_42_method4 = rep(0,B)
res_lap_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x <- sort(rlaplace(n, mean = mu,sd=sigma))
  
  # Laplace Distribution - GW1
  
  # res_lap_pdf_gw_11_method[i] = L1_Distance_calc(x,laplace_grids[[1]],pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grids[[1]],pdf_laplace)
  res_lap_pdf_gw_12_method[i] = err_malc[1]
  smoothed_res_lap_pdf_gw_12_method[i] = err_malc[2]
  res_lap_pdf_gw_12_method2[i] = L2_from_res(x,laplace_grids[[1]],pdf_laplace)
  res_lap_pdf_gw_12_method3[i] = L2_binnednp(x,laplace_grids[[1]],pdf_laplace)
  res_lap_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,laplace_grids[[1]],pdf_laplace)
  
  # Laplace Distribution - GW2
  
  # res_lap_pdf_gw_21_method[i] = L1_Distance_calc(x,laplace_grids[[2]],pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grids[[2]],pdf_laplace)
  res_lap_pdf_gw_22_method[i] = err_malc[1]
  smoothed_res_lap_pdf_gw_22_method[i] = err_malc[2]
  res_lap_pdf_gw_22_method2[i] = L2_from_res(x,laplace_grids[[2]],pdf_laplace)
  res_lap_pdf_gw_22_method3[i] = L2_binnednp(x,laplace_grids[[2]],pdf_laplace)
  res_lap_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,laplace_grids[[2]],pdf_laplace)
  
  # Laplace Distribution - GW3
  
  # res_lap_pdf_gw_31_method[i] = L1_Distance_calc(x,laplace_grids[[3]],pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grids[[3]],pdf_laplace)
  res_lap_pdf_gw_32_method[i] = err_malc[1]
  smoothed_res_lap_pdf_gw_32_method[i] = err_malc[2]
  res_lap_pdf_gw_32_method2[i] = L2_from_res(x,laplace_grids[[3]],pdf_laplace)
  res_lap_pdf_gw_32_method3[i] = L2_binnednp(x,laplace_grids[[3]],pdf_laplace)
  res_lap_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,laplace_grids[[3]],pdf_laplace)
  
  # Laplace Distribution - GW4
  
  # res_lap_pdf_gw_41_method[i] = L1_Distance_calc(x,laplace_grids[[4]],pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grids[[4]],pdf_laplace)
  res_lap_pdf_gw_42_method[i] = err_malc[1]
  smoothed_res_lap_pdf_gw_42_method[i] = err_malc[2]
  res_lap_pdf_gw_42_method2[i] = L2_from_res(x,laplace_grids[[4]],pdf_laplace)
  res_lap_pdf_gw_42_method3[i] = L2_binnednp(x,laplace_grids[[4]],pdf_laplace)
  res_lap_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,laplace_grids[[4]],pdf_laplace)
  
  # Laplace Distribution - GW5
  
  # res_lap_pdf_gw_51_method[i] = L1_Distance_calc(x,laplace_grids[[5]],pdf_laplace,range_laplace)
  err_malc = L2_Distance_calc_both(x,laplace_grids[[5]],pdf_laplace)
  res_lap_pdf_gw_52_method[i] = err_malc[1]
  smoothed_res_lap_pdf_gw_52_method[i] = err_malc[2]
  res_lap_pdf_gw_52_method2[i] = L2_from_res(x,laplace_grids[[5]],pdf_laplace)
  res_lap_pdf_gw_52_method3[i] = L2_binnednp(x,laplace_grids[[5]],pdf_laplace)
  res_lap_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,laplace_grids[[5]],pdf_laplace)
  
  print(i)

}

# Method # 

boxplot(res_lap_pdf_gw_11_method,res_lap_pdf_gw_21_method,res_lap_pdf_gw_31_method,res_lap_pdf_gw_41_method,res_lap_pdf_gw_51_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.7))

boxplot(res_lap_pdf_gw_12_method,res_lap_pdf_gw_22_method,res_lap_pdf_gw_32_method,res_lap_pdf_gw_42_method,res_lap_pdf_gw_52_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.8))

boxplot(res_lap_pdf_gw_12_method2,res_lap_pdf_gw_22_method2,res_lap_pdf_gw_32_method2,res_lap_pdf_gw_42_method2,res_lap_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.8))

boxplot(res_lap_pdf_gw_12_method3,res_lap_pdf_gw_22_method3,res_lap_pdf_gw_32_method3,res_lap_pdf_gw_42_method3,res_lap_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.8))

boxplot(res_lap_pdf_gw_12_method4,res_lap_pdf_gw_22_method4,res_lap_pdf_gw_32_method4,res_lap_pdf_gw_42_method4,res_lap_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.8))

# Save to the specified path

# Chi-Square Distribution - Grid Width Simulation 

MM = qchisq(0.9999999,df) + 3
chisqr_grids = list(seq(0,MM, by = chisqr_deltas[5]),
                    seq(0,MM, by = chisqr_deltas[4]),
                    seq(0,MM, by = chisqr_deltas[3]),
                    seq(0,MM, by = chisqr_deltas[2]),
                    seq(0,MM, by = chisqr_deltas[1]))

# Method # 

res_chisq_pdf_gw_11_method = rep(0,B)
res_chisq_pdf_gw_21_method = rep(0,B)
res_chisq_pdf_gw_31_method = rep(0,B)
res_chisq_pdf_gw_41_method = rep(0,B)
res_chisq_pdf_gw_51_method = rep(0,B)

res_chisq_pdf_gw_12_method = rep(0,B)
smoothed_res_chisq_pdf_gw_12_method = rep(0,B)
res_chisq_pdf_gw_22_method = rep(0,B)
smoothed_res_chisq_pdf_gw_22_method = rep(0,B)
res_chisq_pdf_gw_32_method = rep(0,B)
smoothed_res_chisq_pdf_gw_32_method = rep(0,B)
res_chisq_pdf_gw_42_method = rep(0,B)
smoothed_res_chisq_pdf_gw_42_method = rep(0,B)
res_chisq_pdf_gw_52_method = rep(0,B)
smoothed_res_chisq_pdf_gw_52_method = rep(0,B)

res_chisq_pdf_gw_12_method2 = rep(0,B)
res_chisq_pdf_gw_22_method2 = rep(0,B)
res_chisq_pdf_gw_32_method2 = rep(0,B)
res_chisq_pdf_gw_42_method2 = rep(0,B)
res_chisq_pdf_gw_52_method2 = rep(0,B)

res_chisq_pdf_gw_12_method3 = rep(0,B)
res_chisq_pdf_gw_22_method3 = rep(0,B)
res_chisq_pdf_gw_32_method3 = rep(0,B)
res_chisq_pdf_gw_42_method3 = rep(0,B)
res_chisq_pdf_gw_52_method3 = rep(0,B)

res_chisq_pdf_gw_12_method4 = rep(0,B)
res_chisq_pdf_gw_22_method4 = rep(0,B)
res_chisq_pdf_gw_32_method4 = rep(0,B)
res_chisq_pdf_gw_42_method4 = rep(0,B)
res_chisq_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x <- sort(rchisq(n,df))
  
  # Chi-Square Distribution - GW1
  
  # res_chisq_pdf_gw_11_method[i] = L1_Distance_calc(x,chisqr_grids[[1]],pdf_chisq,range_chisq)
  err_malc = L2_Distance_calc_both(x,chisqr_grids[[1]],pdf_chisq)
  res_chisq_pdf_gw_12_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_gw_12_method[i] = err_malc[2]
  res_chisq_pdf_gw_12_method2[i] = L2_from_res(x,chisqr_grids[[1]],pdf_chisq)
  res_chisq_pdf_gw_12_method3[i] = L2_binnednp(x,chisqr_grids[[1]],pdf_chisq)
  res_chisq_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,chisqr_grids[[1]],pdf_chisq)
  
  # Chi-Square Distribution - GW2
  
  # res_chisq_pdf_gw_21_method[i] = L1_Distance_calc(x,chisqr_grids[[2]],pdf_chisq,range_chisq)
  err_malc = L2_Distance_calc_both(x,chisqr_grids[[2]],pdf_chisq)
  res_chisq_pdf_gw_22_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_gw_22_method[i] = err_malc[2]
  res_chisq_pdf_gw_22_method2[i] = L2_from_res(x,chisqr_grids[[2]],pdf_chisq)
  res_chisq_pdf_gw_22_method3[i] = L2_binnednp(x,chisqr_grids[[2]],pdf_chisq)
  res_chisq_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,chisqr_grids[[2]],pdf_chisq)
  
  # Chi-Square Distribution - GW3
  
  # res_chisq_pdf_gw_31_method[i] = L1_Distance_calc(x,chisqr_grids[[3]],pdf_chisq,range_chisq)
  err_malc = L2_Distance_calc_both(x,chisqr_grids[[3]],pdf_chisq)
  res_chisq_pdf_gw_32_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_gw_32_method[i] = err_malc[2]
  res_chisq_pdf_gw_32_method2[i] = L2_from_res(x,chisqr_grids[[3]],pdf_chisq)
  res_chisq_pdf_gw_32_method3[i] = L2_binnednp(x,chisqr_grids[[3]],pdf_chisq)
  res_chisq_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,chisqr_grids[[3]],pdf_chisq)
  
  # Chi-Square Distribution - GW4
  
  # res_chisq_pdf_gw_41_method[i] = L1_Distance_calc(x,chisqr_grids[[4]],pdf_chisq,range_chisq)
  err_malc = L2_Distance_calc_both(x,chisqr_grids[[4]],pdf_chisq)
  res_chisq_pdf_gw_42_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_gw_42_method[i] = err_malc[2]
  res_chisq_pdf_gw_42_method2[i] = L2_from_res(x,chisqr_grids[[4]],pdf_chisq)
  res_chisq_pdf_gw_42_method3[i] = L2_binnednp(x,chisqr_grids[[4]],pdf_chisq)
  res_chisq_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,chisqr_grids[[4]],pdf_chisq)
  
  # Chi-Square Distribution - GW5
  
  # res_chisq_pdf_gw_51_method[i] = L1_Distance_calc(x,chisqr_grids[[5]],pdf_chisq,range_chisq)
  err_malc = L2_Distance_calc_both(x,chisqr_grids[[5]],pdf_chisq)
  res_chisq_pdf_gw_52_method[i] = err_malc[1]
  smoothed_res_chisq_pdf_gw_52_method[i] = err_malc[2]
  res_chisq_pdf_gw_52_method2[i] = L2_from_res(x,chisqr_grids[[5]],pdf_chisq)
  res_chisq_pdf_gw_52_method3[i] = L2_binnednp(x,chisqr_grids[[5]],pdf_chisq)
  res_chisq_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,chisqr_grids[[5]],pdf_chisq)
  
  print(i)
  
}

# Method # 

boxplot(res_chisq_pdf_gw_11_method,res_chisq_pdf_gw_21_method,res_chisq_pdf_gw_31_method,res_chisq_pdf_gw_41_method,res_chisq_pdf_gw_51_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.9))

boxplot(res_chisq_pdf_gw_12_method,res_chisq_pdf_gw_22_method,res_chisq_pdf_gw_32_method,res_chisq_pdf_gw_42_method,res_chisq_pdf_gw_52_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_gw_12_method2,res_chisq_pdf_gw_22_method2,res_chisq_pdf_gw_32_method2,res_chisq_pdf_gw_42_method2,res_chisq_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_gw_12_method3,res_chisq_pdf_gw_22_method3,res_chisq_pdf_gw_32_method3,res_chisq_pdf_gw_42_method3,res_chisq_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.3))

boxplot(res_chisq_pdf_gw_12_method4,res_chisq_pdf_gw_22_method4,res_chisq_pdf_gw_32_method4,res_chisq_pdf_gw_42_method4,res_chisq_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.3))

# Log-Normal Distribution - Grid Width Simulation 

MM = qlnorm(0.99999999999,mu,sigma) + 30
lnorm_grids = list(seq(0,MM, by = lnorm_deltas[5]),
                   seq(0,MM, by = lnorm_deltas[4]),
                   seq(0,MM, by = lnorm_deltas[3]),
                   seq(0,MM, by = lnorm_deltas[2]),
                   seq(0,MM, by = lnorm_deltas[1]))

# Method # 

res_lnorm_pdf_gw_11_method = rep(0,B)
res_lnorm_pdf_gw_21_method = rep(0,B)
res_lnorm_pdf_gw_31_method = rep(0,B)
res_lnorm_pdf_gw_41_method = rep(0,B)
res_lnorm_pdf_gw_51_method = rep(0,B)

res_lnorm_pdf_gw_12_method = rep(0,B)
smoothed_res_lnorm_pdf_gw_12_method = rep(0,B)
res_lnorm_pdf_gw_22_method = rep(0,B)
smoothed_res_lnorm_pdf_gw_22_method = rep(0,B)
res_lnorm_pdf_gw_32_method = rep(0,B)
smoothed_res_lnorm_pdf_gw_32_method = rep(0,B)
res_lnorm_pdf_gw_42_method = rep(0,B)
smoothed_res_lnorm_pdf_gw_42_method = rep(0,B)
res_lnorm_pdf_gw_52_method = rep(0,B)
smoothed_res_lnorm_pdf_gw_52_method = rep(0,B)

res_lnorm_pdf_gw_12_method2 = rep(0,B)
res_lnorm_pdf_gw_22_method2 = rep(0,B)
res_lnorm_pdf_gw_32_method2 = rep(0,B)
res_lnorm_pdf_gw_42_method2 = rep(0,B)
res_lnorm_pdf_gw_52_method2 = rep(0,B)

res_lnorm_pdf_gw_12_method3 = rep(0,B)
res_lnorm_pdf_gw_22_method3 = rep(0,B)
res_lnorm_pdf_gw_32_method3 = rep(0,B)
res_lnorm_pdf_gw_42_method3 = rep(0,B)
res_lnorm_pdf_gw_52_method3 = rep(0,B)

res_lnorm_pdf_gw_12_method4 = rep(0,B)
res_lnorm_pdf_gw_22_method4 = rep(0,B)
res_lnorm_pdf_gw_32_method4 = rep(0,B)
res_lnorm_pdf_gw_42_method4 = rep(0,B)
res_lnorm_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:100) {
  
  x <- sort(rlnorm(n, mul, sigma))
  
  # GW1
  # res_lnorm_pdf_gw_11_method[i] <- L1_Distance_calc(x, lnorm_grids[[1]], pdf_lnorm, range_lnorm)
  
  err_malc = L2_Distance_calc_both(x,lnorm_grids[[1]],pdf_lnorm)
  res_lnorm_pdf_gw_12_method[i] = err_malc[1]
  smoothed_res_lnorm_pdf_gw_12_method[i] = err_malc[2]
  res_lnorm_pdf_gw_12_method2[i] <- L2_from_res(x,lnorm_grids[[1]],pdf_lnorm)
  res_lnorm_pdf_gw_12_method3[i] <- L2_binnednp(x,lnorm_grids[[1]],pdf_lnorm)
  res_lnorm_pdf_gw_12_method4[i] <- L2_from_kernsmooth(x,lnorm_grids[[1]],pdf_lnorm)
  
  # GW2
  # res_lnorm_pdf_gw_21_method[i] <- L1_Distance_calc(x, lnorm_grids[[2]], pdf_lnorm, range_lnorm)
  
  err_malc = L2_Distance_calc_both(x,lnorm_grids[[2]],pdf_lnorm)
  res_lnorm_pdf_gw_22_method[i] = err_malc[1]
  smoothed_res_lnorm_pdf_gw_22_method[i] = err_malc[2]
  res_lnorm_pdf_gw_22_method2[i] <- L2_from_res(x,lnorm_grids[[2]],pdf_lnorm)
  res_lnorm_pdf_gw_22_method3[i] <- L2_binnednp(x,lnorm_grids[[2]],pdf_lnorm)
  res_lnorm_pdf_gw_22_method4[i] <- L2_from_kernsmooth(x,lnorm_grids[[2]],pdf_lnorm)
  
  # GW3
  # res_lnorm_pdf_gw_31_method[i] <- L1_Distance_calc(x, lnorm_grids[[3]], pdf_lnorm, range_lnorm)
  
  err_malc = L2_Distance_calc_both(x,lnorm_grids[[3]],pdf_lnorm)
  res_lnorm_pdf_gw_32_method[i] = err_malc[1]
  smoothed_res_lnorm_pdf_gw_32_method[i] = err_malc[2]
  res_lnorm_pdf_gw_32_method2[i] <- L2_from_res(x,lnorm_grids[[3]],pdf_lnorm)
  res_lnorm_pdf_gw_32_method3[i] <- L2_binnednp(x,lnorm_grids[[3]],pdf_lnorm)
  res_lnorm_pdf_gw_32_method4[i] <- L2_from_kernsmooth(x,lnorm_grids[[3]],pdf_lnorm)
  
  # GW4
  # res_lnorm_pdf_gw_41_method[i] <- L1_Distance_calc(x, lnorm_grids[[4]], pdf_lnorm, range_lnorm)
  
  err_malc = L2_Distance_calc_both(x,lnorm_grids[[4]],pdf_lnorm)
  res_lnorm_pdf_gw_42_method[i] = err_malc[1]
  smoothed_res_lnorm_pdf_gw_42_method[i] = err_malc[2]
  res_lnorm_pdf_gw_42_method2[i] <- L2_from_res(x,lnorm_grids[[4]],pdf_lnorm)
  res_lnorm_pdf_gw_42_method3[i] <- L2_binnednp(x,lnorm_grids[[4]],pdf_lnorm)
  res_lnorm_pdf_gw_42_method4[i] <- L2_from_kernsmooth(x,lnorm_grids[[4]],pdf_lnorm)
  
  # GW5
  # res_lnorm_pdf_gw_51_method[i] <- L1_Distance_calc(x, lnorm_grids[[5]], pdf_lnorm, range_lnorm)
  
  err_malc = L2_Distance_calc_both(x,lnorm_grids[[5]],pdf_lnorm)
  res_lnorm_pdf_gw_52_method[i] = err_malc[1]
  smoothed_res_lnorm_pdf_gw_52_method[i] = err_malc[2]
  res_lnorm_pdf_gw_52_method2[i] <- L2_from_res(x,lnorm_grids[[5]],pdf_lnorm)
  res_lnorm_pdf_gw_52_method3[i] <- L2_binnednp(x,lnorm_grids[[5]],pdf_lnorm)
  res_lnorm_pdf_gw_52_method4[i] <- L2_from_kernsmooth(x,lnorm_grids[[5]],pdf_lnorm)
  
  print(i)
  
}

# Method # 

boxplot(res_lnorm_pdf_gw_11_method,res_lnorm_pdf_gw_21_method,res_lnorm_pdf_gw_31_method,res_lnorm_pdf_gw_41_method,res_lnorm_pdf_gw_51_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,1.5))

boxplot(res_lnorm_pdf_gw_12_method,res_lnorm_pdf_gw_22_method,res_lnorm_pdf_gw_32_method,res_lnorm_pdf_gw_42_method,res_lnorm_pdf_gw_52_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.45))

boxplot(res_lnorm_pdf_gw_12_method2,res_lnorm_pdf_gw_22_method2,res_lnorm_pdf_gw_32_method2,res_lnorm_pdf_gw_42_method2,res_lnorm_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.45))

boxplot(res_lnorm_pdf_gw_12_method3,res_lnorm_pdf_gw_22_method3,res_lnorm_pdf_gw_32_method3,res_lnorm_pdf_gw_42_method3,res_lnorm_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.45))

boxplot(res_lnorm_pdf_gw_12_method4,res_lnorm_pdf_gw_22_method4,res_lnorm_pdf_gw_32_method4,res_lnorm_pdf_gw_42_method4,res_lnorm_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.45))

# Weibull Distribution - Grid Width Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grids = list(seq(0,MM, by = weibull_deltas[5]),
                     seq(0,MM, by = weibull_deltas[4]),
                     seq(0,MM, by = weibull_deltas[3]),
                     seq(0,MM, by = weibull_deltas[2]),
                     seq(0,MM, by = weibull_deltas[1]))

# Method # 

res_wei_pdf_gw_11_method = rep(0,B)
res_wei_pdf_gw_21_method = rep(0,B)
res_wei_pdf_gw_31_method = rep(0,B)
res_wei_pdf_gw_41_method = rep(0,B)
res_wei_pdf_gw_51_method = rep(0,B)

res_wei_pdf_gw_12_method = rep(0,B)
smoothed_res_wei_pdf_gw_12_method = rep(0,B)
res_wei_pdf_gw_22_method = rep(0,B)
smoothed_res_wei_pdf_gw_22_method = rep(0,B)
res_wei_pdf_gw_32_method = rep(0,B)
smoothed_res_wei_pdf_gw_32_method = rep(0,B)
res_wei_pdf_gw_42_method = rep(0,B)
smoothed_res_wei_pdf_gw_42_method = rep(0,B)
res_wei_pdf_gw_52_method = rep(0,B)
smoothed_res_wei_pdf_gw_52_method = rep(0,B)

res_wei_pdf_gw_12_method2 = rep(0,B)
res_wei_pdf_gw_22_method2 = rep(0,B)
res_wei_pdf_gw_32_method2 = rep(0,B)
res_wei_pdf_gw_42_method2 = rep(0,B)
res_wei_pdf_gw_52_method2 = rep(0,B)

res_wei_pdf_gw_12_method3 = rep(0,B)
res_wei_pdf_gw_22_method3 = rep(0,B)
res_wei_pdf_gw_32_method3 = rep(0,B)
res_wei_pdf_gw_42_method3 = rep(0,B)
res_wei_pdf_gw_52_method3 = rep(0,B)

res_wei_pdf_gw_12_method4 = rep(0,B)
res_wei_pdf_gw_22_method4 = rep(0,B)
res_wei_pdf_gw_32_method4 = rep(0,B)
res_wei_pdf_gw_42_method4 = rep(0,B)
res_wei_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x = sort(rweibull(n,shape))
  
  # Weibull Distribution - GW1
  
  # res_wei_pdf_gw_11_method[i] = L1_Distance_calc(x,weibull_grids[[1]],pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grids[[1]],pdf_weibull)
  res_wei_pdf_gw_12_method[i] = err_malc[1]
  smoothed_res_wei_pdf_gw_12_method[i] = err_malc[2]
  res_wei_pdf_gw_12_method2[i] = L2_from_res(x,weibull_grids[[1]],pdf_weibull)
  res_wei_pdf_gw_12_method3[i] = L2_binnednp(x,weibull_grids[[1]],pdf_weibull)
  res_wei_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,weibull_grids[[1]],pdf_weibull)
  
  # Weibull Distribution - GW2
  
  # res_wei_pdf_gw_21_method[i] = L1_Distance_calc(x,weibull_grids[[2]],pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grids[[2]],pdf_weibull)
  res_wei_pdf_gw_22_method[i] = err_malc[1]
  smoothed_res_wei_pdf_gw_22_method[i] = err_malc[2]
  res_wei_pdf_gw_22_method2[i] = L2_from_res(x,weibull_grids[[2]],pdf_weibull)
  res_wei_pdf_gw_22_method3[i] = L2_binnednp(x,weibull_grids[[2]],pdf_weibull)
  res_wei_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,weibull_grids[[2]],pdf_weibull)
  
  # Weibull Distribution - GW3
  
  # res_wei_pdf_gw_31_method[i] = L1_Distance_calc(x,weibull_grids[[3]],pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grids[[3]],pdf_weibull)
  res_wei_pdf_gw_32_method[i] = err_malc[1]
  smoothed_res_wei_pdf_gw_32_method[i] = err_malc[2]
  res_wei_pdf_gw_32_method2[i] = L2_from_res(x,weibull_grids[[3]],pdf_weibull)
  res_wei_pdf_gw_32_method3[i] = L2_binnednp(x,weibull_grids[[3]],pdf_weibull)
  res_wei_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,weibull_grids[[3]],pdf_weibull)
  
  # Weibull Distribution - GW4
  
  # res_wei_pdf_gw_41_method[i] = L1_Distance_calc(x,weibull_grids[[4]],pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grids[[4]],pdf_weibull)
  res_wei_pdf_gw_42_method[i] = err_malc[1]
  smoothed_res_wei_pdf_gw_42_method[i] = err_malc[2]
  res_wei_pdf_gw_42_method2[i] = L2_from_res(x,weibull_grids[[4]],pdf_weibull)
  res_wei_pdf_gw_42_method3[i] = L2_binnednp(x,weibull_grids[[4]],pdf_weibull)
  res_wei_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,weibull_grids[[4]],pdf_weibull)
  
  # Weibull Distribution - GW5
  
  # res_wei_pdf_gw_51_method[i] = L1_Distance_calc(x,weibull_grids[[5]],pdf_weibull,range_wei)
  err_malc = L2_Distance_calc_both(x,weibull_grids[[5]],pdf_weibull)
  res_wei_pdf_gw_52_method[i] = err_malc[1]
  smoothed_res_wei_pdf_gw_52_method[i] = err_malc[2]
  res_wei_pdf_gw_52_method2[i] = L2_from_res(x,weibull_grids[[5]],pdf_weibull)
  res_wei_pdf_gw_52_method3[i] = L2_binnednp(x,weibull_grids[[5]],pdf_weibull)
  res_wei_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,weibull_grids[[5]],pdf_weibull)
  
  print(i)
  
}

# Method # 

boxplot(res_wei_pdf_gw_11_method,res_wei_pdf_gw_21_method,res_wei_pdf_gw_31_method,res_wei_pdf_gw_41_method,res_wei_pdf_gw_51_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.9))

boxplot(res_wei_pdf_gw_12_method,res_wei_pdf_gw_22_method,res_wei_pdf_gw_32_method,res_wei_pdf_gw_42_method,res_wei_pdf_gw_52_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_wei_pdf_gw_12_method2,res_wei_pdf_gw_22_method2,res_wei_pdf_gw_32_method2,res_wei_pdf_gw_42_method2,res_wei_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_wei_pdf_gw_12_method3,res_wei_pdf_gw_22_method3,res_wei_pdf_gw_32_method3,res_wei_pdf_gw_42_method3,res_wei_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_wei_pdf_gw_12_method4,res_wei_pdf_gw_22_method4,res_wei_pdf_gw_32_method4,res_wei_pdf_gw_42_method4,res_wei_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

# Pareto Distribution - Grid Width Simulation 

MM <- ceiling(qpareto(0.9999,location,shape_p)) * 1.5
pareto_grid <- seq(location,MM, by = pareto_deltas[1])

pareto_grids = list(seq(location,MM, by = pareto_deltas[5]),
                    seq(location,MM, by = pareto_deltas[4]),
                    seq(location,MM, by = pareto_deltas[3]),
                    seq(location,MM, by = pareto_deltas[2]),
                    seq(location,MM, by = pareto_deltas[1]))

# Method # 

res_pareto_pdf_gw_11_method = rep(0,B)
res_pareto_pdf_gw_21_method = rep(0,B)
res_pareto_pdf_gw_31_method = rep(0,B)
res_pareto_pdf_gw_41_method = rep(0,B)
res_pareto_pdf_gw_51_method = rep(0,B)

res_pareto_pdf_gw_12_method = rep(0,B)
smoothed_res_pareto_pdf_gw_12_method = rep(0,B)
res_pareto_pdf_gw_22_method = rep(0,B)
smoothed_res_pareto_pdf_gw_22_method = rep(0,B)
res_pareto_pdf_gw_32_method = rep(0,B)
smoothed_res_pareto_pdf_gw_32_method = rep(0,B)
res_pareto_pdf_gw_42_method = rep(0,B)
smoothed_res_pareto_pdf_gw_42_method = rep(0,B)
res_pareto_pdf_gw_52_method = rep(0,B)
smoothed_res_pareto_pdf_gw_52_method = rep(0,B)

res_pareto_pdf_gw_12_method2 = rep(0,B)
res_pareto_pdf_gw_22_method2 = rep(0,B)
res_pareto_pdf_gw_32_method2 = rep(0,B)
res_pareto_pdf_gw_42_method2 = rep(0,B)
res_pareto_pdf_gw_52_method2 = rep(0,B)

res_pareto_pdf_gw_12_method3 = rep(0,B)
res_pareto_pdf_gw_22_method3 = rep(0,B)
res_pareto_pdf_gw_32_method3 = rep(0,B)
res_pareto_pdf_gw_42_method3 = rep(0,B)
res_pareto_pdf_gw_52_method3 = rep(0,B)

res_pareto_pdf_gw_12_method4 = rep(0,B)
res_pareto_pdf_gw_22_method4 = rep(0,B)
res_pareto_pdf_gw_32_method4 = rep(0,B)
res_pareto_pdf_gw_42_method4 = rep(0,B)
res_pareto_pdf_gw_52_method4 = rep(0,B)


set.seed(10)

for (i in 1:B) {
  
  x <- sort(rpareto(n, location, shape_p))
  
  pareto_grid1 <- make_breaks_cover(x, pareto_grids[[1]])
  pareto_grid2 <- make_breaks_cover(x, pareto_grids[[2]])
  pareto_grid3 <- make_breaks_cover(x, pareto_grids[[3]])
  pareto_grid4 <- make_breaks_cover(x, pareto_grids[[4]])
  pareto_grid5 <- make_breaks_cover(x, pareto_grids[[5]])
  
  # GW1
  # res_pareto_pdf_gw_11_method[i] <- L1_Distance_calc(x, pareto_grids[[1]], pdf_pareto, range_pareto)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_gw_12_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_gw_12_method[i] = err_malc[2]
  res_pareto_pdf_gw_12_method2[i] <- L2_from_res(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_gw_12_method3[i] <- L2_binnednp(x,pareto_grid1,pdf_pareto)
  res_pareto_pdf_gw_12_method4[i] <- L2_from_kernsmooth(x,pareto_grid1,pdf_pareto)
  
  # GW2
  # res_pareto_pdf_gw_21_method[i] <- L1_Distance_calc(x, pareto_grids[[2]], pdf_pareto, range_pareto)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_gw_22_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_gw_22_method[i] = err_malc[2]
  res_pareto_pdf_gw_22_method2[i] <- L2_from_res(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_gw_22_method3[i] <- L2_binnednp(x,pareto_grid2,pdf_pareto)
  res_pareto_pdf_gw_22_method4[i] <- L2_from_kernsmooth(x,pareto_grid2,pdf_pareto)
  
  # GW3
  # res_pareto_pdf_gw_31_method[i] <- L1_Distance_calc(x, pareto_grids[[3]], pdf_pareto, range_pareto)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_gw_32_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_gw_32_method[i] = err_malc[2]
  res_pareto_pdf_gw_32_method2[i] <- L2_from_res(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_gw_32_method3[i] <- L2_binnednp(x,pareto_grid3,pdf_pareto)
  res_pareto_pdf_gw_32_method4[i] <- L2_from_kernsmooth(x,pareto_grid3,pdf_pareto)
  
  # GW4
  # res_pareto_pdf_gw_41_method[i] <- L1_Distance_calc(x, pareto_grids[[4]], pdf_pareto, range_pareto)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_gw_42_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_gw_42_method[i] = err_malc[2]
  res_pareto_pdf_gw_42_method2[i] <- L2_from_res(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_gw_42_method3[i] <- L2_binnednp(x,pareto_grid4,pdf_pareto)
  res_pareto_pdf_gw_42_method4[i] <- L2_from_kernsmooth(x,pareto_grid4,pdf_pareto)
  
  # GW5
  # res_pareto_pdf_gw_51_method[i] <- L1_Distance_calc(x, pareto_grids[[5]], pdf_pareto, range_pareto)
  
  err_malc = L2_Distance_calc_both(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_gw_52_method[i] = err_malc[1]
  smoothed_res_pareto_pdf_gw_52_method[i] = err_malc[2]
  res_pareto_pdf_gw_52_method2[i] <- L2_from_res(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_gw_52_method3[i] <- L2_binnednp(x,pareto_grid5,pdf_pareto)
  res_pareto_pdf_gw_52_method4[i] <- L2_from_kernsmooth(x,pareto_grid5,pdf_pareto)
  
  print(i)
  
}

# Method # 


boxplot(res_pareto_pdf_gw_11_method,res_pareto_pdf_gw_21_method,res_pareto_pdf_gw_31_method,res_pareto_pdf_gw_41_method,res_pareto_pdf_gw_51_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.7))

boxplot(res_pareto_pdf_gw_12_method,res_pareto_pdf_gw_22_method,res_pareto_pdf_gw_32_method,res_pareto_pdf_gw_42_method,res_pareto_pdf_gw_52_method,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4))

boxplot(res_pareto_pdf_gw_12_method2,res_pareto_pdf_gw_22_method2,res_pareto_pdf_gw_32_method2,res_pareto_pdf_gw_42_method2,res_pareto_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4))

boxplot(res_pareto_pdf_gw_12_method3,res_pareto_pdf_gw_22_method3,res_pareto_pdf_gw_32_method3,res_pareto_pdf_gw_42_method3,res_pareto_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4))

boxplot(res_pareto_pdf_gw_12_method4,res_pareto_pdf_gw_22_method4,res_pareto_pdf_gw_32_method4,res_pareto_pdf_gw_42_method4,res_pareto_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4))

# Data Save

# Create the data frame
df_lap_n_method1 <- data.frame(
  `1e2` = res_lap_pdf_n_12_method,
  `1e3` = res_lap_pdf_n_22_method,
  `1e4` = res_lap_pdf_n_32_method,
  `1e5` = res_lap_pdf_n_42_method,
  `1e6` = res_lap_pdf_n_52_method
)

df_lap_n_method2 <- data.frame(
  `1e2` = res_lap_pdf_n_12_method2,
  `1e3` = res_lap_pdf_n_22_method2,
  `1e4` = res_lap_pdf_n_32_method2,
  `1e5` = res_lap_pdf_n_42_method2,
  `1e6` = res_lap_pdf_n_52_method2
)

df_lap_n_method3 <- data.frame(
  `1e2` = res_lap_pdf_n_12_method3,
  `1e3` = res_lap_pdf_n_22_method3,
  `1e4` = res_lap_pdf_n_32_method3,
  `1e5` = res_lap_pdf_n_42_method3,
  `1e6` = res_lap_pdf_n_52_method3
)

df_lap_n_method4 <- data.frame(
  `1e2` = res_lap_pdf_n_12_method4,
  `1e3` = res_lap_pdf_n_22_method4,
  `1e4` = res_lap_pdf_n_32_method4,
  `1e5` = res_lap_pdf_n_42_method4,
  `1e6` = res_lap_pdf_n_52_method4
)

# Create the data frame
df_chisq_n_method1 <- data.frame(
  `1e2` = res_chisq_pdf_n_12_method,
  `1e3` = res_chisq_pdf_n_22_method,
  `1e4` = res_chisq_pdf_n_32_method,
  `1e5` = res_chisq_pdf_n_42_method,
  `1e6` = res_chisq_pdf_n_52_method
)

df_chisq_n_method2 <- data.frame(
  `1e2` = res_chisq_pdf_n_12_method2,
  `1e3` = res_chisq_pdf_n_22_method2,
  `1e4` = res_chisq_pdf_n_32_method2,
  `1e5` = res_chisq_pdf_n_42_method2,
  `1e6` = res_chisq_pdf_n_52_method2
)

df_chisq_n_method3 <- data.frame(
  `1e2` = res_chisq_pdf_n_12_method3,
  `1e3` = res_chisq_pdf_n_22_method3,
  `1e4` = res_chisq_pdf_n_32_method3,
  `1e5` = res_chisq_pdf_n_42_method3,
  `1e6` = res_chisq_pdf_n_52_method3
)

df_chisq_n_method4 <- data.frame(
  `1e2` = res_chisq_pdf_n_12_method4,
  `1e3` = res_chisq_pdf_n_22_method4,
  `1e4` = res_chisq_pdf_n_32_method4,
  `1e5` = res_chisq_pdf_n_42_method4,
  `1e6` = res_chisq_pdf_n_52_method4
)

# Create the data frame
df_lnorm_n_method1 <- data.frame(
  `1e2` = res_logn_pdf_n_12_method,
  `1e3` = res_logn_pdf_n_22_method,
  `1e4` = res_logn_pdf_n_32_method,
  `1e5` = res_logn_pdf_n_42_method,
  `1e6` = res_logn_pdf_n_52_method
)

df_lnorm_n_method2 <- data.frame(
  `1e2` = res_logn_pdf_n_12_method2,
  `1e3` = res_logn_pdf_n_22_method2,
  `1e4` = res_logn_pdf_n_32_method2,
  `1e5` = res_logn_pdf_n_42_method2,
  `1e6` = res_logn_pdf_n_52_method2
)

df_lnorm_n_method3 <- data.frame(
  `1e2` = res_logn_pdf_n_12_method3,
  `1e3` = res_logn_pdf_n_22_method3,
  `1e4` = res_logn_pdf_n_32_method3,
  `1e5` = res_logn_pdf_n_42_method3,
  `1e6` = res_logn_pdf_n_52_method3
)

df_lnorm_n_method4 <- data.frame(
  `1e2` = res_logn_pdf_n_12_method4,
  `1e3` = res_logn_pdf_n_22_method4,
  `1e4` = res_logn_pdf_n_32_method4,
  `1e5` = res_logn_pdf_n_42_method4,
  `1e6` = res_logn_pdf_n_52_method4
)


# Create the data frame
df_wei_n_method1 <- data.frame(
  `1e2` = res_wei_pdf_n_12_method,
  `1e3` = res_wei_pdf_n_22_method,
  `1e4` = res_wei_pdf_n_32_method,
  `1e5` = res_wei_pdf_n_42_method,
  `1e6` = res_wei_pdf_n_52_method
)

df_wei_n_method2 <- data.frame(
  `1e2` = res_wei_pdf_n_12_method2,
  `1e3` = res_wei_pdf_n_22_method2,
  `1e4` = res_wei_pdf_n_32_method2,
  `1e5` = res_wei_pdf_n_42_method2,
  `1e6` = res_wei_pdf_n_52_method2
)

df_wei_n_method3 <- data.frame(
  `1e2` = res_wei_pdf_n_12_method3,
  `1e3` = res_wei_pdf_n_22_method3,
  `1e4` = res_wei_pdf_n_32_method3,
  `1e5` = res_wei_pdf_n_42_method3,
  `1e6` = res_wei_pdf_n_52_method3
)

df_wei_n_method4 <- data.frame(
  `1e2` = res_wei_pdf_n_12_method4,
  `1e3` = res_wei_pdf_n_22_method4,
  `1e4` = res_wei_pdf_n_32_method4,
  `1e5` = res_wei_pdf_n_42_method4,
  `1e6` = res_wei_pdf_n_52_method4
)

# Create the data frame
df_par_n_method1 <- data.frame(
  `1e2` = res_pareto_pdf_n_12_method,
  `1e3` = res_pareto_pdf_n_22_method,
  `1e4` = res_pareto_pdf_n_32_method,
  `1e5` = res_pareto_pdf_n_42_method,
  `1e6` = res_pareto_pdf_n_52_method
)

df_par_n_method2 <- data.frame(
  `1e2` = res_pareto_pdf_n_12_method2,
  `1e3` = res_pareto_pdf_n_22_method2,
  `1e4` = res_pareto_pdf_n_32_method2,
  `1e5` = res_pareto_pdf_n_42_method2,
  `1e6` = res_pareto_pdf_n_52_method2
)

df_par_n_method3 <- data.frame(
  `1e2` = res_pareto_pdf_n_12_method3,
  `1e3` = res_pareto_pdf_n_22_method3,
  `1e4` = res_pareto_pdf_n_32_method3,
  `1e5` = res_pareto_pdf_n_42_method3,
  `1e6` = res_pareto_pdf_n_52_method3
)

df_par_n_method4 <- data.frame(
  `1e2` = res_pareto_pdf_n_12_method4,
  `1e3` = res_pareto_pdf_n_22_method4,
  `1e4` = res_pareto_pdf_n_32_method4,
  `1e5` = res_pareto_pdf_n_42_method4,
  `1e6` = res_pareto_pdf_n_52_method4
)

df_lap_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_lap_pdf_n_12_method,
  `1e3` = smoothed_res_lap_pdf_n_22_method,
  `1e4` = smoothed_res_lap_pdf_n_32_method,
  `1e5` = smoothed_res_lap_pdf_n_42_method,
  `1e6` = smoothed_res_lap_pdf_n_52_method
)

df_chisq_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_chisq_pdf_n_12_method,
  `1e3` = smoothed_res_chisq_pdf_n_22_method,
  `1e4` = smoothed_res_chisq_pdf_n_32_method,
  `1e5` = smoothed_res_chisq_pdf_n_42_method,
  `1e6` = smoothed_res_chisq_pdf_n_52_method
)

df_lnorm_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_logn_pdf_n_12_method,
  `1e3` = smoothed_res_logn_pdf_n_22_method,
  `1e4` = smoothed_res_logn_pdf_n_32_method,
  `1e5` = smoothed_res_logn_pdf_n_42_method,
  `1e6` = smoothed_res_logn_pdf_n_52_method
)

df_wei_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_wei_pdf_n_12_method,
  `1e3` = smoothed_res_wei_pdf_n_22_method,
  `1e4` = smoothed_res_wei_pdf_n_32_method,
  `1e5` = smoothed_res_wei_pdf_n_42_method,
  `1e6` = smoothed_res_wei_pdf_n_52_method
)

df_par_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_pareto_pdf_n_12_method,
  `1e3` = smoothed_res_pareto_pdf_n_22_method,
  `1e4` = smoothed_res_pareto_pdf_n_32_method,
  `1e5` = smoothed_res_pareto_pdf_n_42_method,
  `1e6` = smoothed_res_pareto_pdf_n_52_method
)

write.csv(df_lap_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_chisq_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_lnorm_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_wei_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_par_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method1.csv", row.names = FALSE)

write.csv(df_lap_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_chisq_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_lnorm_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_wei_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_par_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method1_smoothed.csv", row.names = FALSE)

write.csv(df_lap_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_chisq_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_lnorm_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_wei_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_par_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method2.csv", row.names = FALSE)

write.csv(df_lap_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_chisq_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_lnorm_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_wei_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_par_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method3.csv", row.names = FALSE)

write.csv(df_lap_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_chisq_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_lnorm_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_wei_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_par_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method4.csv", row.names = FALSE)

# Create the data frame
df_lap_gw_method1 <- data.frame(
  `2.5` = res_lap_pdf_gw_12_method,
  `2.0` = res_lap_pdf_gw_22_method,
  `1.5` = res_lap_pdf_gw_32_method,
  `1.0` = res_lap_pdf_gw_42_method,
  `0.5` = res_lap_pdf_gw_52_method
)

df_lap_gw_method2 <- data.frame(
  `2.5` = res_lap_pdf_gw_12_method2,
  `2.0` = res_lap_pdf_gw_22_method2,
  `1.5` = res_lap_pdf_gw_32_method2,
  `1.0` = res_lap_pdf_gw_42_method2,
  `0.5` = res_lap_pdf_gw_52_method2
)

df_lap_gw_method3 <- data.frame(
  `2.5` = res_lap_pdf_gw_12_method3,
  `2.0` = res_lap_pdf_gw_22_method3,
  `1.5` = res_lap_pdf_gw_32_method3,
  `1.0` = res_lap_pdf_gw_42_method3,
  `0.5` = res_lap_pdf_gw_52_method3
)

df_lap_gw_method4 <- data.frame(
  `2.5` = res_lap_pdf_gw_12_method4,
  `2.0` = res_lap_pdf_gw_22_method4,
  `1.5` = res_lap_pdf_gw_32_method4,
  `1.0` = res_lap_pdf_gw_42_method4,
  `0.5` = res_lap_pdf_gw_52_method4
)

# Create the data frame
df_chisq_gw_method1 <- data.frame(
  `2.5` = res_chisq_pdf_gw_12_method,
  `2.0` = res_chisq_pdf_gw_22_method,
  `1.5` = res_chisq_pdf_gw_32_method,
  `1.0` = res_chisq_pdf_gw_42_method,
  `0.5` = res_chisq_pdf_gw_52_method
)

df_chisq_gw_method2 <- data.frame(
  `2.5` = res_chisq_pdf_gw_12_method2,
  `2.0` = res_chisq_pdf_gw_22_method2,
  `1.5` = res_chisq_pdf_gw_32_method2,
  `1.0` = res_chisq_pdf_gw_42_method2,
  `0.5` = res_chisq_pdf_gw_52_method2
)

df_chisq_gw_method3 <- data.frame(
  `2.5` = res_chisq_pdf_gw_12_method3,
  `2.0` = res_chisq_pdf_gw_22_method3,
  `1.5` = res_chisq_pdf_gw_32_method3,
  `1.0` = res_chisq_pdf_gw_42_method3,
  `0.5` = res_chisq_pdf_gw_52_method3
)

df_chisq_gw_method4 <- data.frame(
  `2.5` = res_chisq_pdf_gw_12_method4,
  `2.0` = res_chisq_pdf_gw_22_method4,
  `1.5` = res_chisq_pdf_gw_32_method4,
  `1.0` = res_chisq_pdf_gw_42_method4,
  `0.5` = res_chisq_pdf_gw_52_method4
)

# Create the data frame
df_lnorm_gw_method1 <- data.frame(
  `2.5` = res_lnorm_pdf_gw_12_method,
  `2.0` = res_lnorm_pdf_gw_22_method,
  `1.5` = res_lnorm_pdf_gw_32_method,
  `1.0` = res_lnorm_pdf_gw_42_method,
  `0.5` = res_lnorm_pdf_gw_52_method
)

df_lnorm_gw_method2 <- data.frame(
  `2.5` = res_lnorm_pdf_gw_12_method2,
  `2.0` = res_lnorm_pdf_gw_22_method2,
  `1.5` = res_lnorm_pdf_gw_32_method2,
  `1.0` = res_lnorm_pdf_gw_42_method2,
  `0.5` = res_lnorm_pdf_gw_52_method2
)

df_lnorm_gw_method3 <- data.frame(
  `2.5` = res_lnorm_pdf_gw_12_method3,
  `2.0` = res_lnorm_pdf_gw_22_method3,
  `1.5` = res_lnorm_pdf_gw_32_method3,
  `1.0` = res_lnorm_pdf_gw_42_method3,
  `0.5` = res_lnorm_pdf_gw_52_method3
)

df_lnorm_gw_method4 <- data.frame(
  `2.5` = res_lnorm_pdf_gw_12_method4,
  `2.0` = res_lnorm_pdf_gw_22_method4,
  `1.5` = res_lnorm_pdf_gw_32_method4,
  `1.0` = res_lnorm_pdf_gw_42_method4,
  `0.5` = res_lnorm_pdf_gw_52_method4
)

# Create the data frame
df_wei_gw_method1 <- data.frame(
  `2.5` = res_wei_pdf_gw_12_method,
  `2.0` = res_wei_pdf_gw_22_method,
  `1.5` = res_wei_pdf_gw_32_method,
  `1.0` = res_wei_pdf_gw_42_method,
  `0.5` = res_wei_pdf_gw_52_method
)

df_wei_gw_method2 <- data.frame(
  `2.5` = res_wei_pdf_gw_12_method2,
  `2.0` = res_wei_pdf_gw_22_method2,
  `1.5` = res_wei_pdf_gw_32_method2,
  `1.0` = res_wei_pdf_gw_42_method2,
  `0.5` = res_wei_pdf_gw_52_method2
)

df_wei_gw_method3 <- data.frame(
  `2.5` = res_wei_pdf_gw_12_method3,
  `2.0` = res_wei_pdf_gw_22_method3,
  `1.5` = res_wei_pdf_gw_32_method3,
  `1.0` = res_wei_pdf_gw_42_method3,
  `0.5` = res_wei_pdf_gw_52_method3
)

df_wei_gw_method4 <- data.frame(
  `2.5` = res_wei_pdf_gw_12_method4,
  `2.0` = res_wei_pdf_gw_22_method4,
  `1.5` = res_wei_pdf_gw_32_method4,
  `1.0` = res_wei_pdf_gw_42_method4,
  `0.5` = res_wei_pdf_gw_52_method4
)

# Create the data frame
df_par_gw_method1 <- data.frame(
  `2.5` = res_pareto_pdf_gw_12_method,
  `2.0` = res_pareto_pdf_gw_22_method,
  `1.5` = res_pareto_pdf_gw_32_method,
  `1.0` = res_pareto_pdf_gw_42_method,
  `0.5` = res_pareto_pdf_gw_52_method
)

df_par_gw_method2 <- data.frame(
  `2.5` = res_pareto_pdf_gw_12_method2,
  `2.0` = res_pareto_pdf_gw_22_method2,
  `1.5` = res_pareto_pdf_gw_32_method2,
  `1.0` = res_pareto_pdf_gw_42_method2,
  `0.5` = res_pareto_pdf_gw_52_method2
)

df_par_gw_method3 <- data.frame(
  `2.5` = res_pareto_pdf_gw_12_method3,
  `2.0` = res_pareto_pdf_gw_22_method3,
  `1.5` = res_pareto_pdf_gw_32_method3,
  `1.0` = res_pareto_pdf_gw_42_method3,
  `0.5` = res_pareto_pdf_gw_52_method3
)

df_par_gw_method4 <- data.frame(
  `2.5` = res_pareto_pdf_gw_12_method4,
  `2.0` = res_pareto_pdf_gw_22_method4,
  `1.5` = res_pareto_pdf_gw_32_method4,
  `1.0` = res_pareto_pdf_gw_42_method4,
  `0.5` = res_pareto_pdf_gw_52_method4
)

df_lap_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_lap_pdf_gw_12_method,
  `2.0` = smoothed_res_lap_pdf_gw_22_method,
  `1.5` = smoothed_res_lap_pdf_gw_32_method,
  `1.0` = smoothed_res_lap_pdf_gw_42_method,
  `0.5` = smoothed_res_lap_pdf_gw_52_method
)

df_chisq_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_chisq_pdf_gw_12_method,
  `2.0` = smoothed_res_chisq_pdf_gw_22_method,
  `1.5` = smoothed_res_chisq_pdf_gw_32_method,
  `1.0` = smoothed_res_chisq_pdf_gw_42_method,
  `0.5` = smoothed_res_chisq_pdf_gw_52_method
)

df_lnorm_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_lnorm_pdf_gw_12_method,
  `2.0` = smoothed_res_lnorm_pdf_gw_22_method,
  `1.5` = smoothed_res_lnorm_pdf_gw_32_method,
  `1.0` = smoothed_res_lnorm_pdf_gw_42_method,
  `0.5` = smoothed_res_lnorm_pdf_gw_52_method
)

df_wei_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_wei_pdf_gw_12_method,
  `2.0` = smoothed_res_wei_pdf_gw_22_method,
  `1.5` = smoothed_res_wei_pdf_gw_32_method,
  `1.0` = smoothed_res_wei_pdf_gw_42_method,
  `0.5` = smoothed_res_wei_pdf_gw_52_method
)

df_par_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_pareto_pdf_gw_12_method,
  `2.0` = smoothed_res_pareto_pdf_gw_22_method,
  `1.5` = smoothed_res_pareto_pdf_gw_32_method,
  `1.0` = smoothed_res_pareto_pdf_gw_42_method,
  `0.5` = smoothed_res_pareto_pdf_gw_52_method
)

write.csv(df_lap_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_chisq_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_lnorm_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_wei_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_par_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method1.csv", row.names = FALSE)

write.csv(df_lap_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_chisq_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_lnorm_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_wei_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_par_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method1_smoothed.csv", row.names = FALSE)

write.csv(df_lap_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_chisq_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_lnorm_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_wei_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_par_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method2.csv", row.names = FALSE)

write.csv(df_lap_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_chisq_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_lnorm_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_wei_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_par_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method3.csv", row.names = FALSE)

write.csv(df_lap_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_chisq_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_lnorm_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_wei_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_par_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method4.csv", row.names = FALSE)


# Recreating the plots

# Simulations Within the Distributions 

colors6 = rep("#FDB462",5)
colors7 = rep("#B3DE69",5)
colors8 = rep("#BC80BD",5)
colors9 = rep("cyan", 5)
colors10 = rep("#FB8072",5)

# Load new data -----------------------------------------------------------

df_lap_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method1.csv")
df_lap_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method2.csv")
df_lap_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method3.csv")
df_lap_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method4.csv")
df_lap_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_sample_sim_data_method1_smoothed.csv")

df_chisq_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method1.csv")
df_chisq_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method2.csv")
df_chisq_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method3.csv")
df_chisq_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method4.csv")
df_chisq_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_sample_sim_data_method1_smoothed.csv")

df_lnorm_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method1.csv")
df_lnorm_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method2.csv")
df_lnorm_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method3.csv")
df_lnorm_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method4.csv")
df_lnorm_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_sample_sim_data_method1_smoothed.csv")

df_wei_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method1.csv")
df_wei_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method2.csv")
df_wei_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method3.csv")
df_wei_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method4.csv")
df_wei_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_sample_sim_data_method1_smoothed.csv")

df_par_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method1.csv")
df_par_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method2.csv")
df_par_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method3.csv")
df_par_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method4.csv")
df_par_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_sample_sim_data_method1_smoothed.csv")

df_lap_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method1.csv")
df_lap_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method2.csv")
df_lap_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method3.csv")
df_lap_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method4.csv")
df_lap_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/laplace_grid_sim_data_method1_smoothed.csv")

df_chisq_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method1.csv")
df_chisq_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method2.csv")
df_chisq_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method3.csv")
df_chisq_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method4.csv")
df_chisq_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/chisqr_grid_sim_data_method1_smoothed.csv")

df_lnorm_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method1.csv")
df_lnorm_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method2.csv")
df_lnorm_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method3.csv")
df_lnorm_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method4.csv")
df_lnorm_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/lnorm_grid_sim_data_method1_smoothed.csv")

df_wei_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method1.csv")
df_wei_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method2.csv")
df_wei_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method3.csv")
df_wei_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method4.csv")
df_wei_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/weibull_grid_sim_data_method1_smoothed.csv")

df_par_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method1.csv")
df_par_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method2.csv")
df_par_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method3.csv")
df_par_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method4.csv")
df_par_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/pareto_grid_sim_data_method1_smoothed.csv")

plot_method1_smoothed_n <- function(df, file, cols, y_limits) {
  grDevices::pdf(file, width = 6, height = 6)
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
      font.axis=1, cex.main=2)
  boxplot(df$X1e2,
          df$X1e3,
          df$X1e4,
          df$X1e5,
          df$X1e6,
          names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
          xlab = '', col = cols, ylim = y_limits)
  dev.off()
}

plot_method1_smoothed_gw <- function(df, file, cols, y_limits) {
  grDevices::pdf(file, width = 6, height = 6)
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
      font.axis=1, cex.main=2)
  boxplot(df$X2.5,
          df$X2.0,
          df$X1.5,
          df$X1.0,
          df$X0.5,
          names = c("2.5","2.0","1.5","1.0","0.5"),
          xlab = '', col = cols, ylim = y_limits)
  dev.off()
}

############################
# METHOD 1 - SMOOTHED
############################

plot_method1_smoothed_n(
  df_lap_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Laplace-n-Method1-Smoothed.pdf",
  colors6,
  c(0,0.3)
)

plot_method1_smoothed_gw(
  df_lap_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Laplace-gw-Method1-Smoothed.pdf",
  colors6,
  c(0,0.8)
)

plot_method1_smoothed_n(
  df_chisq_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Chisq-n-Method1-Smoothed.pdf",
  colors7,
  c(0,0.2)
)

plot_method1_smoothed_gw(
  df_chisq_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Chisq-gw-Method1-Smoothed.pdf",
  colors7,
  c(0,0.4)
)

plot_method1_smoothed_n(
  df_lnorm_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Lnorm-n-Method1-Smoothed.pdf",
  colors8,
  c(0,0.4)
)

plot_method1_smoothed_gw(
  df_lnorm_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Lnorm-gw-Method1-Smoothed.pdf",
  colors8,
  c(0,0.4)
)

plot_method1_smoothed_n(
  df_wei_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Wei-n-Method1-Smoothed.pdf",
  colors9,
  c(0,0.2)
)

plot_method1_smoothed_gw(
  df_wei_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Wei-gw-Method1-Smoothed.pdf",
  colors9,
  c(0,0.5)
)

plot_method1_smoothed_n(
  df_par_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Pareto-n-Method1-Smoothed.pdf",
  colors10,
  c(0,0.7)
)

plot_method1_smoothed_gw(
  df_par_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Pareto-gw-Method1-Smoothed.pdf",
  colors10,
  c(0,4)
)

# Laplace -----------------------------------------------------------------

# Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Laplace-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_n_method1$X1e2, 
        df_lap_n_method1$X1e3, 
        df_lap_n_method1$X1e4, 
        df_lap_n_method1$X1e5, 
        df_lap_n_method1$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors6, ylim=c(0,0.3))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Laplace-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_n_method2$X1e2, 
        df_lap_n_method2$X1e3, 
        df_lap_n_method2$X1e4, 
        df_lap_n_method2$X1e5, 
        df_lap_n_method2$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors6, ylim=c(0,0.3))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Laplace-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_n_method3$X1e2, 
        df_lap_n_method3$X1e3, 
        df_lap_n_method3$X1e4, 
        df_lap_n_method3$X1e5, 
        df_lap_n_method3$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors6, ylim=c(0,0.3))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Laplace-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_n_method4$X1e2, 
        df_lap_n_method4$X1e3, 
        df_lap_n_method4$X1e4, 
        df_lap_n_method4$X1e5, 
        df_lap_n_method4$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors6, ylim=c(0,0.3))
dev.off()

# Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Laplace-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_gw_method1$X2.5, 
        df_lap_gw_method1$X2.0, 
        df_lap_gw_method1$X1.5, 
        df_lap_gw_method1$X1.0, 
        df_lap_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors6, ylim=c(0,0.8))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Laplace-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_gw_method2$X2.5, 
        df_lap_gw_method2$X2.0, 
        df_lap_gw_method2$X1.5, 
        df_lap_gw_method2$X1.0, 
        df_lap_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors6, ylim=c(0,0.8))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Laplace-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_gw_method3$X2.5, 
        df_lap_gw_method3$X2.0, 
        df_lap_gw_method3$X1.5, 
        df_lap_gw_method3$X1.0, 
        df_lap_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors6, ylim=c(0,0.8))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Laplace-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lap_gw_method4$X2.5, 
        df_lap_gw_method4$X2.0, 
        df_lap_gw_method4$X1.5, 
        df_lap_gw_method4$X1.0, 
        df_lap_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors6, ylim=c(0,0.8))
dev.off()


# Chi-square --------------------------------------------------------------

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Chisq-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_n_method1$X1e2, 
        df_chisq_n_method1$X1e3, 
        df_chisq_n_method1$X1e4, 
        df_chisq_n_method1$X1e5, 
        df_chisq_n_method1$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors7, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Chisq-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_n_method2$X1e2, 
        df_chisq_n_method2$X1e3, 
        df_chisq_n_method2$X1e4, 
        df_chisq_n_method2$X1e5, 
        df_chisq_n_method2$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors7, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Chisq-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_n_method3$X1e2, 
        df_chisq_n_method3$X1e3, 
        df_chisq_n_method3$X1e4, 
        df_chisq_n_method3$X1e5, 
        df_chisq_n_method3$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors7, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Chisq-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_n_method4$X1e2, 
        df_chisq_n_method4$X1e3, 
        df_chisq_n_method4$X1e4, 
        df_chisq_n_method4$X1e5, 
        df_chisq_n_method4$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors7, ylim=c(0,0.2))
dev.off()



grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Chisq-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_gw_method1$X2.5, 
        df_chisq_gw_method1$X2.0, 
        df_chisq_gw_method1$X1.5, 
        df_chisq_gw_method1$X1.0, 
        df_chisq_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors7, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Chisq-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_gw_method2$X2.5, 
        df_chisq_gw_method2$X2.0, 
        df_chisq_gw_method2$X1.5, 
        df_chisq_gw_method2$X1.0, 
        df_chisq_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors7, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Chisq-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_gw_method3$X2.5, 
        df_chisq_gw_method3$X2.0, 
        df_chisq_gw_method3$X1.5, 
        df_chisq_gw_method3$X1.0, 
        df_chisq_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors7, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Chisq-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_chisq_gw_method4$X2.5, 
        df_chisq_gw_method4$X2.0, 
        df_chisq_gw_method4$X1.5, 
        df_chisq_gw_method4$X1.0, 
        df_chisq_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors7, ylim=c(0,0.4))
dev.off()

# Log-normal --------------------------------------------------------------

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Lnorm-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_n_method1$X1e2, 
        df_lnorm_n_method1$X1e3, 
        df_lnorm_n_method1$X1e4, 
        df_lnorm_n_method1$X1e5, 
        df_lnorm_n_method1$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Lnorm-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_n_method2$X1e2, 
        df_lnorm_n_method2$X1e3, 
        df_lnorm_n_method2$X1e4, 
        df_lnorm_n_method2$X1e5, 
        df_lnorm_n_method2$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Lnorm-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_n_method3$X1e2, 
        df_lnorm_n_method3$X1e3, 
        df_lnorm_n_method3$X1e4, 
        df_lnorm_n_method3$X1e5, 
        df_lnorm_n_method3$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Lnorm-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_n_method4$X1e2, 
        df_lnorm_n_method4$X1e3, 
        df_lnorm_n_method4$X1e4, 
        df_lnorm_n_method4$X1e5, 
        df_lnorm_n_method4$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Lnorm-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_gw_method1$X2.5, 
        df_lnorm_gw_method1$X2.0, 
        df_lnorm_gw_method1$X1.5, 
        df_lnorm_gw_method1$X1.0, 
        df_lnorm_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Lnorm-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_gw_method2$X2.5, 
        df_lnorm_gw_method2$X2.0, 
        df_lnorm_gw_method2$X1.5, 
        df_lnorm_gw_method2$X1.0, 
        df_lnorm_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Lnorm-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_gw_method3$X2.5, 
        df_lnorm_gw_method3$X2.0, 
        df_lnorm_gw_method3$X1.5, 
        df_lnorm_gw_method3$X1.0, 
        df_lnorm_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors8, ylim=c(0,0.4))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Lnorm-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_lnorm_gw_method4$X2.5, 
        df_lnorm_gw_method4$X2.0, 
        df_lnorm_gw_method4$X1.5, 
        df_lnorm_gw_method4$X1.0, 
        df_lnorm_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors8, ylim=c(0,0.40))
dev.off()

# Weibull -----------------------------------------------------------------

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Wei-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_n_method1$X1e2, 
        df_wei_n_method1$X1e3, 
        df_wei_n_method1$X1e4, 
        df_wei_n_method1$X1e5, 
        df_wei_n_method1$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors9, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Wei-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_n_method2$X1e2, 
        df_wei_n_method2$X1e3, 
        df_wei_n_method2$X1e4, 
        df_wei_n_method2$X1e5, 
        df_wei_n_method2$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors9, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Wei-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_n_method3$X1e2, 
        df_wei_n_method3$X1e3, 
        df_wei_n_method3$X1e4, 
        df_wei_n_method3$X1e5, 
        df_wei_n_method3$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors9, ylim=c(0,0.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Wei-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_n_method4$X1e2, 
        df_wei_n_method4$X1e3, 
        df_wei_n_method4$X1e4, 
        df_wei_n_method4$X1e5, 
        df_wei_n_method4$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors9, ylim=c(0,0.2))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Wei-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_gw_method1$X2.5, 
        df_wei_gw_method1$X2.0, 
        df_wei_gw_method1$X1.5,
        df_wei_gw_method1$X1.0, 
        df_wei_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors9, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Wei-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_gw_method2$X2.5, 
        df_wei_gw_method2$X2.0, 
        df_wei_gw_method2$X1.5,
        df_wei_gw_method2$X1.0, 
        df_wei_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors9, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Wei-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_gw_method3$X2.5, 
        df_wei_gw_method3$X2.0, 
        df_wei_gw_method3$X1.5,
        df_wei_gw_method3$X1.0, 
        df_wei_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors9, ylim=c(0,0.4))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Wei-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_wei_gw_method4$X2.5, 
        df_wei_gw_method4$X2.0, 
        df_wei_gw_method4$X1.5,
        df_wei_gw_method4$X1.0, 
        df_wei_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors9, ylim=c(0,0.4))
dev.off()


# Pareto ------------------------------------------------------------------

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Pareto-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_n_method1$X1e2, 
        df_par_n_method1$X1e3, 
        df_par_n_method1$X1e4, 
        df_par_n_method1$X1e5, 
        df_par_n_method1$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors10, ylim=c(0,0.7))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Pareto-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_n_method2$X1e2, 
        df_par_n_method2$X1e3, 
        df_par_n_method2$X1e4, 
        df_par_n_method2$X1e5, 
        df_par_n_method2$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors10, ylim=c(0,0.7))
dev.off()


grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Pareto-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_n_method3$X1e2, 
        df_par_n_method3$X1e3, 
        df_par_n_method3$X1e4, 
        df_par_n_method3$X1e5, 
        df_par_n_method3$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors10, ylim=c(0,0.7))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Pareto-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_n_method4$X1e2, 
        df_par_n_method4$X1e3, 
        df_par_n_method4$X1e4, 
        df_par_n_method4$X1e5, 
        df_par_n_method4$X1e6,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        xlab = '', col = colors10, ylim=c(0,0.7))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Pareto-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_gw_method1$X2.5, 
        df_par_gw_method1$X2.0, 
        df_par_gw_method1$X1.5, 
        df_par_gw_method1$X1.0, 
        df_par_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors10, ylim=c(0,1.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Pareto-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_gw_method2$X2.5, 
        df_par_gw_method2$X2.0, 
        df_par_gw_method2$X1.5, 
        df_par_gw_method2$X1.0, 
        df_par_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors10, ylim=c(0,1.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Pareto-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_gw_method3$X2.5, 
        df_par_gw_method3$X2.0, 
        df_par_gw_method3$X1.5, 
        df_par_gw_method3$X1.0, 
        df_par_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors10, ylim=c(0,1.2))
dev.off()

grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Pareto-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T, cex=5, bty="l", font.lab=1, mfrow=c(1,1), cex.lab=2, cex.axis=2,
    font.axis=1, cex.main=2)
boxplot(df_par_gw_method4$X2.5, 
        df_par_gw_method4$X2.0, 
        df_par_gw_method4$X1.5, 
        df_par_gw_method4$X1.0, 
        df_par_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '', col = colors10, ylim=c(0,1.2))
dev.off()
