# Parameters

# Change the density calculation + 
# Run Normal and Gamma and send it + 
# Make a Package for all
# Refer in the paper that you wrote the code +
# Redo figure 4 and change the color according to the Figure 1 + 
# Check if they are multivariate or not. Check if they can have non-uniform grid or not + 


alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
mul = 1;
df=3;
shape_p = 4;

# Ranges

range_norm = c(-Inf,Inf)
range_beta = c(0,1)
range_gamma = c(0,Inf)
range_logistic = c(-Inf,Inf)
range_t = c(-Inf,Inf)

# Pdfs

pdf_beta <- function(x) {
  dbeta(x,alpha,beta)
}

pdf_norm <- function(x) {
  dnorm(x, mean = mu, sd = sigma)
}

pdf_gamma <- function(x) {
  dgamma(x, alpha, beta)
}

pdf_logistic <- function(x) {
  dlogis(x)
}

pdf_t <- function(x) {
  dt(x,df)
}

# Standardized Delta Values #

norm_deltas  = c(0.5, 1.0, 1.5, 2.0, 2.5)
beta_deltas  = c(0.0799, 0.1597, 0.2396, 0.3194, 0.3993)
gamma_deltas = c(0.1414, 0.2828, 0.4243, 0.5657, 0.7071)
logistic_deltas = c(0.9069, 1.8138, 2.7207, 3.6276, 4.5345)
t_deltas = c(0.8660, 1.7321, 2.5981, 3.4641, 4.3301)

# Corresponding k values to Standardized Delta Values #

norm_ks  = c(12, 6, 4, 3, 2)
beta_ks  = c(22, 11, 7, 5, 4)
gamma_ks = c(56, 28, 18, 14, 11)
logistic_ks = c(35, 17, 11, 8, 7)
t_ks = c(66, 33, 22, 16, 13)

# Simulations Within the Distribution

library(colorspace)

# Colors

colors1 <- rep("#92CAEB",5)
colors2 <- rep("#F3CF70",5)
colors3 <- rep("#66D7A5",5)
colors4 <- rep("#E6A4C6",5)
colors5 <- rep("#D55E00",5)

# Normal Distribution - Sample Size Simulation

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grid <- seq(-MM,MM, by = norm_deltas[1])

# Method-1 Residuals #

res_norm_pdf_n_11_method1 = rep(0,B)
res_norm_pdf_n_21_method1 = rep(0,B)
res_norm_pdf_n_31_method1 = rep(0,B)
res_norm_pdf_n_41_method1 = rep(0,B)
res_norm_pdf_n_51_method1 = rep(0,B)

res_norm_pdf_n_12_method1 = rep(0,B)
res_norm_pdf_n_22_method1 = rep(0,B)
res_norm_pdf_n_32_method1 = rep(0,B)
res_norm_pdf_n_42_method1 = rep(0,B)
res_norm_pdf_n_52_method1 = rep(0,B)


smoothed_res_norm_pdf_n_11_method1 = rep(0,B)
smoothed_res_norm_pdf_n_21_method1 = rep(0,B)
smoothed_res_norm_pdf_n_31_method1 = rep(0,B)
smoothed_res_norm_pdf_n_41_method1 = rep(0,B)
smoothed_res_norm_pdf_n_51_method1 = rep(0,B)

smoothed_res_norm_pdf_n_12_method1 = rep(0,B)
smoothed_res_norm_pdf_n_22_method1 = rep(0,B)
smoothed_res_norm_pdf_n_32_method1 = rep(0,B)
smoothed_res_norm_pdf_n_42_method1 = rep(0,B)
smoothed_res_norm_pdf_n_52_method1 = rep(0,B)

# Method-2 Residuals #

res_norm_pdf_n_12_method2 = rep(0,B)
res_norm_pdf_n_22_method2 = rep(0,B)
res_norm_pdf_n_32_method2 = rep(0,B)
res_norm_pdf_n_42_method2 = rep(0,B)
res_norm_pdf_n_52_method2 = rep(0,B)

# Method-3 Residuals #

res_norm_pdf_n_12_method3 = rep(0,B)
res_norm_pdf_n_22_method3 = rep(0,B)
res_norm_pdf_n_32_method3 = rep(0,B)
res_norm_pdf_n_42_method3 = rep(0,B)
res_norm_pdf_n_52_method3 = rep(0,B)

res_norm_pdf_n_12_method32 = rep(0,B)
res_norm_pdf_n_22_method32 = rep(0,B)
res_norm_pdf_n_32_method32 = rep(0,B)
res_norm_pdf_n_42_method32 = rep(0,B)
res_norm_pdf_n_52_method32 = rep(0,B)

# Method-4 Residuals #

res_norm_pdf_n_12_method4 = rep(0,B)
res_norm_pdf_n_22_method4 = rep(0,B)
res_norm_pdf_n_32_method4 = rep(0,B)
res_norm_pdf_n_42_method4 = rep(0,B)
res_norm_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Normal Distribution - n1
  
  x <- rnorm(n[1], mean = mu,sd=sigma)
  
  err_malc = L2_Distance_calc_both(x,norm_grid,pdf_norm)
  res_norm_pdf_n_12_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_n_12_method1[i] = err_malc[2]
  
  res_norm_pdf_n_12_method2[i] = L2_from_res(x,norm_grid,pdf_norm)
  res_norm_pdf_n_12_method3[i] = L2_binnednp(x,norm_grid,pdf_norm)
  res_norm_pdf_n_12_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm)
  
  # Normal Distribution - n2
  
  x <- rnorm(n[2], mean = mu,sd=sigma)
  
  err_malc = L2_Distance_calc_both(x,norm_grid,pdf_norm)
  res_norm_pdf_n_22_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_n_22_method1[i] = err_malc[2]
  
  res_norm_pdf_n_22_method2[i] = L2_from_res(x,norm_grid,pdf_norm)
  res_norm_pdf_n_22_method3[i] = L2_binnednp(x,norm_grid,pdf_norm)
  res_norm_pdf_n_22_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm)
  
  # Normal Distribution - n3
  
  x <- rnorm(n[3], mean = mu,sd=sigma)
  
  err_malc = L2_Distance_calc_both(x,norm_grid,pdf_norm)
  res_norm_pdf_n_32_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_n_32_method1[i] = err_malc[2]
  
  res_norm_pdf_n_32_method2[i] = L2_from_res(x,norm_grid,pdf_norm)
  res_norm_pdf_n_32_method3[i] = L2_binnednp(x,norm_grid,pdf_norm)
  res_norm_pdf_n_32_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm)
  
  # Normal Distribution - n4
  
  x <- rnorm(n[4], mean = mu,sd=sigma)
  
  err_malc = L2_Distance_calc_both(x,norm_grid,pdf_norm)
  res_norm_pdf_n_42_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_n_42_method1[i] = err_malc[2]
  
  res_norm_pdf_n_42_method2[i] = L2_from_res(x,norm_grid,pdf_norm)
  res_norm_pdf_n_42_method3[i] = L2_binnednp(x,norm_grid,pdf_norm)
  res_norm_pdf_n_42_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm)
  
  # Normal Distribution - n5
  
  x <- rnorm(n[5], mean = mu,sd=sigma)
  
  err_malc = L2_Distance_calc_both(x,norm_grid,pdf_norm)
  res_norm_pdf_n_52_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_n_52_method1[i] = err_malc[2]
  
  res_norm_pdf_n_52_method2[i] = L2_from_res(x,norm_grid,pdf_norm)
  res_norm_pdf_n_52_method3[i] = L2_binnednp(x,norm_grid,pdf_norm)
  res_norm_pdf_n_52_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm)
  
  print(i)
  
}

# Method-1 #

par(mfrow = c(3,2))

boxplot(res_norm_pdf_n_12_method1,
        res_norm_pdf_n_22_method1,
        res_norm_pdf_n_32_method1,
        res_norm_pdf_n_42_method1,
        res_norm_pdf_n_52_method1,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        col = colors1,ylim=c(0,0.2))

boxplot(smoothed_res_norm_pdf_n_12_method1,
        smoothed_res_norm_pdf_n_22_method1,
        smoothed_res_norm_pdf_n_32_method1,
        smoothed_res_norm_pdf_n_42_method1,
        smoothed_res_norm_pdf_n_52_method1,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        col = colors1,ylim=c(0,0.2))

boxplot(res_norm_pdf_n_12_method2,
        res_norm_pdf_n_22_method2,
        res_norm_pdf_n_32_method2,
        res_norm_pdf_n_42_method2,
        res_norm_pdf_n_52_method2,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        col = colors1,ylim=c(0,0.2))

boxplot(res_norm_pdf_n_12_method3,
        res_norm_pdf_n_22_method3,
        res_norm_pdf_n_32_method3,
        res_norm_pdf_n_42_method3,
        res_norm_pdf_n_52_method3,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        col = colors1,ylim=c(0,0.2))

boxplot(res_norm_pdf_n_12_method4,
        res_norm_pdf_n_22_method4,
        res_norm_pdf_n_32_method4,
        res_norm_pdf_n_42_method4,
        res_norm_pdf_n_52_method4,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        col = colors1,ylim=c(0,0.2))


# Beta Distribution - Sample Size Simulation

beta_grid <- seq(0,1.04, by = beta_deltas[1])

# Method 1 #

res_beta_pdf_n_11_method1 = rep(0,B)
res_beta_pdf_n_21_method1 = rep(0,B)
res_beta_pdf_n_31_method1 = rep(0,B)
res_beta_pdf_n_41_method1 = rep(0,B)
res_beta_pdf_n_51_method1 = rep(0,B)

res_beta_pdf_n_12_method1 = rep(0,B)
res_beta_pdf_n_22_method1 = rep(0,B)
res_beta_pdf_n_32_method1 = rep(0,B)
res_beta_pdf_n_42_method1 = rep(0,B)
res_beta_pdf_n_52_method1 = rep(0,B)

smoothed_res_beta_pdf_n_11_method1 = rep(0,B)
smoothed_res_beta_pdf_n_21_method1 = rep(0,B)
smoothed_res_beta_pdf_n_31_method1 = rep(0,B)
smoothed_res_beta_pdf_n_41_method1 = rep(0,B)
smoothed_res_beta_pdf_n_51_method1 = rep(0,B)

smoothed_res_beta_pdf_n_12_method1 = rep(0,B)
smoothed_res_beta_pdf_n_22_method1 = rep(0,B)
smoothed_res_beta_pdf_n_32_method1 = rep(0,B)
smoothed_res_beta_pdf_n_42_method1 = rep(0,B)
smoothed_res_beta_pdf_n_52_method1 = rep(0,B)

# Method-2 #

res_beta_pdf_n_12_method2 = rep(0,B)
res_beta_pdf_n_22_method2 = rep(0,B)
res_beta_pdf_n_32_method2 = rep(0,B)
res_beta_pdf_n_42_method2 = rep(0,B)
res_beta_pdf_n_52_method2 = rep(0,B)

# Method-3 #

res_beta_pdf_n_12_method3 = rep(0,B)
res_beta_pdf_n_22_method3 = rep(0,B)
res_beta_pdf_n_32_method3 = rep(0,B)
res_beta_pdf_n_42_method3 = rep(0,B)
res_beta_pdf_n_52_method3 = rep(0,B)

# Method-4 #

res_beta_pdf_n_12_method4 = rep(0,B)
res_beta_pdf_n_22_method4 = rep(0,B)
res_beta_pdf_n_32_method4 = rep(0,B)
res_beta_pdf_n_42_method4 = rep(0,B)
res_beta_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Beta Distribution - n1
  
  x   <- rbeta(n[1],alpha,beta)
  
  err_malc = L2_Distance_calc_both(x,beta_grid,pdf_beta)
  res_beta_pdf_n_12_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_n_12_method1[i] = err_malc[2]
  
  res_beta_pdf_n_12_method2[i] = L2_from_res(x,beta_grid,pdf_beta)
  res_beta_pdf_n_12_method3[i] = L2_binnednp(x,beta_grid,pdf_beta)
  res_beta_pdf_n_12_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta)
  
  # Beta Distribution - n2
  
  x   <- rbeta(n[2],alpha,beta)
  
  err_malc = L2_Distance_calc_both(x,beta_grid,pdf_beta)
  res_beta_pdf_n_22_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_n_22_method1[i] = err_malc[2]
  
  res_beta_pdf_n_22_method2[i] = L2_from_res(x,beta_grid,pdf_beta)
  res_beta_pdf_n_22_method3[i] = L2_binnednp(x,beta_grid,pdf_beta)
  res_beta_pdf_n_22_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta)
  
  # Beta Distribution - n3
  
  x   <- rbeta(n[3],alpha,beta)
  
  err_malc = L2_Distance_calc_both(x,beta_grid,pdf_beta)
  res_beta_pdf_n_32_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_n_32_method1[i] = err_malc[2]
  
  res_beta_pdf_n_32_method2[i] = L2_from_res(x,beta_grid,pdf_beta)
  res_beta_pdf_n_32_method3[i] = L2_binnednp(x,beta_grid,pdf_beta)
  res_beta_pdf_n_32_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta)
  
  # Beta Distribution - n4
  
  x   <- rbeta(n[4],alpha,beta)
  
  err_malc = L2_Distance_calc_both(x,beta_grid,pdf_beta)
  res_beta_pdf_n_42_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_n_42_method1[i] = err_malc[2]
  
  res_beta_pdf_n_42_method2[i] = L2_from_res(x,beta_grid,pdf_beta)
  res_beta_pdf_n_42_method3[i] = L2_binnednp(x,beta_grid,pdf_beta)
  res_beta_pdf_n_42_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta)
  
  # Beta Distribution - n5
  
  x   <- rbeta(n[5],alpha,beta)
  
  err_malc = L2_Distance_calc_both(x,beta_grid,pdf_beta)
  res_beta_pdf_n_52_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_n_52_method1[i] = err_malc[2]
  
  res_beta_pdf_n_52_method2[i] = L2_from_res(x,beta_grid,pdf_beta)
  res_beta_pdf_n_52_method3[i] = L2_binnednp(x,beta_grid,pdf_beta)
  res_beta_pdf_n_52_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta)
  
  print(i)
  
}

# Method #

boxplot(res_beta_pdf_n_12_method1,res_beta_pdf_n_22_method1,res_beta_pdf_n_32_method1,res_beta_pdf_n_42_method1,res_beta_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))

boxplot(smoothed_res_beta_pdf_n_12_method1,
        smoothed_res_beta_pdf_n_22_method1,
        smoothed_res_beta_pdf_n_32_method1,
        smoothed_res_beta_pdf_n_42_method1,
        smoothed_res_beta_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))

boxplot(res_beta_pdf_n_12_method2,res_beta_pdf_n_22_method2,res_beta_pdf_n_32_method2,res_beta_pdf_n_42_method2,res_beta_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))

boxplot(res_beta_pdf_n_12_method3,res_beta_pdf_n_22_method3,res_beta_pdf_n_32_method3,res_beta_pdf_n_42_method3,res_beta_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1))

boxplot(res_beta_pdf_n_12_method4,res_beta_pdf_n_22_method4,res_beta_pdf_n_32_method4,res_beta_pdf_n_42_method4,res_beta_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))

# Gamma Distribution - Sample Size Simulation

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grid <- seq(0,MM, by = gamma_deltas[1])

# Method-1 #

res_gamma_pdf_n_11_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_11_method1 = rep(0,B)
res_gamma_pdf_n_21_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_21_method1 = rep(0,B)
res_gamma_pdf_n_31_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_31_method1 = rep(0,B)
res_gamma_pdf_n_41_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_41_method1 = rep(0,B)
res_gamma_pdf_n_51_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_51_method1 = rep(0,B)

res_gamma_pdf_n_12_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_12_method1 = rep(0,B)
res_gamma_pdf_n_22_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_22_method1 = rep(0,B)
res_gamma_pdf_n_32_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_32_method1 = rep(0,B)
res_gamma_pdf_n_42_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_42_method1 = rep(0,B)
res_gamma_pdf_n_52_method1 = rep(0,B)
smoothed_res_gamma_pdf_n_52_method1 = rep(0,B)

# Method-2 #

res_gamma_pdf_n_12_method2 = rep(0,B)
res_gamma_pdf_n_22_method2 = rep(0,B)
res_gamma_pdf_n_32_method2 = rep(0,B)
res_gamma_pdf_n_42_method2 = rep(0,B)
res_gamma_pdf_n_52_method2 = rep(0,B)

# Method-3 #

res_gamma_pdf_n_12_method3 = rep(0,B)
res_gamma_pdf_n_22_method3 = rep(0,B)
res_gamma_pdf_n_32_method3 = rep(0,B)
res_gamma_pdf_n_42_method3 = rep(0,B)
res_gamma_pdf_n_52_method3 = rep(0,B)

# Method-4 #

res_gamma_pdf_n_12_method4 = rep(0,B)
res_gamma_pdf_n_22_method4 = rep(0,B)
res_gamma_pdf_n_32_method4 = rep(0,B)
res_gamma_pdf_n_42_method4 = rep(0,B)
res_gamma_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Gamma Distribution - n1
  
  x <- (rgamma(n[1],alpha,beta))
  
  # res_gamma_pdf_n_11_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_12_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_n_12_method1[i] = err_malc[2]
  res_gamma_pdf_n_12_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_12_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_12_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma)
  
  # Gamma Distribution - n2
  
  x   <- (rgamma(n[2],alpha,beta))
  
  # res_gamma_pdf_n_21_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_22_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_n_22_method1[i] = err_malc[2]
  res_gamma_pdf_n_22_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_22_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_22_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma)
  
  # Gamma Distribution - n3
  
  x   <- (rgamma(n[3],alpha,beta))
  
  # res_gamma_pdf_n_31_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_32_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_n_32_method1[i] = err_malc[2]
  res_gamma_pdf_n_32_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_32_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_32_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma)
  
  # Gamma Distribution - n4
  
  x   <- (rgamma(n[4],alpha,beta))
  
  # res_gamma_pdf_n_41_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_42_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_n_42_method1[i] = err_malc[2]
  res_gamma_pdf_n_42_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_42_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_42_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma)
  
  # Gamma Distribution - n5
  
  x   <- (rgamma(n[5],alpha,beta))
  
  # res_gamma_pdf_n_51_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_52_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_n_52_method1[i] = err_malc[2]
  res_gamma_pdf_n_52_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_52_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma)
  res_gamma_pdf_n_52_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma)
  
  print(i)
  
}

# Method-1 #

boxplot(smoothed_res_gamma_pdf_n_12_method1,smoothed_res_gamma_pdf_n_22_method1,
        smoothed_res_gamma_pdf_n_32_method1,smoothed_res_gamma_pdf_n_42_method1,
        smoothed_res_gamma_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_n_12_method1,res_gamma_pdf_n_22_method1,res_gamma_pdf_n_32_method1,res_gamma_pdf_n_42_method1,res_gamma_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_n_12_method2,res_gamma_pdf_n_22_method2,res_gamma_pdf_n_32_method2,res_gamma_pdf_n_42_method2,res_gamma_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_n_12_method3,res_gamma_pdf_n_22_method3,res_gamma_pdf_n_32_method3,res_gamma_pdf_n_42_method3,res_gamma_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_n_12_method4,res_gamma_pdf_n_22_method4,res_gamma_pdf_n_32_method4,res_gamma_pdf_n_42_method4,res_gamma_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))


# LaTeX summary tables for L2 errors: Normal and Gamma
# Assumes all res_* objects already exist in your R session.

make_summary_table <- function(method_vectors, stat = c("mean", "median"), digits = 4) {
  stat <- match.arg(stat)
  
  stat_fun <- switch(
    stat,
    mean = function(z) mean(z, na.rm = TRUE),
    median = function(z) median(z, na.rm = TRUE)
  )
  
  out <- sapply(method_vectors, function(vs) {
    sapply(vs, stat_fun)
  })
  
  out <- t(out)
  out <- round(out, digits)
  out
}

print_latex_table <- function(tab, caption = "", label = "", digits = 4) {
  n_cols <- colnames(tab)
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{l", paste(rep("c", length(n_cols)), collapse = ""), "}\n", sep = "")
  cat("\\hline\n")
  cat("Method & ", paste(n_cols, collapse = " & "), " \\\\\n", sep = "")
  cat("\\hline\n")
  
  for (i in seq_len(nrow(tab))) {
    vals <- formatC(tab[i, ], format = "f", digits = digits)
    cat(rownames(tab)[i], " & ", paste(vals, collapse = " & "), " \\\\\n", sep = "")
  }
  
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  if (nzchar(caption)) cat("\\caption{", caption, "}\n", sep = "")
  if (nzchar(label)) cat("\\label{", label, "}\n", sep = "")
  cat("\\end{table}\n\n")
}

n_names <- c("$10^1$", "$10^2$", "$10^3$", "$10^4$", "$10^5$")

normal_methods <- list(
  "MALC" = list(
    res_norm_pdf_n_12_method1,
    res_norm_pdf_n_22_method1,
    res_norm_pdf_n_32_method1,
    res_norm_pdf_n_42_method1,
    res_norm_pdf_n_52_method1
  ),
  "MALC Smoothed" = list(
    smoothed_res_norm_pdf_n_12_method1,
    smoothed_res_norm_pdf_n_22_method1,
    smoothed_res_norm_pdf_n_32_method1,
    smoothed_res_norm_pdf_n_42_method1,
    smoothed_res_norm_pdf_n_52_method1
  ),
  "Method 2" = list(
    res_norm_pdf_n_12_method2,
    res_norm_pdf_n_22_method2,
    res_norm_pdf_n_32_method2,
    res_norm_pdf_n_42_method2,
    res_norm_pdf_n_52_method2
  ),
  "Method 3" = list(
    res_norm_pdf_n_12_method3,
    res_norm_pdf_n_22_method3,
    res_norm_pdf_n_32_method3,
    res_norm_pdf_n_42_method3,
    res_norm_pdf_n_52_method3
  ),
  "Method 4" = list(
    res_norm_pdf_n_12_method4,
    res_norm_pdf_n_22_method4,
    res_norm_pdf_n_32_method4,
    res_norm_pdf_n_42_method4,
    res_norm_pdf_n_52_method4
  )
)

gamma_methods <- list(
  "MALC" = list(
    res_gamma_pdf_n_12_method1,
    res_gamma_pdf_n_22_method1,
    res_gamma_pdf_n_32_method1,
    res_gamma_pdf_n_42_method1,
    res_gamma_pdf_n_52_method1
  ),
  "MALC Smoothed" = list(
    smoothed_res_gamma_pdf_n_12_method1,
    smoothed_res_gamma_pdf_n_22_method1,
    smoothed_res_gamma_pdf_n_32_method1,
    smoothed_res_gamma_pdf_n_42_method1,
    smoothed_res_gamma_pdf_n_52_method1
  ),
  "Method 2" = list(
    res_gamma_pdf_n_12_method2,
    res_gamma_pdf_n_22_method2,
    res_gamma_pdf_n_32_method2,
    res_gamma_pdf_n_42_method2,
    res_gamma_pdf_n_52_method2
  ),
  "Method 3" = list(
    res_gamma_pdf_n_12_method3,
    res_gamma_pdf_n_22_method3,
    res_gamma_pdf_n_32_method3,
    res_gamma_pdf_n_42_method3,
    res_gamma_pdf_n_52_method3
  ),
  "Method 4" = list(
    res_gamma_pdf_n_12_method4,
    res_gamma_pdf_n_22_method4,
    res_gamma_pdf_n_32_method4,
    res_gamma_pdf_n_42_method4,
    res_gamma_pdf_n_52_method4
  )
)

normal_mean <- make_summary_table(normal_methods, stat = "mean")
normal_median <- make_summary_table(normal_methods, stat = "median")
gamma_mean <- make_summary_table(gamma_methods, stat = "mean")
gamma_median <- make_summary_table(gamma_methods, stat = "median")

colnames(normal_mean) <- n_names
colnames(normal_median) <- n_names
colnames(gamma_mean) <- n_names
colnames(gamma_median) <- n_names

print_latex_table(
  normal_mean,
  caption = "Mean L2 error for the normal distribution.",
  label = "tab:normal_mean_l2"
)

print_latex_table(
  normal_median,
  caption = "Median L2 error for the normal distribution.",
  label = "tab:normal_median_l2"
)

print_latex_table(
  gamma_mean,
  caption = "Mean L2 error for the gamma distribution.",
  label = "tab:gamma_mean_l2"
)

print_latex_table(
  gamma_median,
  caption = "Median L2 error for the gamma distribution.",
  label = "tab:gamma_median_l2"
)

# Logistic Distribution - Sample Size Simulation

MM = qlogis(0.9999999,mu,sigma) + 10
logistic_grid <- seq(-MM,MM, by = logistic_deltas[1])

# Method-1 #

res_logis_pdf_n_11_method1 = rep(0,B)
smoothed_res_logis_pdf_n_11_method1 = rep(0,B)
res_logis_pdf_n_21_method1 = rep(0,B)
smoothed_res_logis_pdf_n_21_method1 = rep(0,B)
res_logis_pdf_n_31_method1 = rep(0,B)
smoothed_res_logis_pdf_n_31_method1 = rep(0,B)
res_logis_pdf_n_41_method1 = rep(0,B)
smoothed_res_logis_pdf_n_41_method1 = rep(0,B)
res_logis_pdf_n_51_method1 = rep(0,B)
smoothed_res_logis_pdf_n_51_method1 = rep(0,B)

res_logis_pdf_n_12_method1 = rep(0,B)
smoothed_res_logis_pdf_n_12_method1 = rep(0,B)
res_logis_pdf_n_22_method1 = rep(0,B)
smoothed_res_logis_pdf_n_22_method1 = rep(0,B)
res_logis_pdf_n_32_method1 = rep(0,B)
smoothed_res_logis_pdf_n_32_method1 = rep(0,B)
res_logis_pdf_n_42_method1 = rep(0,B)
smoothed_res_logis_pdf_n_42_method1 = rep(0,B)
res_logis_pdf_n_52_method1 = rep(0,B)
smoothed_res_logis_pdf_n_52_method1 = rep(0,B)

# Method-2 #

res_logis_pdf_n_12_method2 = rep(0,B)
res_logis_pdf_n_22_method2 = rep(0,B)
res_logis_pdf_n_32_method2 = rep(0,B)
res_logis_pdf_n_42_method2 = rep(0,B)
res_logis_pdf_n_52_method2 = rep(0,B)

# Method-3 #

res_logis_pdf_n_12_method3 = rep(0,B)
res_logis_pdf_n_22_method3 = rep(0,B)
res_logis_pdf_n_32_method3 = rep(0,B)
res_logis_pdf_n_42_method3 = rep(0,B)
res_logis_pdf_n_52_method3 = rep(0,B)

# Method-4 #

res_logis_pdf_n_12_method4 = rep(0,B)
res_logis_pdf_n_22_method4 = rep(0,B)
res_logis_pdf_n_32_method4 = rep(0,B)
res_logis_pdf_n_42_method4 = rep(0,B)
res_logis_pdf_n_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  # Logistic Distribution - n1
  
  x = (rlogis(n[1]))
  
  # res_logis_pdf_n_11_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_12_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_n_12_method1[i] = err_malc[2]
  res_logis_pdf_n_12_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_12_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_12_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic)
  
  # Logistic Distribution - n2
  
  x = (rlogis(n[2]))
  
  # res_logis_pdf_n_21_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_22_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_n_22_method1[i] = err_malc[2]
  res_logis_pdf_n_22_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_22_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_22_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic)
  
  # Logistic Distribution - n3
  
  x = sort(rlogis(n[3]))
  
  # res_logis_pdf_n_31_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_32_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_n_32_method1[i] = err_malc[2]
  res_logis_pdf_n_32_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_32_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_32_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic)
  
  # Logistic Distribution - n4
  
  x = (rlogis(n[4]))
  
  # res_logis_pdf_n_41_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_42_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_n_42_method1[i] = err_malc[2]
  res_logis_pdf_n_42_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_42_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_42_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic)
  
  # Logistic Distribution - n5
  
  x = (rlogis(n[5]))
  
  # res_logis_pdf_n_51_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_52_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_n_52_method1[i] = err_malc[2]
  res_logis_pdf_n_52_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_52_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic)
  res_logis_pdf_n_52_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic)
  
  print(i)
  
}

# Method #

boxplot(res_logis_pdf_n_11_method1,res_logis_pdf_n_21_method1,res_logis_pdf_n_31_method1,res_logis_pdf_n_41_method1,res_logis_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))

boxplot(res_logis_pdf_n_12_method1,res_logis_pdf_n_22_method1,res_logis_pdf_n_32_method1,res_logis_pdf_n_42_method1,res_logis_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_logis_pdf_n_12_method2,res_logis_pdf_n_22_method2,res_logis_pdf_n_32_method2,res_logis_pdf_n_42_method2,res_logis_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_logis_pdf_n_12_method3,res_logis_pdf_n_22_method3,res_logis_pdf_n_32_method3,res_logis_pdf_n_42_method3,res_logis_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))

boxplot(res_logis_pdf_n_12_method4,res_logis_pdf_n_22_method4,res_logis_pdf_n_32_method4,res_logis_pdf_n_42_method4,res_logis_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))

# Student's t Distribution - Sample Size Simulation

make_breaks_cover <- function(x, grid) {
  # grid is your proposed breaks (sorted, equally spaced)
  step <- grid[2] - grid[1]
  lo <- min(grid[1],  floor(min(x)/step) * step)
  hi <- max(tail(grid,1), ceiling(max(x)/step) * step)
  seq(lo-step, hi+step, by = step)
}


MM1 <- ceiling(qt(0.9999, df)) * 1
MM2 <- ceiling(qt(0.9999, df)) * 1.5
MM3 <- ceiling(qt(0.9999, df)) * 2
MM4 <- ceiling(qt(0.9999, df)) * 2.5
MM5 <- ceiling(qt(0.9999, df)) * 3

t_grid1 <- seq(-MM1,MM1, by = t_deltas[1])
t_grid2 <- seq(-MM2,MM2, by = t_deltas[1])
t_grid3 <- seq(-MM3,MM3, by = t_deltas[1])
t_grid4 <- seq(-MM4,MM4, by = t_deltas[1])
t_grid5 <- seq(-MM5,MM5, by = t_deltas[1])

# Method-1 #

res_t_pdf_n_11_method1 = rep(0,B)
smoothed_res_t_pdf_n_11_method1 = rep(0,B)
res_t_pdf_n_21_method1 = rep(0,B)
smoothed_res_t_pdf_n_21_method1 = rep(0,B)
res_t_pdf_n_31_method1 = rep(0,B)
smoothed_res_t_pdf_n_31_method1 = rep(0,B)
res_t_pdf_n_41_method1 = rep(0,B)
smoothed_res_t_pdf_n_41_method1 = rep(0,B)
res_t_pdf_n_51_method1 = rep(0,B)
smoothed_res_t_pdf_n_51_method1 = rep(0,B)

res_t_pdf_n_12_method1 = rep(0,B)
smoothed_res_t_pdf_n_12_method1 = rep(0,B)
res_t_pdf_n_22_method1 = rep(0,B)
smoothed_res_t_pdf_n_22_method1 = rep(0,B)
res_t_pdf_n_32_method1 = rep(0,B)
smoothed_res_t_pdf_n_32_method1 = rep(0,B)
res_t_pdf_n_42_method1 = rep(0,B)
smoothed_res_t_pdf_n_42_method1 = rep(0,B)
res_t_pdf_n_52_method1 = rep(0,B)
smoothed_res_t_pdf_n_52_method1 = rep(0,B)

# Method-5 #

res_t_pdf_n_11_method5 = rep(0,B)
res_t_pdf_n_21_method5 = rep(0,B)
res_t_pdf_n_31_method5 = rep(0,B)
res_t_pdf_n_41_method5 = rep(0,B)
res_t_pdf_n_51_method5 = rep(0,B)

res_t_pdf_n_12_method5 = rep(0,B)
res_t_pdf_n_22_method5 = rep(0,B)
res_t_pdf_n_32_method5 = rep(0,B)
res_t_pdf_n_42_method5 = rep(0,B)
res_t_pdf_n_52_method5 = rep(0,B)

# Method-2 #

res_t_pdf_n_12_method2 = rep(0,B)
res_t_pdf_n_22_method2 = rep(0,B)
res_t_pdf_n_32_method2 = rep(0,B)
res_t_pdf_n_42_method2 = rep(0,B)
res_t_pdf_n_52_method2 = rep(0,B)

# Method-3 #

res_t_pdf_n_12_method3 = rep(0,B)
res_t_pdf_n_22_method3 = rep(0,B)
res_t_pdf_n_32_method3 = rep(0,B)
res_t_pdf_n_42_method3 = rep(0,B)
res_t_pdf_n_52_method3 = rep(0,B)

# Method-4 #

res_t_pdf_n_12_method4 = rep(0,B)
res_t_pdf_n_22_method4 = rep(0,B)
res_t_pdf_n_32_method4 = rep(0,B)
res_t_pdf_n_42_method4 = rep(0,B)
res_t_pdf_n_52_method4 = rep(0,B)

set.seed(10)


for (i in 1:B) {
  
  # t Distribution - n1
  
  x = (rt(n[1],df = df))
  t_grid1 <- make_breaks_cover(x, t_grid1)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_11_method1[i] = L1_Distance_calc(x,t_grid,pdf_t)
  err_malc = L2_Distance_calc_both(x,t_grid1,pdf_t)
  res_t_pdf_n_12_method1[i] = err_malc[1]
  smoothed_res_t_pdf_n_12_method1[i] = err_malc[2]
  res_t_pdf_n_12_method2[i] = L2_from_res(x,t_grid1,pdf_t)
  res_t_pdf_n_12_method3[i] = L2_binnednp(x,t_grid1,pdf_t)
  res_t_pdf_n_12_method4[i] = L2_from_kernsmooth(x,t_grid1,pdf_t)
  
  # t Distribution - n2
  
  x = (rt(n[2],df = df))
  t_grid2 <- make_breaks_cover(x, t_grid2)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_21_method1[i] = L1_Distance_calc(x,t_grid,pdf_t)
  err_malc = L2_Distance_calc_both(x,t_grid2,pdf_t)
  res_t_pdf_n_22_method1[i] = err_malc[1]
  smoothed_res_t_pdf_n_22_method1[i] = err_malc[2]
  res_t_pdf_n_22_method2[i] = L2_from_res(x,t_grid2,pdf_t)
  res_t_pdf_n_22_method3[i] = L2_binnednp(x,t_grid2,pdf_t)
  res_t_pdf_n_22_method4[i] = L2_from_kernsmooth(x,t_grid2,pdf_t)
  
  # t Distribution - n3
  
  x = (rt(n[3],df = df))
  t_grid3 <- make_breaks_cover(x, t_grid3)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_31_method1[i] = L1_Distance_calc(x,t_grid,pdf_t)
  err_malc = L2_Distance_calc_both(x,t_grid3,pdf_t)
  res_t_pdf_n_32_method1[i] = err_malc[1]
  smoothed_res_t_pdf_n_32_method1[i] = err_malc[2]
  res_t_pdf_n_32_method2[i] = L2_from_res(x,t_grid3,pdf_t)
  res_t_pdf_n_32_method3[i] = L2_binnednp(x,t_grid3,pdf_t)
  res_t_pdf_n_32_method4[i] = L2_from_kernsmooth(x,t_grid3,pdf_t)
  
  # t Distribution - n4
  
  x = (rt(n[4],df = df))
  t_grid4 <- make_breaks_cover(x, t_grid4)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_41_method1[i] = L1_Distance_calc(x,t_grid,pdf_t)
  err_malc = L2_Distance_calc_both(x,t_grid4,pdf_t)
  res_t_pdf_n_42_method1[i] = err_malc[1]
  smoothed_res_t_pdf_n_42_method1[i] = err_malc[2]
  res_t_pdf_n_42_method2[i] = L2_from_res(x,t_grid4,pdf_t)
  res_t_pdf_n_42_method3[i] = L2_binnednp(x,t_grid4,pdf_t)
  res_t_pdf_n_42_method4[i] = L2_from_kernsmooth(x,t_grid4,pdf_t)
  
  # t Distribution - n5
  
  x = (rt(n[5],df = df))
  t_grid5 <- make_breaks_cover(x, t_grid5)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_51_method1[i] = L1_Distance_calc(x,t_grid,pdf_t)
  err_malc = L2_Distance_calc_both(x,t_grid5,pdf_t)
  res_t_pdf_n_52_method1[i] = err_malc[1]
  smoothed_res_t_pdf_n_52_method1[i] = err_malc[2]
  res_t_pdf_n_52_method2[i] = L2_from_res(x,t_grid5,pdf_t)
  res_t_pdf_n_52_method3[i] = L2_binnednp(x,t_grid5,pdf_t)
  res_t_pdf_n_52_method4[i] = L2_from_kernsmooth(x,t_grid5,pdf_t)
  
  print(i)
  
}


# Method-1 #

boxplot(res_t_pdf_n_12_method1,
        res_t_pdf_n_22_method1,
        res_t_pdf_n_32_method1,
        res_t_pdf_n_42_method1,
        res_t_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.2))


boxplot(smoothed_res_t_pdf_n_12_method1,
        smoothed_res_t_pdf_n_22_method1,
        smoothed_res_t_pdf_n_32_method1,
        smoothed_res_t_pdf_n_42_method1,
        smoothed_res_t_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
dev.off()

boxplot(res_t_pdf_n_12_method2,res_t_pdf_n_22_method2,res_t_pdf_n_32_method2,res_t_pdf_n_42_method2,res_t_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))

boxplot(res_t_pdf_n_12_method3,res_t_pdf_n_22_method3,res_t_pdf_n_32_method3,res_t_pdf_n_42_method3,res_t_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))

boxplot(res_t_pdf_n_12_method4,res_t_pdf_n_22_method4,res_t_pdf_n_32_method4,res_t_pdf_n_42_method4,res_t_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))

n = c(10^3)
B = 100

# Normal Distribution - Grid-Width Simulation

MM = qnorm(0.9999999,mu,sigma) + 10
norm_grids = list(seq(-MM,MM, by = norm_deltas[5]),
                  seq(-MM,MM, by = norm_deltas[4]),
                  seq(-MM,MM, by = norm_deltas[3]),
                  seq(-MM,MM, by = norm_deltas[2]),
                  seq(-MM,MM, by = norm_deltas[1]))

# Method-1 #

res_norm_pdf_gw_11_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_11_method1 = rep(0,B)
res_norm_pdf_gw_21_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_21_method1 = rep(0,B)
res_norm_pdf_gw_31_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_31_method1 = rep(0,B)
res_norm_pdf_gw_41_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_41_method1 = rep(0,B)
res_norm_pdf_gw_51_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_51_method1 = rep(0,B)

res_norm_pdf_gw_12_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_12_method1 = rep(0,B)
res_norm_pdf_gw_22_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_22_method1 = rep(0,B)
res_norm_pdf_gw_32_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_32_method1 = rep(0,B)
res_norm_pdf_gw_42_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_42_method1 = rep(0,B)
res_norm_pdf_gw_52_method1 = rep(0,B)
smoothed_res_norm_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_norm_pdf_gw_12_method2 = rep(0,B)
res_norm_pdf_gw_22_method2 = rep(0,B)
res_norm_pdf_gw_32_method2 = rep(0,B)
res_norm_pdf_gw_42_method2 = rep(0,B)
res_norm_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_norm_pdf_gw_12_method3 = rep(0,B)
res_norm_pdf_gw_22_method3 = rep(0,B)
res_norm_pdf_gw_32_method3 = rep(0,B)
res_norm_pdf_gw_42_method3 = rep(0,B)
res_norm_pdf_gw_52_method3 = rep(0,B)

# Method-4 #

res_norm_pdf_gw_12_method4 = rep(0,B)
res_norm_pdf_gw_22_method4 = rep(0,B)
res_norm_pdf_gw_32_method4 = rep(0,B)
res_norm_pdf_gw_42_method4 = rep(0,B)
res_norm_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  
  # Normal Distribution - GW1
  
  # res_norm_pdf_gw_11_method1[i] = L1_Distance_calc(x,norm_grids[[1]],pdf_norm)
  err_malc = L2_Distance_calc_both(x,norm_grids[[1]],pdf_norm)
  res_norm_pdf_gw_12_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_gw_12_method1[i] = err_malc[2]
  res_norm_pdf_gw_12_method2[i] = L2_from_res(x,norm_grids[[1]],pdf_norm)
  res_norm_pdf_gw_12_method3[i] = L2_binnednp(x,norm_grids[[1]],pdf_norm)
  res_norm_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,norm_grids[[1]],pdf_norm)
  
  # Normal Distribution - GW2
  
  # res_norm_pdf_gw_21_method1[i] = L1_Distance_calc(x,norm_grids[[2]],pdf_norm)
  err_malc = L2_Distance_calc_both(x,norm_grids[[2]],pdf_norm)
  res_norm_pdf_gw_22_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_gw_22_method1[i] = err_malc[2]
  res_norm_pdf_gw_22_method2[i] = L2_from_res(x,norm_grids[[2]],pdf_norm)
  res_norm_pdf_gw_22_method3[i] = L2_binnednp(x,norm_grids[[2]],pdf_norm)
  res_norm_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,norm_grids[[2]],pdf_norm)
  
  # Normal Distribution - GW3
  
  # res_norm_pdf_gw_31_method1[i] = L1_Distance_calc(x,norm_grids[[3]],pdf_norm)
  err_malc = L2_Distance_calc_both(x,norm_grids[[3]],pdf_norm)
  res_norm_pdf_gw_32_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_gw_32_method1[i] = err_malc[2]
  res_norm_pdf_gw_32_method2[i] = L2_from_res(x,norm_grids[[3]],pdf_norm)
  res_norm_pdf_gw_32_method3[i] = L2_binnednp(x,norm_grids[[3]],pdf_norm)
  res_norm_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,norm_grids[[3]],pdf_norm)
  
  # Normal Distribution - GW4
  
  # res_norm_pdf_gw_41_method1[i] = L1_Distance_calc(x,norm_grids[[4]],pdf_norm)
  err_malc = L2_Distance_calc_both(x,norm_grids[[4]],pdf_norm)
  res_norm_pdf_gw_42_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_gw_42_method1[i] = err_malc[2]
  res_norm_pdf_gw_42_method2[i] = L2_from_res(x,norm_grids[[4]],pdf_norm)
  res_norm_pdf_gw_42_method3[i] = L2_binnednp(x,norm_grids[[4]],pdf_norm)
  res_norm_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,norm_grids[[4]],pdf_norm)
  
  # Normal Distribution - GW5
  
  # res_norm_pdf_gw_51_method1[i] = L1_Distance_calc(x,norm_grids[[5]],pdf_norm)
  err_malc = L2_Distance_calc_both(x,norm_grids[[5]],pdf_norm)
  res_norm_pdf_gw_52_method1[i] = err_malc[1]
  smoothed_res_norm_pdf_gw_52_method1[i] = err_malc[2]
  res_norm_pdf_gw_52_method2[i] = L2_from_res(x,norm_grids[[5]],pdf_norm)
  res_norm_pdf_gw_52_method3[i] = L2_binnednp(x,norm_grids[[5]],pdf_norm)
  res_norm_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,norm_grids[[5]],pdf_norm)
  
  print(i)
  
}

# Method-1 #

boxplot(res_norm_pdf_gw_11_method1,res_norm_pdf_gw_21_method1,res_norm_pdf_gw_31_method1,res_norm_pdf_gw_41_method1,res_norm_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))

boxplot(res_norm_pdf_gw_12_method1,res_norm_pdf_gw_22_method1,res_norm_pdf_gw_32_method1,res_norm_pdf_gw_42_method1,res_norm_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))

boxplot(res_norm_pdf_gw_12_method2,res_norm_pdf_gw_22_method2,res_norm_pdf_gw_32_method2,res_norm_pdf_gw_42_method2,res_norm_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))

boxplot(res_norm_pdf_gw_12_method3,res_norm_pdf_gw_22_method3,res_norm_pdf_gw_32_method3,res_norm_pdf_gw_42_method3,res_norm_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))

boxplot(res_norm_pdf_gw_12_method4,res_norm_pdf_gw_22_method4,res_norm_pdf_gw_32_method4,res_norm_pdf_gw_42_method4,res_norm_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))

# Beta Distribution - Grid Width Simulation

beta_grids = list(seq(0,1.2, by = beta_deltas[5]),
                  seq(0,1.3, by = beta_deltas[4]),
                  seq(0,1.2, by = beta_deltas[3]),
                  seq(0,1.2, by = beta_deltas[2]),
                  seq(0,1.1, by = beta_deltas[1]))

# Method-1 #

res_beta_pdf_gw_11_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_11_method1 = rep(0,B)
res_beta_pdf_gw_21_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_21_method1 = rep(0,B)
res_beta_pdf_gw_31_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_31_method1 = rep(0,B)
res_beta_pdf_gw_41_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_41_method1 = rep(0,B)
res_beta_pdf_gw_51_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_51_method1 = rep(0,B)

res_beta_pdf_gw_12_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_12_method1 = rep(0,B)
res_beta_pdf_gw_22_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_22_method1 = rep(0,B)
res_beta_pdf_gw_32_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_32_method1 = rep(0,B)
res_beta_pdf_gw_42_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_42_method1 = rep(0,B)
res_beta_pdf_gw_52_method1 = rep(0,B)
smoothed_res_beta_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_beta_pdf_gw_12_method2 = rep(0,B)
res_beta_pdf_gw_22_method2 = rep(0,B)
res_beta_pdf_gw_32_method2 = rep(0,B)
res_beta_pdf_gw_42_method2 = rep(0,B)
res_beta_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_beta_pdf_gw_12_method3 = rep(0,B)
res_beta_pdf_gw_22_method3 = rep(0,B)
res_beta_pdf_gw_32_method3 = rep(0,B)
res_beta_pdf_gw_42_method3 = rep(0,B)
res_beta_pdf_gw_52_method3 = rep(0,B)

# Method-4 #

res_beta_pdf_gw_12_method4 = rep(0,B)
res_beta_pdf_gw_22_method4 = rep(0,B)
res_beta_pdf_gw_32_method4 = rep(0,B)
res_beta_pdf_gw_42_method4 = rep(0,B)
res_beta_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x <- sort(rbeta(n,alpha,beta))
  
  # Beta Distribution - GW1
  
  # res_beta_pdf_gw_11_method1[i] = L1_Distance_calc(x,beta_grids[[1]],pdf_beta)
  err_malc = L2_Distance_calc_both(x,beta_grids[[1]],pdf_beta)
  res_beta_pdf_gw_12_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_gw_12_method1[i] = err_malc[2]
  res_beta_pdf_gw_12_method2[i] = L2_from_res(x,beta_grids[[1]],pdf_beta)
  res_beta_pdf_gw_12_method3[i] = L2_binnednp(x,beta_grids[[1]],pdf_beta)
  res_beta_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,beta_grids[[1]],pdf_beta)
  
  
  # Beta Distribution - GW2
  
  # res_beta_pdf_gw_21_method1[i] = L1_Distance_calc(x,beta_grids[[2]],pdf_beta)
  err_malc = L2_Distance_calc_both(x,beta_grids[[2]],pdf_beta)
  res_beta_pdf_gw_22_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_gw_22_method1[i] = err_malc[2]
  res_beta_pdf_gw_22_method2[i] = L2_from_res(x,beta_grids[[2]],pdf_beta)
  res_beta_pdf_gw_22_method3[i] = L2_binnednp(x,beta_grids[[2]],pdf_beta)
  res_beta_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,beta_grids[[2]],pdf_beta)
  
  # Beta Distribution - GW3
  
  # res_beta_pdf_gw_31_method1[i] = L1_Distance_calc(x,beta_grids[[3]],pdf_beta)
  err_malc = L2_Distance_calc_both(x,beta_grids[[3]],pdf_beta)
  res_beta_pdf_gw_32_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_gw_32_method1[i] = err_malc[2]
  res_beta_pdf_gw_32_method2[i] = L2_from_res(x,beta_grids[[3]],pdf_beta)
  res_beta_pdf_gw_32_method3[i] = L2_binnednp(x,beta_grids[[3]],pdf_beta)
  res_beta_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,beta_grids[[3]],pdf_beta)
  
  # Beta Distribution - GW4
  
  # res_beta_pdf_gw_41_method1[i] = L1_Distance_calc(x,beta_grids[[4]],pdf_beta)
  err_malc = L2_Distance_calc_both(x,beta_grids[[4]],pdf_beta)
  res_beta_pdf_gw_42_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_gw_42_method1[i] = err_malc[2]
  res_beta_pdf_gw_42_method2[i] = L2_from_res(x,beta_grids[[4]],pdf_beta)
  res_beta_pdf_gw_42_method3[i] = L2_binnednp(x,beta_grids[[4]],pdf_beta)
  res_beta_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,beta_grids[[4]],pdf_beta)
  
  
  # Beta Distribution - GW5
  
  # res_beta_pdf_gw_51_method1[i] = L1_Distance_calc(x,beta_grids[[5]],pdf_beta)
  err_malc = L2_Distance_calc_both(x,beta_grids[[5]],pdf_beta)
  res_beta_pdf_gw_52_method1[i] = err_malc[1]
  smoothed_res_beta_pdf_gw_52_method1[i] = err_malc[2]
  res_beta_pdf_gw_52_method2[i] = L2_from_res(x,beta_grids[[5]],pdf_beta)
  res_beta_pdf_gw_52_method3[i] = L2_binnednp(x,beta_grids[[5]],pdf_beta)
  res_beta_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,beta_grids[[5]],pdf_beta)
  
  print(i)
  
}

# Method-1 #

boxplot(res_beta_pdf_gw_11_method1,res_beta_pdf_gw_21_method1,res_beta_pdf_gw_31_method1,res_beta_pdf_gw_41_method1,res_beta_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))

boxplot(res_beta_pdf_gw_12_method1,res_beta_pdf_gw_22_method1,res_beta_pdf_gw_32_method1,res_beta_pdf_gw_42_method1,res_beta_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))

boxplot(res_beta_pdf_gw_12_method2,res_beta_pdf_gw_22_method2,res_beta_pdf_gw_32_method2,res_beta_pdf_gw_42_method2,res_beta_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))

boxplot(res_beta_pdf_gw_12_method3,res_beta_pdf_gw_22_method3,res_beta_pdf_gw_32_method3,res_beta_pdf_gw_42_method3,res_beta_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))

boxplot(res_beta_pdf_gw_12_method4,res_beta_pdf_gw_22_method4,res_beta_pdf_gw_32_method4,res_beta_pdf_gw_42_method4,res_beta_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))


# Gamma Distribution - Grid Width Simulation

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grids = list(seq(0,MM, by = gamma_deltas[5]),
                   seq(0,MM, by = gamma_deltas[4]),
                   seq(0,MM, by = gamma_deltas[3]),
                   seq(0,MM, by = gamma_deltas[2]),
                   seq(0,MM, by = gamma_deltas[1]))

# Method-1 #

res_gamma_pdf_gw_11_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_11_method1 = rep(0,B)
res_gamma_pdf_gw_21_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_21_method1 = rep(0,B)
res_gamma_pdf_gw_31_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_31_method1 = rep(0,B)
res_gamma_pdf_gw_41_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_41_method1 = rep(0,B)
res_gamma_pdf_gw_51_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_51_method1 = rep(0,B)

res_gamma_pdf_gw_12_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_12_method1 = rep(0,B)
res_gamma_pdf_gw_22_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_22_method1 = rep(0,B)
res_gamma_pdf_gw_32_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_32_method1 = rep(0,B)
res_gamma_pdf_gw_42_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_42_method1 = rep(0,B)
res_gamma_pdf_gw_52_method1 = rep(0,B)
smoothed_res_gamma_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_gamma_pdf_gw_12_method2 = rep(0,B)
res_gamma_pdf_gw_22_method2 = rep(0,B)
res_gamma_pdf_gw_32_method2 = rep(0,B)
res_gamma_pdf_gw_42_method2 = rep(0,B)
res_gamma_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_gamma_pdf_gw_12_method3 = rep(0,B)
res_gamma_pdf_gw_22_method3 = rep(0,B)
res_gamma_pdf_gw_32_method3 = rep(0,B)
res_gamma_pdf_gw_42_method3 = rep(0,B)
res_gamma_pdf_gw_52_method3 = rep(0,B)

# Method-4 #

res_gamma_pdf_gw_12_method4 = rep(0,B)
res_gamma_pdf_gw_22_method4 = rep(0,B)
res_gamma_pdf_gw_32_method4 = rep(0,B)
res_gamma_pdf_gw_42_method4 = rep(0,B)
res_gamma_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x <- sort(rgamma(n,alpha,beta))
  
  # Gamma Distribution - GW1
  
  # res_gamma_pdf_gw_11_method1[i] = L1_Distance_calc(x,gamma_grids[[1]],pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grids[[1]],pdf_gamma)
  res_gamma_pdf_gw_12_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_gw_12_method1[i] = err_malc[2]
  res_gamma_pdf_gw_12_method2[i] = L2_from_res(x,gamma_grids[[1]],pdf_gamma)
  res_gamma_pdf_gw_12_method3[i] = L2_binnednp(x,gamma_grids[[1]],pdf_gamma)
  res_gamma_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,gamma_grids[[1]],pdf_gamma)
  
  # Gamma Distribution - GW2
  
  # res_gamma_pdf_gw_21_method1[i] = L1_Distance_calc(x,gamma_grids[[2]],pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grids[[2]],pdf_gamma)
  res_gamma_pdf_gw_22_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_gw_22_method1[i] = err_malc[2]
  res_gamma_pdf_gw_22_method2[i] = L2_from_res(x,gamma_grids[[2]],pdf_gamma)
  res_gamma_pdf_gw_22_method3[i] = L2_binnednp(x,gamma_grids[[2]],pdf_gamma)
  res_gamma_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,gamma_grids[[2]],pdf_gamma)
  
  # Gamma Distribution - GW3
  
  # res_gamma_pdf_gw_31_method1[i] = L1_Distance_calc(x,gamma_grids[[3]],pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grids[[3]],pdf_gamma)
  res_gamma_pdf_gw_32_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_gw_32_method1[i] = err_malc[2]
  res_gamma_pdf_gw_32_method2[i] = L2_from_res(x,gamma_grids[[3]],pdf_gamma)
  res_gamma_pdf_gw_32_method3[i] = L2_binnednp(x,gamma_grids[[3]],pdf_gamma)
  res_gamma_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,gamma_grids[[3]],pdf_gamma)
  
  # Gamma Distribution - GW4
  
  # res_gamma_pdf_gw_41_method1[i] = L1_Distance_calc(x,gamma_grids[[4]],pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grids[[4]],pdf_gamma)
  res_gamma_pdf_gw_42_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_gw_42_method1[i] = err_malc[2]
  res_gamma_pdf_gw_42_method2[i] = L2_from_res(x,gamma_grids[[4]],pdf_gamma)
  res_gamma_pdf_gw_42_method3[i] = L2_binnednp(x,gamma_grids[[4]],pdf_gamma)
  res_gamma_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,gamma_grids[[4]],pdf_gamma)
  
  # Gamma Distribution - GW5
  
  # res_gamma_pdf_gw_51_method1[i] = L1_Distance_calc(x,gamma_grids[[5]],pdf_gamma)
  err_malc = L2_Distance_calc_both(x,gamma_grids[[5]],pdf_gamma)
  res_gamma_pdf_gw_52_method1[i] = err_malc[1]
  smoothed_res_gamma_pdf_gw_52_method1[i] = err_malc[2]
  res_gamma_pdf_gw_52_method2[i] = L2_from_res(x,gamma_grids[[5]],pdf_gamma)
  res_gamma_pdf_gw_52_method3[i] = L2_binnednp(x,gamma_grids[[5]],pdf_gamma)
  res_gamma_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,gamma_grids[[5]],pdf_gamma)
  
  print(i)
  
}

# Method #

boxplot(res_gamma_pdf_gw_11_method1,res_gamma_pdf_gw_21_method1,res_gamma_pdf_gw_31_method1,res_gamma_pdf_gw_41_method1,res_gamma_pdf_gw_51_method1,
        names = c("1 (7)", "0.8 (8)", "0.6 (11)", "0.4 (16)", "0.2 (31)"),
        col = colors3,ylim=c(0,0.9))

boxplot(res_gamma_pdf_gw_12_method1,res_gamma_pdf_gw_22_method1,res_gamma_pdf_gw_32_method1,res_gamma_pdf_gw_42_method1,res_gamma_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_gw_12_method2,res_gamma_pdf_gw_22_method2,res_gamma_pdf_gw_32_method2,res_gamma_pdf_gw_42_method2,res_gamma_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_gw_12_method3,res_gamma_pdf_gw_22_method3,res_gamma_pdf_gw_32_method3,res_gamma_pdf_gw_42_method3,res_gamma_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))

boxplot(res_gamma_pdf_gw_12_method4,res_gamma_pdf_gw_22_method4,res_gamma_pdf_gw_32_method4,res_gamma_pdf_gw_42_method4,res_gamma_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))


# Logistic Distribution - Grid Width Simulation

MM = qlogis(0.9999999,mu,sigma) + 3
logistic_grids = list(seq(-MM,MM, by = logistic_deltas[5]),
                      seq(-MM,MM, by = logistic_deltas[4]),
                      seq(-MM,MM, by = logistic_deltas[3]),
                      seq(-MM,MM, by = logistic_deltas[2]),
                      seq(-MM,MM, by = logistic_deltas[1]))


# Method-1 #

res_logis_pdf_gw_11_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_11_method1 = rep(0,B)
res_logis_pdf_gw_21_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_21_method1 = rep(0,B)
res_logis_pdf_gw_31_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_31_method1 = rep(0,B)
res_logis_pdf_gw_41_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_41_method1 = rep(0,B)
res_logis_pdf_gw_51_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_51_method1 = rep(0,B)

res_logis_pdf_gw_12_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_12_method1 = rep(0,B)
res_logis_pdf_gw_22_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_22_method1 = rep(0,B)
res_logis_pdf_gw_32_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_32_method1 = rep(0,B)
res_logis_pdf_gw_42_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_42_method1 = rep(0,B)
res_logis_pdf_gw_52_method1 = rep(0,B)
smoothed_res_logis_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_logis_pdf_gw_12_method2 = rep(0,B)
res_logis_pdf_gw_22_method2 = rep(0,B)
res_logis_pdf_gw_32_method2 = rep(0,B)
res_logis_pdf_gw_42_method2 = rep(0,B)
res_logis_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_logis_pdf_gw_12_method3 = rep(0,B)
res_logis_pdf_gw_22_method3 = rep(0,B)
res_logis_pdf_gw_32_method3 = rep(0,B)
res_logis_pdf_gw_42_method3 = rep(0,B)
res_logis_pdf_gw_52_method3 = rep(0,B)

# Method-4 #

res_logis_pdf_gw_12_method4 = rep(0,B)
res_logis_pdf_gw_22_method4 = rep(0,B)
res_logis_pdf_gw_32_method4 = rep(0,B)
res_logis_pdf_gw_42_method4 = rep(0,B)
res_logis_pdf_gw_52_method4 = rep(0,B)

library(R.utils)

safe_L2_Distance_calc_both <- function(x, grid, pdf, timeout = 120) {
  withTimeout(
    L2_Distance_calc_both(x, grid, pdf),
    timeout = timeout,
    onTimeout = "error"
  )
}

set.seed(10)

i = 1

while (i<=B) {
  
  x = sort(rlogis(n))
  
  # Logistic Distribution - GW1
  
  # res_logis_pdf_gw_11_method1[i] = L1_Distance_calc(x,logistic_grids[[1]],pdf_logistic)
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, logistic_grids[[1]], pdf_logistic),
                         timeout = 120, onTimeout = "error"),
    silent = TRUE
  )
  
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  
  res_logis_pdf_gw_12_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_gw_12_method1[i] = err_malc[2]
  res_logis_pdf_gw_12_method2[i] = L2_from_res(x,logistic_grids[[1]],pdf_logistic)
  res_logis_pdf_gw_12_method3[i] = L2_binnednp(x,logistic_grids[[1]],pdf_logistic)
  res_logis_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,logistic_grids[[1]],pdf_logistic)
  
  # Logistic Distribution - GW2
  
  # res_logis_pdf_gw_21_method1[i] = L1_Distance_calc(x,logistic_grids[[2]],pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grids[[2]],pdf_logistic)
  res_logis_pdf_gw_22_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_gw_22_method1[i] = err_malc[2]
  res_logis_pdf_gw_22_method2[i] = L2_from_res(x,logistic_grids[[2]],pdf_logistic)
  res_logis_pdf_gw_22_method3[i] = L2_binnednp(x,logistic_grids[[2]],pdf_logistic)
  res_logis_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,logistic_grids[[2]],pdf_logistic)
  
  # Logistic Distribution - GW3
  
  # res_logis_pdf_gw_31_method1[i] = L1_Distance_calc(x,logistic_grids[[3]],pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grids[[3]],pdf_logistic)
  res_logis_pdf_gw_32_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_gw_32_method1[i] = err_malc[2]
  res_logis_pdf_gw_32_method2[i] = L2_from_res(x,logistic_grids[[3]],pdf_logistic)
  res_logis_pdf_gw_32_method3[i] = L2_binnednp(x,logistic_grids[[3]],pdf_logistic)
  res_logis_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,logistic_grids[[3]],pdf_logistic)
  
  # Logistic Distribution - GW4
  
  # res_logis_pdf_gw_41_method1[i] = L1_Distance_calc(x,logistic_grids[[4]],pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grids[[4]],pdf_logistic)
  res_logis_pdf_gw_42_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_gw_42_method1[i] = err_malc[2]
  res_logis_pdf_gw_42_method2[i] = L2_from_res(x,logistic_grids[[4]],pdf_logistic)
  res_logis_pdf_gw_42_method3[i] = L2_binnednp(x,logistic_grids[[4]],pdf_logistic)
  res_logis_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,logistic_grids[[4]],pdf_logistic)
  
  # Logistic Distribution - GW5
  
  # res_logis_pdf_gw_51_method1[i] = L1_Distance_calc(x,logistic_grids[[5]],pdf_logistic)
  err_malc = L2_Distance_calc_both(x,logistic_grids[[5]],pdf_logistic)
  res_logis_pdf_gw_52_method1[i] = err_malc[1]
  smoothed_res_logis_pdf_gw_52_method1[i] = err_malc[2]
  res_logis_pdf_gw_52_method2[i] = L2_from_res(x,logistic_grids[[5]],pdf_logistic)
  res_logis_pdf_gw_52_method3[i] = L2_binnednp(x,logistic_grids[[5]],pdf_logistic)
  res_logis_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,logistic_grids[[5]],pdf_logistic)
  
  i = i + 1
  print(i)
  
}

# Method #

boxplot(res_logis_pdf_gw_11_method1,res_logis_pdf_gw_21_method1,res_logis_pdf_gw_31_method1,res_logis_pdf_gw_41_method1,res_logis_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))

boxplot(res_logis_pdf_gw_12_method1,res_logis_pdf_gw_22_method1,res_logis_pdf_gw_32_method1,res_logis_pdf_gw_42_method1,res_logis_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_logis_pdf_gw_12_method2,res_logis_pdf_gw_22_method2,res_logis_pdf_gw_32_method2,res_logis_pdf_gw_42_method2,res_logis_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_logis_pdf_gw_12_method3,res_logis_pdf_gw_22_method3,res_logis_pdf_gw_32_method3,res_logis_pdf_gw_42_method3,res_logis_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

boxplot(res_logis_pdf_gw_12_method4,res_logis_pdf_gw_22_method4,res_logis_pdf_gw_32_method4,res_logis_pdf_gw_42_method4,res_logis_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))

# Student's t Distribution - Grid Width Simulation

MM = qt(0.9999999,df) + 3
t_grids = list(seq(-MM,MM, by = t_deltas[5]),
               seq(-MM,MM, by = t_deltas[4]),
               seq(-MM,MM, by = t_deltas[3]),
               seq(-MM,MM, by = t_deltas[2]),
               seq(-MM,MM, by = t_deltas[1]))

# Method-1 #

res_t_pdf_gw_11_method1 = rep(0,B)
smoothed_res_t_pdf_gw_11_method1 = rep(0,B)
res_t_pdf_gw_21_method1 = rep(0,B)
smoothed_res_t_pdf_gw_21_method1 = rep(0,B)
res_t_pdf_gw_31_method1 = rep(0,B)
smoothed_res_t_pdf_gw_31_method1 = rep(0,B)
res_t_pdf_gw_41_method1 = rep(0,B)
smoothed_res_t_pdf_gw_41_method1 = rep(0,B)
res_t_pdf_gw_51_method1 = rep(0,B)
smoothed_res_t_pdf_gw_51_method1 = rep(0,B)

res_t_pdf_gw_12_method1 = rep(0,B)
smoothed_res_t_pdf_gw_12_method1 = rep(0,B)
res_t_pdf_gw_22_method1 = rep(0,B)
smoothed_res_t_pdf_gw_22_method1 = rep(0,B)
res_t_pdf_gw_32_method1 = rep(0,B)
smoothed_res_t_pdf_gw_32_method1 = rep(0,B)
res_t_pdf_gw_42_method1 = rep(0,B)
smoothed_res_t_pdf_gw_42_method1 = rep(0,B)
res_t_pdf_gw_52_method1 = rep(0,B)
smoothed_res_t_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_t_pdf_gw_12_method2 = rep(0,B)
res_t_pdf_gw_22_method2 = rep(0,B)
res_t_pdf_gw_32_method2 = rep(0,B)
res_t_pdf_gw_42_method2 = rep(0,B)
res_t_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_t_pdf_gw_12_method3 = rep(0,B)
res_t_pdf_gw_22_method3 = rep(0,B)
res_t_pdf_gw_32_method3 = rep(0,B)
res_t_pdf_gw_42_method3 = rep(0,B)
res_t_pdf_gw_52_method3 = rep(0,B)

# Method-4 #

res_t_pdf_gw_12_method4 = rep(0,B)
res_t_pdf_gw_22_method4 = rep(0,B)
res_t_pdf_gw_32_method4 = rep(0,B)
res_t_pdf_gw_42_method4 = rep(0,B)
res_t_pdf_gw_52_method4 = rep(0,B)

set.seed(10)

i = 1 

while (i <= B) {
  
  x = sort(rt(n,df = df))
  # t_grids[[1]] <- make_breaks_cover(x, t_grids[[1]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[2]] <- make_breaks_cover(x, t_grids[[2]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[3]] <- make_breaks_cover(x, t_grids[[3]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[4]] <- make_breaks_cover(x, t_grids[[4]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[5]] <- make_breaks_cover(x, t_grids[[5]])   # use the matching grid (t_grid1...t_grid5)
  
  # t Distribution - GW1
  
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, t_grids[[1]], pdf_t),
                         timeout = 60, onTimeout = "error"),
    silent = TRUE
  )
  
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  
  res_t_pdf_gw_12_method1[i] = err_malc[1]
  smoothed_res_t_pdf_gw_12_method1[i] = err_malc[2]
  res_t_pdf_gw_12_method2[i] = L2_from_res(x,t_grids[[1]],pdf_t)
  res_t_pdf_gw_12_method3[i] = L2_binnednp(x,t_grids[[1]],pdf_t)
  res_t_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,t_grids[[1]],pdf_t)
  
  # t Distribution - GW2
  
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, t_grids[[2]], pdf_t),
                         timeout = 60, onTimeout = "error"),
    silent = TRUE
  )
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  
  res_t_pdf_gw_22_method1[i] = err_malc[1]
  smoothed_res_t_pdf_gw_22_method1[i] = err_malc[2]
  res_t_pdf_gw_22_method2[i] = L2_from_res(x,t_grids[[2]],pdf_t)
  res_t_pdf_gw_22_method3[i] = L2_binnednp(x,t_grids[[2]],pdf_t)
  res_t_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,t_grids[[2]],pdf_t)
  
  
  # t Distribution - GW3
  
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, t_grids[[3]], pdf_t),
                         timeout = 60, onTimeout = "error"),
    silent = TRUE
  )
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  
  res_t_pdf_gw_32_method1[i] = err_malc[1]
  smoothed_res_t_pdf_gw_32_method1[i] = err_malc[2]
  res_t_pdf_gw_32_method2[i] = L2_from_res(x,t_grids[[3]],pdf_t)
  res_t_pdf_gw_32_method3[i] = L2_binnednp(x,t_grids[[3]],pdf_t)
  res_t_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,t_grids[[3]],pdf_t)
  
  # t Distribution - GW4
  
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, t_grids[[4]], pdf_t),
                         timeout = 60, onTimeout = "error"),
    silent = TRUE
  )
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  
  res_t_pdf_gw_42_method1[i] = err_malc[1]
  smoothed_res_t_pdf_gw_42_method1[i] = err_malc[2]
  res_t_pdf_gw_42_method2[i] = L2_from_res(x,t_grids[[4]],pdf_t)
  res_t_pdf_gw_42_method3[i] = L2_binnednp(x,t_grids[[4]],pdf_t)
  res_t_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,t_grids[[4]],pdf_t)
  
  # t Distribution - GW5
  
  err_malc <- try(
    R.utils::withTimeout(L2_Distance_calc_both(x, t_grids[[5]], pdf_t),
                         timeout = 60, onTimeout = "error"),
    silent = TRUE
  )
  if (inherits(err_malc, "try-error")) next  # draw new x, same i
  res_t_pdf_gw_52_method1[i] = err_malc[1]
  smoothed_res_t_pdf_gw_52_method1[i] = err_malc[2]
  res_t_pdf_gw_52_method2[i] = L2_from_res(x,t_grids[[5]],pdf_t)
  res_t_pdf_gw_52_method3[i] = L2_binnednp(x,t_grids[[5]],pdf_t)
  res_t_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,t_grids[[5]],pdf_t)
  
  print(i)
  
  i = i + 1 
  
}

# Method #

boxplot(res_t_pdf_gw_11_method1,res_t_pdf_gw_21_method1,res_t_pdf_gw_31_method1,res_t_pdf_gw_41_method1,res_t_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5),xlab('Grid Width'))

boxplot(res_t_pdf_gw_12_method1,res_t_pdf_gw_22_method1,res_t_pdf_gw_32_method1,res_t_pdf_gw_42_method1,res_t_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))

boxplot(res_t_pdf_gw_12_method2,res_t_pdf_gw_22_method2,res_t_pdf_gw_32_method2,res_t_pdf_gw_42_method2,res_t_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))

boxplot(res_t_pdf_gw_12_method3,res_t_pdf_gw_22_method3,res_t_pdf_gw_32_method3,res_t_pdf_gw_42_method3,res_t_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))

boxplot(res_t_pdf_gw_12_method4,res_t_pdf_gw_22_method4,res_t_pdf_gw_32_method4,res_t_pdf_gw_42_method4,res_t_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))

# Data Save

# Create the data frame
df_norm_n_method1 <- data.frame(
  `1e2` = res_norm_pdf_n_12_method1,
  `1e3` = res_norm_pdf_n_22_method1,
  `1e4` = res_norm_pdf_n_32_method1,
  `1e5` = res_norm_pdf_n_42_method1,
  `1e6` = res_norm_pdf_n_52_method1
)

df_norm_n_method2 <- data.frame(
  `1e2` = res_norm_pdf_n_12_method2,
  `1e3` = res_norm_pdf_n_22_method2,
  `1e4` = res_norm_pdf_n_32_method2,
  `1e5` = res_norm_pdf_n_42_method2,
  `1e6` = res_norm_pdf_n_52_method2
)

df_norm_n_method3 <- data.frame(
  `1e2` = res_norm_pdf_n_12_method3,
  `1e3` = res_norm_pdf_n_22_method3,
  `1e4` = res_norm_pdf_n_32_method3,
  `1e5` = res_norm_pdf_n_42_method3,
  `1e6` = res_norm_pdf_n_52_method3
)

df_norm_n_method4 <- data.frame(
  `1e2` = res_norm_pdf_n_12_method4,
  `1e3` = res_norm_pdf_n_22_method4,
  `1e4` = res_norm_pdf_n_32_method4,
  `1e5` = res_norm_pdf_n_42_method4,
  `1e6` = res_norm_pdf_n_52_method4
)


# Create the data frame
df_beta_n_method1 <- data.frame(
  `1e2` = res_beta_pdf_n_12_method1,
  `1e3` = res_beta_pdf_n_22_method1,
  `1e4` = res_beta_pdf_n_32_method1,
  `1e5` = res_beta_pdf_n_42_method1,
  `1e6` = res_beta_pdf_n_52_method1
)

df_beta_n_method2 <- data.frame(
  `1e2` = res_beta_pdf_n_12_method2,
  `1e3` = res_beta_pdf_n_22_method2,
  `1e4` = res_beta_pdf_n_32_method2,
  `1e5` = res_beta_pdf_n_42_method2,
  `1e6` = res_beta_pdf_n_52_method2
)

df_beta_n_method3 <- data.frame(
  `1e2` = res_beta_pdf_n_12_method3,
  `1e3` = res_beta_pdf_n_22_method3,
  `1e4` = res_beta_pdf_n_32_method3,
  `1e5` = res_beta_pdf_n_42_method3,
  `1e6` = res_beta_pdf_n_52_method3
)

df_beta_n_method4 <- data.frame(
  `1e2` = res_beta_pdf_n_12_method4,
  `1e3` = res_beta_pdf_n_22_method4,
  `1e4` = res_beta_pdf_n_32_method4,
  `1e5` = res_beta_pdf_n_42_method4,
  `1e6` = res_beta_pdf_n_52_method4
)

# Create the data frame
df_gamma_n_method1 <- data.frame(
  `1e2` = res_gamma_pdf_n_12_method1,
  `1e3` = res_gamma_pdf_n_22_method1,
  `1e4` = res_gamma_pdf_n_32_method1,
  `1e5` = res_gamma_pdf_n_42_method1,
  `1e6` = res_gamma_pdf_n_52_method1
)

df_gamma_n_method2 <- data.frame(
  `1e2` = res_gamma_pdf_n_12_method2,
  `1e3` = res_gamma_pdf_n_22_method2,
  `1e4` = res_gamma_pdf_n_32_method2,
  `1e5` = res_gamma_pdf_n_42_method2,
  `1e6` = res_gamma_pdf_n_52_method2
)

df_gamma_n_method3 <- data.frame(
  `1e2` = res_gamma_pdf_n_12_method3,
  `1e3` = res_gamma_pdf_n_22_method3,
  `1e4` = res_gamma_pdf_n_32_method3,
  `1e5` = res_gamma_pdf_n_42_method3,
  `1e6` = res_gamma_pdf_n_52_method3
)

df_gamma_n_method4 <- data.frame(
  `1e2` = res_gamma_pdf_n_12_method4,
  `1e3` = res_gamma_pdf_n_22_method4,
  `1e4` = res_gamma_pdf_n_32_method4,
  `1e5` = res_gamma_pdf_n_42_method4,
  `1e6` = res_gamma_pdf_n_52_method4
)

# Create the data frame
df_logis_n_method1 <- data.frame(
  `1e2` = res_logis_pdf_n_12_method1,
  `1e3` = res_logis_pdf_n_22_method1,
  `1e4` = res_logis_pdf_n_32_method1,
  `1e5` = res_logis_pdf_n_42_method1,
  `1e6` = res_logis_pdf_n_52_method1
)

df_logis_n_method2 <- data.frame(
  `1e2` = res_logis_pdf_n_12_method2,
  `1e3` = res_logis_pdf_n_22_method2,
  `1e4` = res_logis_pdf_n_32_method2,
  `1e5` = res_logis_pdf_n_42_method2,
  `1e6` = res_logis_pdf_n_52_method2
)

df_logis_n_method3 <- data.frame(
  `1e2` = res_logis_pdf_n_12_method3,
  `1e3` = res_logis_pdf_n_22_method3,
  `1e4` = res_logis_pdf_n_32_method3,
  `1e5` = res_logis_pdf_n_42_method3,
  `1e6` = res_logis_pdf_n_52_method3
)

df_logis_n_method4 <- data.frame(
  `1e2` = res_logis_pdf_n_12_method4,
  `1e3` = res_logis_pdf_n_22_method4,
  `1e4` = res_logis_pdf_n_32_method4,
  `1e5` = res_logis_pdf_n_42_method4,
  `1e6` = res_logis_pdf_n_52_method4
)


# Create the data frame
df_t_n_method1 <- data.frame(
  `1e2` = res_t_pdf_n_12_method1,
  `1e3` = res_t_pdf_n_22_method1,
  `1e4` = res_t_pdf_n_32_method1,
  `1e5` = res_t_pdf_n_42_method1,
  `1e6` = res_t_pdf_n_52_method1
)

df_t_n_method2 <- data.frame(
  `1e2` = res_t_pdf_n_12_method2,
  `1e3` = res_t_pdf_n_22_method2,
  `1e4` = res_t_pdf_n_32_method2,
  `1e5` = res_t_pdf_n_42_method2,
  `1e6` = res_t_pdf_n_52_method2
)

df_t_n_method3 <- data.frame(
  `1e2` = res_t_pdf_n_12_method3,
  `1e3` = res_t_pdf_n_22_method3,
  `1e4` = res_t_pdf_n_32_method3,
  `1e5` = res_t_pdf_n_42_method3,
  `1e6` = res_t_pdf_n_52_method3
)

df_t_n_method4 <- data.frame(
  `1e2` = res_t_pdf_n_12_method4,
  `1e3` = res_t_pdf_n_22_method4,
  `1e4` = res_t_pdf_n_32_method4,
  `1e5` = res_t_pdf_n_42_method4,
  `1e6` = res_t_pdf_n_52_method4
)

df_norm_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_norm_pdf_n_12_method1,
  `1e3` = smoothed_res_norm_pdf_n_22_method1,
  `1e4` = smoothed_res_norm_pdf_n_32_method1,
  `1e5` = smoothed_res_norm_pdf_n_42_method1,
  `1e6` = smoothed_res_norm_pdf_n_52_method1
)

df_beta_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_beta_pdf_n_12_method1,
  `1e3` = smoothed_res_beta_pdf_n_22_method1,
  `1e4` = smoothed_res_beta_pdf_n_32_method1,
  `1e5` = smoothed_res_beta_pdf_n_42_method1,
  `1e6` = smoothed_res_beta_pdf_n_52_method1
)

df_gamma_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_gamma_pdf_n_12_method1,
  `1e3` = smoothed_res_gamma_pdf_n_22_method1,
  `1e4` = smoothed_res_gamma_pdf_n_32_method1,
  `1e5` = smoothed_res_gamma_pdf_n_42_method1,
  `1e6` = smoothed_res_gamma_pdf_n_52_method1
)

df_logis_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_logis_pdf_n_12_method1,
  `1e3` = smoothed_res_logis_pdf_n_22_method1,
  `1e4` = smoothed_res_logis_pdf_n_32_method1,
  `1e5` = smoothed_res_logis_pdf_n_42_method1,
  `1e6` = smoothed_res_logis_pdf_n_52_method1
)

df_t_n_method1_smoothed <- data.frame(
  `1e2` = smoothed_res_t_pdf_n_12_method1,
  `1e3` = smoothed_res_t_pdf_n_22_method1,
  `1e4` = smoothed_res_t_pdf_n_32_method1,
  `1e5` = smoothed_res_t_pdf_n_42_method1,
  `1e6` = smoothed_res_t_pdf_n_52_method1
)

write.csv(df_norm_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_beta_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_gamma_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_logis_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1.csv", row.names = FALSE)
write.csv(df_t_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1.csv", row.names = FALSE)

write.csv(df_norm_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_beta_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_gamma_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_logis_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_t_n_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1_smoothed.csv", row.names = FALSE)

write.csv(df_norm_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_beta_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_gamma_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_logis_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method2.csv", row.names = FALSE)
write.csv(df_t_n_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method2.csv", row.names = FALSE)

write.csv(df_norm_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_beta_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_gamma_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_logis_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method3.csv", row.names = FALSE)
write.csv(df_t_n_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method3.csv", row.names = FALSE)

write.csv(df_norm_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_beta_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_gamma_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_logis_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method4.csv", row.names = FALSE)
write.csv(df_t_n_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method4.csv", row.names = FALSE)

# Create the data frame
df_norm_gw_method1 <- data.frame(
  `2.5` = res_norm_pdf_gw_12_method1,
  `2.0` = res_norm_pdf_gw_22_method1,
  `1.5` = res_norm_pdf_gw_32_method1,
  `1.0` = res_norm_pdf_gw_42_method1,
  `0.5` = res_norm_pdf_gw_52_method1
)

df_norm_gw_method2 <- data.frame(
  `2.5` = res_norm_pdf_gw_12_method2,
  `2.0` = res_norm_pdf_gw_22_method2,
  `1.5` = res_norm_pdf_gw_32_method2,
  `1.0` = res_norm_pdf_gw_42_method2,
  `0.5` = res_norm_pdf_gw_52_method2
)

df_norm_gw_method3 <- data.frame(
  `2.5` = res_norm_pdf_gw_12_method3,
  `2.0` = res_norm_pdf_gw_22_method3,
  `1.5` = res_norm_pdf_gw_32_method3,
  `1.0` = res_norm_pdf_gw_42_method3,
  `0.5` = res_norm_pdf_gw_52_method3
)

df_norm_gw_method4 <- data.frame(
  `2.5` = res_norm_pdf_gw_12_method4,
  `2.0` = res_norm_pdf_gw_22_method4,
  `1.5` = res_norm_pdf_gw_32_method4,
  `1.0` = res_norm_pdf_gw_42_method4,
  `0.5` = res_norm_pdf_gw_52_method4
)

# Create the data frame
df_beta_gw_method1 <- data.frame(
  `2.5` = res_beta_pdf_gw_12_method1,
  `2.0` = res_beta_pdf_gw_22_method1,
  `1.5` = res_beta_pdf_gw_32_method1,
  `1.0` = res_beta_pdf_gw_42_method1,
  `0.5` = res_beta_pdf_gw_52_method1
)

df_beta_gw_method2 <- data.frame(
  `2.5` = res_beta_pdf_gw_12_method2,
  `2.0` = res_beta_pdf_gw_22_method2,
  `1.5` = res_beta_pdf_gw_32_method2,
  `1.0` = res_beta_pdf_gw_42_method2,
  `0.5` = res_beta_pdf_gw_52_method2
)

df_beta_gw_method3 <- data.frame(
  `2.5` = res_beta_pdf_gw_12_method3,
  `2.0` = res_beta_pdf_gw_22_method3,
  `1.5` = res_beta_pdf_gw_32_method3,
  `1.0` = res_beta_pdf_gw_42_method3,
  `0.5` = res_beta_pdf_gw_52_method3
)

df_beta_gw_method4 <- data.frame(
  `2.5` = res_beta_pdf_gw_12_method4,
  `2.0` = res_beta_pdf_gw_22_method4,
  `1.5` = res_beta_pdf_gw_32_method4,
  `1.0` = res_beta_pdf_gw_42_method4,
  `0.5` = res_beta_pdf_gw_52_method4
)

# Create the data frame
df_gamma_gw_method1 <- data.frame(
  `2.5` = res_gamma_pdf_gw_12_method1,
  `2.0` = res_gamma_pdf_gw_22_method1,
  `1.5` = res_gamma_pdf_gw_32_method1,
  `1.0` = res_gamma_pdf_gw_42_method1,
  `0.5` = res_gamma_pdf_gw_52_method1
)

df_gamma_gw_method2 <- data.frame(
  `2.5` = res_gamma_pdf_gw_12_method2,
  `2.0` = res_gamma_pdf_gw_22_method2,
  `1.5` = res_gamma_pdf_gw_32_method2,
  `1.0` = res_gamma_pdf_gw_42_method2,
  `0.5` = res_gamma_pdf_gw_52_method2
)

df_gamma_gw_method3 <- data.frame(
  `2.5` = res_gamma_pdf_gw_12_method3,
  `2.0` = res_gamma_pdf_gw_22_method3,
  `1.5` = res_gamma_pdf_gw_32_method3,
  `1.0` = res_gamma_pdf_gw_42_method3,
  `0.5` = res_gamma_pdf_gw_52_method3
)

df_gamma_gw_method4 <- data.frame(
  `2.5` = res_gamma_pdf_gw_12_method4,
  `2.0` = res_gamma_pdf_gw_22_method4,
  `1.5` = res_gamma_pdf_gw_32_method4,
  `1.0` = res_gamma_pdf_gw_42_method4,
  `0.5` = res_gamma_pdf_gw_52_method4
)

# Create the data frame
df_logis_gw_method1 <- data.frame(
  `2.5` = res_logis_pdf_gw_12_method1,
  `2.0` = res_logis_pdf_gw_22_method1,
  `1.5` = res_logis_pdf_gw_32_method1,
  `1.0` = res_logis_pdf_gw_42_method1,
  `0.5` = res_logis_pdf_gw_52_method1
)

df_logis_gw_method2 <- data.frame(
  `2.5` = res_logis_pdf_gw_12_method2,
  `2.0` = res_logis_pdf_gw_22_method2,
  `1.5` = res_logis_pdf_gw_32_method2,
  `1.0` = res_logis_pdf_gw_42_method2,
  `0.5` = res_logis_pdf_gw_52_method2
)

df_logis_gw_method3 <- data.frame(
  `2.5` = res_logis_pdf_gw_12_method3,
  `2.0` = res_logis_pdf_gw_22_method3,
  `1.5` = res_logis_pdf_gw_32_method3,
  `1.0` = res_logis_pdf_gw_42_method3,
  `0.5` = res_logis_pdf_gw_52_method3
)

df_logis_gw_method4 <- data.frame(
  `2.5` = res_logis_pdf_gw_12_method4,
  `2.0` = res_logis_pdf_gw_22_method4,
  `1.5` = res_logis_pdf_gw_32_method4,
  `1.0` = res_logis_pdf_gw_42_method4,
  `0.5` = res_logis_pdf_gw_52_method4
)

# Create the data frame
df_t_gw_method1 <- data.frame(
  `2.5` = res_t_pdf_gw_12_method1,
  `2.0` = res_t_pdf_gw_22_method1,
  `1.5` = res_t_pdf_gw_32_method1,
  `1.0` = res_t_pdf_gw_42_method1,
  `0.5` = res_t_pdf_gw_52_method1
)

df_t_gw_method2 <- data.frame(
  `2.5` = res_t_pdf_gw_12_method2,
  `2.0` = res_t_pdf_gw_22_method2,
  `1.5` = res_t_pdf_gw_32_method2,
  `1.0` = res_t_pdf_gw_42_method2,
  `0.5` = res_t_pdf_gw_52_method2
)

df_t_gw_method3 <- data.frame(
  `2.5` = res_t_pdf_gw_12_method3,
  `2.0` = res_t_pdf_gw_22_method3,
  `1.5` = res_t_pdf_gw_32_method3,
  `1.0` = res_t_pdf_gw_42_method3,
  `0.5` = res_t_pdf_gw_52_method3
)

df_t_gw_method4 <- data.frame(
  `2.5` = res_t_pdf_gw_12_method4,
  `2.0` = res_t_pdf_gw_22_method4,
  `1.5` = res_t_pdf_gw_32_method4,
  `1.0` = res_t_pdf_gw_42_method4,
  `0.5` = res_t_pdf_gw_52_method4
)

df_norm_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_norm_pdf_gw_12_method1,
  `2.0` = smoothed_res_norm_pdf_gw_22_method1,
  `1.5` = smoothed_res_norm_pdf_gw_32_method1,
  `1.0` = smoothed_res_norm_pdf_gw_42_method1,
  `0.5` = smoothed_res_norm_pdf_gw_52_method1
)

df_beta_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_beta_pdf_gw_12_method1,
  `2.0` = smoothed_res_beta_pdf_gw_22_method1,
  `1.5` = smoothed_res_beta_pdf_gw_32_method1,
  `1.0` = smoothed_res_beta_pdf_gw_42_method1,
  `0.5` = smoothed_res_beta_pdf_gw_52_method1
)

df_gamma_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_gamma_pdf_gw_12_method1,
  `2.0` = smoothed_res_gamma_pdf_gw_22_method1,
  `1.5` = smoothed_res_gamma_pdf_gw_32_method1,
  `1.0` = smoothed_res_gamma_pdf_gw_42_method1,
  `0.5` = smoothed_res_gamma_pdf_gw_52_method1
)

df_logis_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_logis_pdf_gw_12_method1,
  `2.0` = smoothed_res_logis_pdf_gw_22_method1,
  `1.5` = smoothed_res_logis_pdf_gw_32_method1,
  `1.0` = smoothed_res_logis_pdf_gw_42_method1,
  `0.5` = smoothed_res_logis_pdf_gw_52_method1
)

df_t_gw_method1_smoothed <- data.frame(
  `2.5` = smoothed_res_t_pdf_gw_12_method1,
  `2.0` = smoothed_res_t_pdf_gw_22_method1,
  `1.5` = smoothed_res_t_pdf_gw_32_method1,
  `1.0` = smoothed_res_t_pdf_gw_42_method1,
  `0.5` = smoothed_res_t_pdf_gw_52_method1
)

write.csv(df_norm_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_beta_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_gamma_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_logis_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1.csv", row.names = FALSE)
write.csv(df_t_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1.csv", row.names = FALSE)

write.csv(df_norm_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_beta_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_gamma_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_logis_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1_smoothed.csv", row.names = FALSE)
write.csv(df_t_gw_method1_smoothed, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1_smoothed.csv", row.names = FALSE)

write.csv(df_norm_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_beta_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_gamma_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_logis_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method2.csv", row.names = FALSE)
write.csv(df_t_gw_method2, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method2.csv", row.names = FALSE)

write.csv(df_norm_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_beta_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_gamma_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_logis_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method3.csv", row.names = FALSE)
write.csv(df_t_gw_method3, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method3.csv", row.names = FALSE)

write.csv(df_norm_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_beta_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_gamma_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_logis_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method4.csv", row.names = FALSE)
write.csv(df_t_gw_method4, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method4.csv", row.names = FALSE)


# Recreating the plots

colors1 <- rep("#92CAEB",5)
colors2 <- rep("#F3CF70",5)
colors3 <- rep("#66D7A5",5)
colors4 <- rep("#E6A4C6",5)
colors5 <- rep("#D55E00",5)

# METHOD 1 #
df_norm_n_method1  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1.csv")
df_beta_n_method1  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1.csv")
df_gamma_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1.csv")
df_logis_n_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1.csv")
df_t_n_method1     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1.csv")

df_norm_n_method1_smoothed  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1_smoothed.csv")
df_beta_n_method1_smoothed  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1_smoothed.csv")
df_gamma_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1_smoothed.csv")
df_logis_n_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1_smoothed.csv")
df_t_n_method1_smoothed     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1_smoothed.csv")

df_norm_gw_method1  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1.csv")
df_beta_gw_method1  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1.csv")
df_gamma_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1.csv")
df_logis_gw_method1 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1.csv")
df_t_gw_method1     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1.csv")

df_norm_gw_method1_smoothed  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1_smoothed.csv")
df_beta_gw_method1_smoothed  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1_smoothed.csv")
df_gamma_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1_smoothed.csv")
df_logis_gw_method1_smoothed <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1_smoothed.csv")
df_t_gw_method1_smoothed     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1_smoothed.csv")

# METHOD 2 #
df_norm_n_method2  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method2.csv")
df_beta_n_method2  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method2.csv")
df_gamma_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method2.csv")
df_logis_n_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method2.csv")
df_t_n_method2     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method2.csv")

df_norm_gw_method2  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method2.csv")
df_beta_gw_method2  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method2.csv")
df_gamma_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method2.csv")
df_logis_gw_method2 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method2.csv")
df_t_gw_method2     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method2.csv")

# METHOD 3 #
df_norm_n_method3  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method3.csv")
df_beta_n_method3  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method3.csv")
df_gamma_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method3.csv")
df_logis_n_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method3.csv")
df_t_n_method3     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method3.csv")

df_norm_gw_method3  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method3.csv")
df_beta_gw_method3  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method3.csv")
df_gamma_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method3.csv")
df_logis_gw_method3 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method3.csv")
df_t_gw_method3     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method3.csv")

# METHOD 4 #
df_norm_n_method4  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method4.csv")
df_beta_n_method4  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method4.csv")
df_gamma_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method4.csv")
df_logis_n_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method4.csv")
df_t_n_method4     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method4.csv")

df_norm_gw_method4  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method4.csv")
df_beta_gw_method4  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method4.csv")
df_gamma_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method4.csv")
df_logis_gw_method4 <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method4.csv")
df_t_gw_method4     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method4.csv")


############################
# METHOD 1
############################

# Normal - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=TRUE,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_n_method1$X1e2,
        df_norm_n_method1$X1e3,
        df_norm_n_method1$X1e4,
        df_norm_n_method1$X1e5,
        df_norm_n_method1$X1e6,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        xlab = '',
        col = colors1, ylim=c(0,0.2))
dev.off()

# Normal - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_gw_method1$X2.5, df_norm_gw_method1$X2.0, df_norm_gw_method1$X1.5,
        df_norm_gw_method1$X1.0, df_norm_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors1, ylim=c(0,0.4))
dev.off()

# Beta - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_n_method1$X1e2,
        df_beta_n_method1$X1e3,
        df_beta_n_method1$X1e4,
        df_beta_n_method1$X1e5,
        df_beta_n_method1$X1e6,
        xlab = '',
        names = c(expression(10^2),expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2, ylim=c(0,0.45))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Beta - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_gw_method1$X2.5, df_beta_gw_method1$X2.0, df_beta_gw_method1$X1.5,
        df_beta_gw_method1$X1.0, df_beta_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors2, ylim=c(0,1))
dev.off()

# Gamma - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_n_method1$X1e2,
        df_gamma_n_method1$X1e3,
        df_gamma_n_method1$X1e4,
        df_gamma_n_method1$X1e5,
        df_gamma_n_method1$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3, ylim=c(0,0.45))
dev.off()

# Gamma - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_gw_method1$X2.5, df_gamma_gw_method1$X2.0, df_gamma_gw_method1$X1.5,
        df_gamma_gw_method1$X1.0, df_gamma_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors3, ylim=c(0,1))
dev.off()

# Logistic - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_n_method1$X1e2,
        df_logis_n_method1$X1e3,
        df_logis_n_method1$X1e4,
        df_logis_n_method1$X1e5,
        df_logis_n_method1$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4, ylim=c(0,0.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Logistic - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_gw_method1$X2.5, df_logis_gw_method1$X2.0, df_logis_gw_method1$X1.5,
        df_logis_gw_method1$X1.0, df_logis_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors4, ylim=c(0,0.4))
dev.off()

# t - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_n_method1$X1e2,
        df_t_n_method1$X1e3,
        df_t_n_method1$X1e4,
        df_t_n_method1$X1e5,
        df_t_n_method1$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5, ylim=c(0,0.2))
dev.off()

# t - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method1.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_gw_method1$X2.5, df_t_gw_method1$X2.0, df_t_gw_method1$X1.5,
        df_t_gw_method1$X1.0, df_t_gw_method1$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors5, ylim=c(0,0.5))
dev.off()

plot_method1_smoothed_n <- function(df, file, cols, y_limits) {
  grDevices::pdf(file, width = 6, height = 6)
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
      font.axis=1,cex.main=2)
  boxplot(df$X1e2,
          df$X1e3,
          df$X1e4,
          df$X1e5,
          df$X1e6,
          xlab = '',
          names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
          col = cols, ylim = y_limits)
  dev.off()
}

plot_method1_smoothed_gw <- function(df, file, cols, y_limits) {
  grDevices::pdf(file, width = 6, height = 6)
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
      font.axis=1,cex.main=2)
  boxplot(df$X2.5,
          df$X2.0,
          df$X1.5,
          df$X1.0,
          df$X0.5,
          names = c("2.5","2.0","1.5","1.0","0.5"),
          xlab = '',
          col = cols, ylim = y_limits)
  dev.off()
}

############################
# METHOD 1 - SMOOTHED
############################

plot_method1_smoothed_n(
  df_norm_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method1-Smoothed.pdf",
  colors1,
  c(0,0.2)
)

plot_method1_smoothed_gw(
  df_norm_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method1-Smoothed.pdf",
  colors1,
  c(0,0.4)
)

plot_method1_smoothed_n(
  df_beta_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method1-Smoothed.pdf",
  colors2,
  c(0,0.45)
)

plot_method1_smoothed_gw(
  df_beta_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method1-Smoothed.pdf",
  colors2,
  c(0,1)
)

plot_method1_smoothed_n(
  df_gamma_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method1-Smoothed.pdf",
  colors3,
  c(0,0.45)
)

plot_method1_smoothed_gw(
  df_gamma_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method1-Smoothed.pdf",
  colors3,
  c(0,1)
)

plot_method1_smoothed_n(
  df_logis_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method1-Smoothed.pdf",
  colors4,
  c(0,0.2)
)

plot_method1_smoothed_gw(
  df_logis_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method1-Smoothed.pdf",
  colors4,
  c(0,0.4)
)

plot_method1_smoothed_n(
  df_t_n_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method1-Smoothed.pdf",
  colors5,
  c(0,0.2)
)

plot_method1_smoothed_gw(
  df_t_gw_method1_smoothed,
  "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method1-Smoothed.pdf",
  colors5,
  c(0,0.5)
)

############################
# METHOD 2
############################

# Normal - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_n_method2$X1e2,
        df_norm_n_method2$X1e3,
        df_norm_n_method2$X1e4,
        df_norm_n_method2$X1e5,
        df_norm_n_method2$X1e6,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        xlab = '',
        col = colors1, ylim=c(0,0.2))
dev.off()

# Normal - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_gw_method2$X2.5, df_norm_gw_method2$X2.0, df_norm_gw_method2$X1.5,
        df_norm_gw_method2$X1.0, df_norm_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors1, ylim=c(0,0.4))
dev.off()

# Beta - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_n_method2$X1e2,
        df_beta_n_method2$X1e3,
        df_beta_n_method2$X1e4,
        df_beta_n_method2$X1e5,
        df_beta_n_method2$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2, ylim=c(0,0.45))
dev.off()

# Beta - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_gw_method2$X2.5, df_beta_gw_method2$X2.0, df_beta_gw_method2$X1.5,
        df_beta_gw_method2$X1.0, df_beta_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors2, ylim=c(0,1))
dev.off()

# Gamma - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_n_method2$X1e2,
        df_gamma_n_method2$X1e3,
        df_gamma_n_method2$X1e4,
        df_gamma_n_method2$X1e5,
        df_gamma_n_method2$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3, ylim=c(0,0.45))
dev.off()

# Gamma - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_gw_method2$X2.5, df_gamma_gw_method2$X2.0, df_gamma_gw_method2$X1.5,
        df_gamma_gw_method2$X1.0, df_gamma_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors3, ylim=c(0,1))
dev.off()

# Logistic - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_n_method2$X1e2,
        df_logis_n_method2$X1e3,
        df_logis_n_method2$X1e4,
        df_logis_n_method2$X1e5,
        df_logis_n_method2$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4, ylim=c(0,0.15))
dev.off()

# Logistic - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_gw_method2$X2.5, df_logis_gw_method2$X2.0, df_logis_gw_method2$X1.5,
        df_logis_gw_method2$X1.0, df_logis_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors4, ylim=c(0,0.4))
dev.off()

# t - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_n_method2$X1e2,
        df_t_n_method2$X1e3,
        df_t_n_method2$X1e4,
        df_t_n_method2$X1e5,
        df_t_n_method2$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5, ylim=c(0,0.2))
dev.off()

# t - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method2.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_gw_method2$X2.5, df_t_gw_method2$X2.0, df_t_gw_method2$X1.5,
        df_t_gw_method2$X1.0, df_t_gw_method2$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors5, ylim=c(0,0.5))
dev.off()

############################
# METHOD 3
############################

# Normal - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_n_method3$X1e2,
        df_norm_n_method3$X1e3,
        df_norm_n_method3$X1e4,
        df_norm_n_method3$X1e5,
        df_norm_n_method3$X1e6,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        xlab = '',
        col = colors1, ylim=c(0,0.2))
dev.off()

# Normal - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_gw_method3$X2.5, df_norm_gw_method3$X2.0, df_norm_gw_method3$X1.5,
        df_norm_gw_method3$X1.0, df_norm_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors1, ylim=c(0,0.4))
dev.off()

# Beta - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_n_method3$X1e2,
        df_beta_n_method3$X1e3,
        df_beta_n_method3$X1e4,
        df_beta_n_method3$X1e5,
        df_beta_n_method3$X1e6,
        xlab = '',
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2, ylim=c(0,0.45))
dev.off()

# Beta - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_gw_method3$X2.5, df_beta_gw_method3$X2.0, df_beta_gw_method3$X1.5,
        df_beta_gw_method3$X1.0, df_beta_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors2, ylim=c(0,1))
dev.off()

# Gamma - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_n_method3$X1e2,
        df_gamma_n_method3$X1e3,
        df_gamma_n_method3$X1e4,
        df_gamma_n_method3$X1e5,
        df_gamma_n_method3$X1e6,
        xlab = '',
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3, ylim=c(0,0.45))
dev.off()

# Gamma - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_gw_method3$X2.5, df_gamma_gw_method3$X2.0, df_gamma_gw_method3$X1.5,
        df_gamma_gw_method3$X1.0, df_gamma_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors3, ylim=c(0,1))
dev.off()

# Logistic - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_n_method3$X1e2,
        df_logis_n_method3$X1e3,
        df_logis_n_method3$X1e4,
        df_logis_n_method3$X1e5,
        df_logis_n_method3$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4, ylim=c(0,0.15))
dev.off()

# Logistic - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_gw_method3$X2.5, df_logis_gw_method3$X2.0, df_logis_gw_method3$X1.5,
        df_logis_gw_method3$X1.0, df_logis_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors4, ylim=c(0,0.4))
dev.off()

# t - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_n_method3$X1e2,
        df_t_n_method3$X1e3,
        df_t_n_method3$X1e4,
        df_t_n_method3$X1e5,
        df_t_n_method3$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5, ylim=c(0,0.2))
dev.off()

# t - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method3.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_gw_method3$X2.5, df_t_gw_method3$X2.0, df_t_gw_method3$X1.5,
        df_t_gw_method3$X1.0, df_t_gw_method3$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors5, ylim=c(0,0.5))
dev.off()

############################
# METHOD 4
############################

# Normal - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_n_method4$X1e2,
        df_norm_n_method4$X1e3,
        df_norm_n_method4$X1e4,
        df_norm_n_method4$X1e5,
        df_norm_n_method4$X1e6,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        xlab = '',
        col = colors1, ylim=c(0,0.2))
dev.off()

# Normal - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_gw_method4$X2.5, df_norm_gw_method4$X2.0, df_norm_gw_method4$X1.5,
        df_norm_gw_method4$X1.0, df_norm_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors1, ylim=c(0,0.4))
dev.off()

# Beta - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_n_method4$X1e2,
        df_beta_n_method4$X1e3,
        df_beta_n_method4$X1e4,
        df_beta_n_method4$X1e5,
        df_beta_n_method4$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2, ylim=c(0,0.45))
dev.off()

# Beta - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_gw_method4$X2.5, df_beta_gw_method4$X2.0, df_beta_gw_method4$X1.5,
        df_beta_gw_method4$X1.0, df_beta_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors2, ylim=c(0,1))
dev.off()

# Gamma - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_n_method4$X1e2,
        df_gamma_n_method4$X1e3,
        df_gamma_n_method4$X1e4,
        df_gamma_n_method4$X1e5,
        df_gamma_n_method4$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3, ylim=c(0,0.45))
dev.off()

# Gamma - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_gw_method4$X2.5, df_gamma_gw_method4$X2.0, df_gamma_gw_method4$X1.5,
        df_gamma_gw_method4$X1.0, df_gamma_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors3, ylim=c(0,1))
dev.off()

# Logistic - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_n_method4$X1e2,
        df_logis_n_method4$X1e3,
        df_logis_n_method4$X1e4,
        df_logis_n_method4$X1e5,
        df_logis_n_method4$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4, ylim=c(0,0.15))
dev.off()

# Logistic - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_gw_method4$X2.5, df_logis_gw_method4$X2.0, df_logis_gw_method4$X1.5,
        df_logis_gw_method4$X1.0, df_logis_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors4, ylim=c(0,0.4))
dev.off()

# t - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_n_method4$X1e2,
        df_t_n_method4$X1e3,
        df_t_n_method4$X1e4,
        df_t_n_method4$X1e5,
        df_t_n_method4$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5, ylim=c(0,0.2))
dev.off()

# t - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method4.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_gw_method4$X2.5, df_t_gw_method4$X2.0, df_t_gw_method4$X1.5,
        df_t_gw_method4$X1.0, df_t_gw_method4$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors5, ylim=c(0,0.5))
dev.off()

