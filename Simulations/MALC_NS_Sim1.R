# Parameters

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
  
  # res_norm_pdf_n_11_method1[i] = L1_Distance_calc(x,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_12_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_12_method2[i] = L2_from_res(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_12_method3[i] = L2_binnednp(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_12_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm,range_norm)
  
  # Normal Distribution - n2
  
  x <- rnorm(n[2], mean = mu,sd=sigma)
  
  # res_norm_pdf_n_21_method1[i] = L1_Distance_calc(x,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_22_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_22_method2[i] = L2_from_res(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_22_method3[i] = L2_binnednp(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_22_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm,range_norm)
  
  # Normal Distribution - n3
  
  x <- rnorm(n[3], mean = mu,sd=sigma)
  
  # res_norm_pdf_n_31_method1[i] = L1_Distance_calc(x,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_32_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_32_method2[i] = L2_from_res(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_32_method3[i] = L2_binnednp(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_32_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm,range_norm)
  
  # Normal Distribution - n4
  
  x <- rnorm(n[4], mean = mu,sd=sigma)
  
  # res_norm_pdf_n_41_method1[i] = L1_Distance_calc(x,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_42_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_42_method2[i] = L2_from_res(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_42_method3[i] = L2_binnednp(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_42_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm,range_norm)
  
  # Normal Distribution - n5
  
  x <- rnorm(n[5], mean = mu,sd=sigma)
  
  # res_norm_pdf_n_51_method1[i] = L1_Distance_calc(x,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_52_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_52_method2[i] = L2_from_res(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_52_method3[i] = L2_binnednp(x,norm_grid,pdf_norm,range_norm)
  # res_norm_pdf_n_52_method4[i] = L2_from_kernsmooth(x,norm_grid,pdf_norm,range_norm)
  
  print(i)
  
}

# Method-1 # 

par(mfrow = c(2,2))

par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(3,2),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
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
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
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
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
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

dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
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
dev.off()


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
  # res_beta_pdf_n_11_method1[i] = L1_Distance_calc(x,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_12_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_12_method2[i] = L2_from_res(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_12_method3[i] = L2_binnednp(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_12_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta,range_beta)
  
  # Beta Distribution - n2
  
  x   <- rbeta(n[2],alpha,beta)
  # res_beta_pdf_n_21_method1[i] = L1_Distance_calc(x,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_22_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_22_method2[i] = L2_from_res(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_22_method3[i] = L2_binnednp(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_22_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta,range_beta)
  
  # Beta Distribution - n3
  
  x   <- rbeta(n[3],alpha,beta)
  # res_beta_pdf_n_31_method1[i] = L1_Distance_calc(x,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_32_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_32_method2[i] = L2_from_res(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_32_method3[i] = L2_binnednp(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_32_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta,range_beta)
  
  # Beta Distribution - n4
  
  x   <- rbeta(n[4],alpha,beta)
  # res_beta_pdf_n_41_method1[i] = L1_Distance_calc(x,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_42_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_42_method2[i] = L2_from_res(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_42_method3[i] = L2_binnednp(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_42_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta,range_beta)
  
  # Beta Distribution - n5
  
  x   <- rbeta(n[5],alpha,beta)
  # res_beta_pdf_n_51_method1[i] = L1_Distance_calc(x,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_52_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_52_method2[i] = L2_from_res(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_52_method3[i] = L2_binnednp(x,beta_grid,pdf_beta,range_beta)
  # res_beta_pdf_n_52_method4[i] = L2_from_kernsmooth(x,beta_grid,pdf_beta,range_beta)
  
  print(i)
  
}

# Method # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-n/L1-PDF-1-Beta-n-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11_method1,res_beta_pdf_n_21_method1,res_beta_pdf_n_31_method1,res_beta_pdf_n_41_method1,res_beta_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

par(mfrow = c(2,2))
grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method1,res_beta_pdf_n_22_method1,res_beta_pdf_n_32_method1,res_beta_pdf_n_42_method1,res_beta_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method2,res_beta_pdf_n_22_method2,res_beta_pdf_n_32_method2,res_beta_pdf_n_42_method2,res_beta_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method3,res_beta_pdf_n_22_method3,res_beta_pdf_n_32_method3,res_beta_pdf_n_42_method3,res_beta_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method4,res_beta_pdf_n_22_method4,res_beta_pdf_n_32_method4,res_beta_pdf_n_42_method4,res_beta_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))
dev.off()


# Gamma Distribution - Sample Size Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grid <- seq(0,MM, by = gamma_deltas[1])

# Method-1 #

res_gamma_pdf_n_11_method1 = rep(0,B)
res_gamma_pdf_n_21_method1 = rep(0,B)
res_gamma_pdf_n_31_method1 = rep(0,B)
res_gamma_pdf_n_41_method1 = rep(0,B)
res_gamma_pdf_n_51_method1 = rep(0,B)

res_gamma_pdf_n_12_method1 = rep(0,B)
res_gamma_pdf_n_22_method1 = rep(0,B)
res_gamma_pdf_n_32_method1 = rep(0,B)
res_gamma_pdf_n_42_method1 = rep(0,B)
res_gamma_pdf_n_52_method1 = rep(0,B)

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
  
  # res_gamma_pdf_n_11_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_12_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_12_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_12_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_12_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n2
  
  x   <- (rgamma(n[2],alpha,beta))
  
  # res_gamma_pdf_n_21_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_22_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_22_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_22_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_22_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n3
  
  x   <- (rgamma(n[3],alpha,beta))
  
  # res_gamma_pdf_n_31_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_32_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_32_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_32_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_32_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n4
  
  x   <- (rgamma(n[4],alpha,beta))
  
  # res_gamma_pdf_n_41_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_42_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_42_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_42_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_42_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n5
  
  x   <- (rgamma(n[5],alpha,beta))
  
  # res_gamma_pdf_n_51_method1[i] = L1_Distance_calc(x,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_52_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_52_method2[i] = L2_from_res(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_52_method3[i] = L2_binnednp(x,gamma_grid,pdf_gamma,range_gamma)
  # res_gamma_pdf_n_52_method4[i] = L2_from_kernsmooth(x,gamma_grid,pdf_gamma,range_gamma)
  
  print(i)
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-n/L1-PDF-1-Gamma-n-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11_method1,res_gamma_pdf_n_21_method1,res_gamma_pdf_n_31_method1,res_gamma_pdf_n_41_method1,res_gamma_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method1,res_gamma_pdf_n_22_method1,res_gamma_pdf_n_32_method1,res_gamma_pdf_n_42_method1,res_gamma_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.3))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method2,res_gamma_pdf_n_22_method2,res_gamma_pdf_n_32_method2,res_gamma_pdf_n_42_method2,res_gamma_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.3))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method3,res_gamma_pdf_n_22_method3,res_gamma_pdf_n_32_method3,res_gamma_pdf_n_42_method3,res_gamma_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.3))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method4,res_gamma_pdf_n_22_method4,res_gamma_pdf_n_32_method4,res_gamma_pdf_n_42_method4,res_gamma_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.3))
dev.off()

# Logistic Distribution - Sample Size Simulation 

MM = qlogis(0.9999999,mu,sigma) + 10
logistic_grid <- seq(-MM,MM, by = logistic_deltas[1])

# Method-1 # 

res_logis_pdf_n_11_method1 = rep(0,B)
res_logis_pdf_n_21_method1 = rep(0,B)
res_logis_pdf_n_31_method1 = rep(0,B)
res_logis_pdf_n_41_method1 = rep(0,B)
res_logis_pdf_n_51_method1 = rep(0,B)

res_logis_pdf_n_12_method1 = rep(0,B)
res_logis_pdf_n_22_method1 = rep(0,B)
res_logis_pdf_n_32_method1 = rep(0,B)
res_logis_pdf_n_42_method1 = rep(0,B)
res_logis_pdf_n_52_method1 = rep(0,B)

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
  
  # res_logis_pdf_n_11_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_12_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_12_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_12_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_12_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n2
  
  x = (rlogis(n[2]))
  
  # res_logis_pdf_n_21_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_22_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_22_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_22_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_22_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n3
  
  x = sort(rlogis(n[3]))
  
  # res_logis_pdf_n_31_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_32_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_32_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_32_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_32_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n4
  
  x = (rlogis(n[4]))
  
  # res_logis_pdf_n_41_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_42_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_42_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_42_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_42_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n5
  
  x = (rlogis(n[5]))
  
  # res_logis_pdf_n_51_method1[i] = L1_Distance_calc(x,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_52_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_52_method2[i] = L2_from_res(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_52_method3[i] = L2_binnednp(x,logistic_grid,pdf_logistic,range_logistic)
  # res_logis_pdf_n_52_method4[i] = L2_from_kernsmooth(x,logistic_grid,pdf_logistic,range_logistic)
  
  print(i)
  
}

# Method #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-n/L1-PDF-1-Logistic-n-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11_method1,res_logis_pdf_n_21_method1,res_logis_pdf_n_31_method1,res_logis_pdf_n_41_method1,res_logis_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method1,res_logis_pdf_n_22_method1,res_logis_pdf_n_32_method1,res_logis_pdf_n_42_method1,res_logis_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method2,res_logis_pdf_n_22_method2,res_logis_pdf_n_32_method2,res_logis_pdf_n_42_method2,res_logis_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method3,res_logis_pdf_n_22_method3,res_logis_pdf_n_32_method3,res_logis_pdf_n_42_method3,res_logis_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method4,res_logis_pdf_n_22_method4,res_logis_pdf_n_32_method4,res_logis_pdf_n_42_method4,res_logis_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.2))
dev.off()

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
res_t_pdf_n_21_method1 = rep(0,B)
res_t_pdf_n_31_method1 = rep(0,B)
res_t_pdf_n_41_method1 = rep(0,B)
res_t_pdf_n_51_method1 = rep(0,B)

res_t_pdf_n_12_method1 = rep(0,B)
res_t_pdf_n_22_method1 = rep(0,B)
res_t_pdf_n_32_method1 = rep(0,B)
res_t_pdf_n_42_method1 = rep(0,B)
res_t_pdf_n_52_method1 = rep(0,B)

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
  
  # res_t_pdf_n_11_method1[i] = L1_Distance_calc(x,t_grid,pdf_t,range_t)
  res_t_pdf_n_12_method1[i] = L2_Distance_calc_nonlogconc(x,t_grid1,pdf_t,range_t)
  # res_t_pdf_n_12_method2[i] = L2_from_res(x,t_grid1,pdf_t,range_t)
  # res_t_pdf_n_12_method3[i] = L2_binnednp(x,t_grid1,pdf_t,range_t)
  # res_t_pdf_n_12_method4[i] = L2_from_kernsmooth(x,t_grid1,pdf_t,range_t)
  
  # t Distribution - n2
  
  x = (rt(n[2],df = df))
  t_grid2 <- make_breaks_cover(x, t_grid2)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_21_method1[i] = L1_Distance_calc(x,t_grid,pdf_t,range_t)
  res_t_pdf_n_22_method1[i] = L2_Distance_calc_nonlogconc(x,t_grid2,pdf_t,range_t)
  # res_t_pdf_n_22_method2[i] = L2_from_res(x,t_grid2,pdf_t,range_t)
  # res_t_pdf_n_22_method3[i] = L2_binnednp(x,t_grid2,pdf_t,range_t)
  # res_t_pdf_n_22_method4[i] = L2_from_kernsmooth(x,t_grid2,pdf_t,range_t)
  
  # t Distribution - n3
  
  x = (rt(n[3],df = df))
  t_grid3 <- make_breaks_cover(x, t_grid3)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_31_method1[i] = L1_Distance_calc(x,t_grid,pdf_t,range_t)
  res_t_pdf_n_32_method1[i] = L2_Distance_calc_nonlogconc(x,t_grid3,pdf_t,range_t)
  # res_t_pdf_n_32_method2[i] = L2_from_res(x,t_grid3,pdf_t,range_t)
  # res_t_pdf_n_32_method3[i] = L2_binnednp(x,t_grid3,pdf_t,range_t)
  # res_t_pdf_n_32_method4[i] = L2_from_kernsmooth(x,t_grid3,pdf_t,range_t)
  
  # t Distribution - n4
  
  x = (rt(n[4],df = df))
  t_grid4 <- make_breaks_cover(x, t_grid4)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_41_method1[i] = L1_Distance_calc(x,t_grid,pdf_t,range_t)
  res_t_pdf_n_42_method1[i] = L2_Distance_calc_nonlogconc(x,t_grid4,pdf_t,range_t)
  # res_t_pdf_n_42_method2[i] = L2_from_res(x,t_grid4,pdf_t,range_t)
  # res_t_pdf_n_42_method3[i] = L2_binnednp(x,t_grid4,pdf_t,range_t)
  # res_t_pdf_n_42_method4[i] = L2_from_kernsmooth(x,t_grid4,pdf_t,range_t)
  
  # t Distribution - n5
  
  x = (rt(n[5],df = df))
  t_grid5 <- make_breaks_cover(x, t_grid5)   # use the matching grid (t_grid1...t_grid5)
  
  # res_t_pdf_n_51_method1[i] = L1_Distance_calc(x,t_grid,pdf_t,range_t)
  res_t_pdf_n_52_method1[i] = L2_Distance_calc_nonlogconc(x,t_grid5,pdf_t,range_t)
  # res_t_pdf_n_52_method2[i] = L2_from_res(x,t_grid5,pdf_t,range_t)
  # res_t_pdf_n_52_method3[i] = L2_binnednp(x,t_grid5,pdf_t,range_t)
  # res_t_pdf_n_52_method4[i] = L2_from_kernsmooth(x,t_grid5,pdf_t,range_t)
  
  print(i)
  
}


# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-n/L1-PDF-1-t-n-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11_method1,res_t_pdf_n_21_method1,res_t_pdf_n_31_method1,res_t_pdf_n_41_method1,res_t_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method1,res_t_pdf_n_22_method1,res_t_pdf_n_32_method1,res_t_pdf_n_42_method1,res_t_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method2,res_t_pdf_n_22_method2,res_t_pdf_n_32_method2,res_t_pdf_n_42_method2,res_t_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method3,res_t_pdf_n_22_method3,res_t_pdf_n_32_method3,res_t_pdf_n_42_method3,res_t_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method4,res_t_pdf_n_22_method4,res_t_pdf_n_32_method4,res_t_pdf_n_42_method4,res_t_pdf_n_52_method4,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
dev.off()

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
res_norm_pdf_gw_21_method1 = rep(0,B)
res_norm_pdf_gw_31_method1 = rep(0,B)
res_norm_pdf_gw_41_method1 = rep(0,B)
res_norm_pdf_gw_51_method1 = rep(0,B)

res_norm_pdf_gw_12_method1 = rep(0,B)
res_norm_pdf_gw_22_method1 = rep(0,B)
res_norm_pdf_gw_32_method1 = rep(0,B)
res_norm_pdf_gw_42_method1 = rep(0,B)
res_norm_pdf_gw_52_method1 = rep(0,B)

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
  
  # res_norm_pdf_gw_11_method1[i] = L1_Distance_calc(x,norm_grids[[1]],pdf_norm,range_norm)
  res_norm_pdf_gw_12_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grids[[1]],pdf_norm,range_norm)
  # res_norm_pdf_gw_12_method2[i] = L2_from_res(x,norm_grids[[1]],pdf_norm,range_norm)
  # res_norm_pdf_gw_12_method3[i] = L2_binnednp(x,norm_grids[[1]],pdf_norm,range_norm)
  # res_norm_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,norm_grids[[1]],pdf_norm,range_norm)
  
  # Normal Distribution - GW2
  
  # res_norm_pdf_gw_21_method1[i] = L1_Distance_calc(x,norm_grids[[2]],pdf_norm,range_norm)
  res_norm_pdf_gw_22_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grids[[2]],pdf_norm,range_norm)
  # res_norm_pdf_gw_22_method2[i] = L2_from_res(x,norm_grids[[2]],pdf_norm,range_norm)
  # res_norm_pdf_gw_22_method3[i] = L2_binnednp(x,norm_grids[[2]],pdf_norm,range_norm)
  # res_norm_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,norm_grids[[2]],pdf_norm,range_norm)
  
  # Normal Distribution - GW3
  
  # res_norm_pdf_gw_31_method1[i] = L1_Distance_calc(x,norm_grids[[3]],pdf_norm,range_norm)
  res_norm_pdf_gw_32_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grids[[3]],pdf_norm,range_norm)
  # res_norm_pdf_gw_32_method2[i] = L2_from_res(x,norm_grids[[3]],pdf_norm,range_norm)
  # res_norm_pdf_gw_32_method3[i] = L2_binnednp(x,norm_grids[[3]],pdf_norm,range_norm)
  # res_norm_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,norm_grids[[3]],pdf_norm,range_norm)
  
  # Normal Distribution - GW4
  
  # res_norm_pdf_gw_41_method1[i] = L1_Distance_calc(x,norm_grids[[4]],pdf_norm,range_norm)
  res_norm_pdf_gw_42_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grids[[4]],pdf_norm,range_norm)
  # res_norm_pdf_gw_42_method2[i] = L2_from_res(x,norm_grids[[4]],pdf_norm,range_norm)
  # res_norm_pdf_gw_42_method3[i] = L2_binnednp(x,norm_grids[[4]],pdf_norm,range_norm)
  # res_norm_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,norm_grids[[4]],pdf_norm,range_norm)
  
  # Normal Distribution - GW5
  
  # res_norm_pdf_gw_51_method1[i] = L1_Distance_calc(x,norm_grids[[5]],pdf_norm,range_norm)
  res_norm_pdf_gw_52_method1[i] = L2_Distance_calc_nonlogconc(x,norm_grids[[5]],pdf_norm,range_norm)
  # res_norm_pdf_gw_52_method2[i] = L2_from_res(x,norm_grids[[5]],pdf_norm,range_norm)
  # res_norm_pdf_gw_52_method3[i] = L2_binnednp(x,norm_grids[[5]],pdf_norm,range_norm)
  # res_norm_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,norm_grids[[5]],pdf_norm,range_norm)
  
  print(i)
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-GW/L1-PDF-1-Normal-gw-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11_method1,res_norm_pdf_gw_21_method1,res_norm_pdf_gw_31_method1,res_norm_pdf_gw_41_method1,res_norm_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method1,res_norm_pdf_gw_22_method1,res_norm_pdf_gw_32_method1,res_norm_pdf_gw_42_method1,res_norm_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method2,res_norm_pdf_gw_22_method2,res_norm_pdf_gw_32_method2,res_norm_pdf_gw_42_method2,res_norm_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method3,res_norm_pdf_gw_22_method3,res_norm_pdf_gw_32_method3,res_norm_pdf_gw_42_method3,res_norm_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method4,res_norm_pdf_gw_22_method4,res_norm_pdf_gw_32_method4,res_norm_pdf_gw_42_method4,res_norm_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()


# Beta Distribution - Grid Width Simulation 

beta_grids = list(seq(0,1.2, by = beta_deltas[5]),
                  seq(0,1.3, by = beta_deltas[4]),
                  seq(0,1.2, by = beta_deltas[3]),
                  seq(0,1.2, by = beta_deltas[2]),
                  seq(0,1.1, by = beta_deltas[1]))

# Method-1 #

res_beta_pdf_gw_11_method1 = rep(0,B)
res_beta_pdf_gw_21_method1 = rep(0,B)
res_beta_pdf_gw_31_method1 = rep(0,B)
res_beta_pdf_gw_41_method1 = rep(0,B)
res_beta_pdf_gw_51_method1 = rep(0,B)

res_beta_pdf_gw_12_method1 = rep(0,B)
res_beta_pdf_gw_22_method1 = rep(0,B)
res_beta_pdf_gw_32_method1 = rep(0,B)
res_beta_pdf_gw_42_method1 = rep(0,B)
res_beta_pdf_gw_52_method1 = rep(0,B)

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
  
  # res_beta_pdf_gw_11_method1[i] = L1_Distance_calc(x,beta_grids[[1]],pdf_beta,range_beta)
  res_beta_pdf_gw_12_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grids[[1]],pdf_beta,range_beta)
  # res_beta_pdf_gw_12_method2[i] = L2_from_res(x,beta_grids[[1]],pdf_beta,range_beta)
  # res_beta_pdf_gw_12_method3[i] = L2_binnednp(x,beta_grids[[1]],pdf_beta,range_beta)
  # res_beta_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,beta_grids[[1]],pdf_beta,range_beta)
  
  
  # Beta Distribution - GW2
  
  # res_beta_pdf_gw_21_method1[i] = L1_Distance_calc(x,beta_grids[[2]],pdf_beta,range_beta)
  res_beta_pdf_gw_22_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grids[[2]],pdf_beta,range_beta)
  # res_beta_pdf_gw_22_method2[i] = L2_from_res(x,beta_grids[[2]],pdf_beta,range_beta)
  # res_beta_pdf_gw_22_method3[i] = L2_binnednp(x,beta_grids[[2]],pdf_beta,range_beta)
  # res_beta_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,beta_grids[[2]],pdf_beta,range_beta)
  
  # Beta Distribution - GW3
  
  # res_beta_pdf_gw_31_method1[i] = L1_Distance_calc(x,beta_grids[[3]],pdf_beta,range_beta)
  res_beta_pdf_gw_32_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grids[[3]],pdf_beta,range_beta)
  # res_beta_pdf_gw_32_method2[i] = L2_from_res(x,beta_grids[[3]],pdf_beta,range_beta)
  # res_beta_pdf_gw_32_method3[i] = L2_binnednp(x,beta_grids[[3]],pdf_beta,range_beta)
  # res_beta_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,beta_grids[[3]],pdf_beta,range_beta)
  
  # Beta Distribution - GW4
  
  # res_beta_pdf_gw_41_method1[i] = L1_Distance_calc(x,beta_grids[[4]],pdf_beta,range_beta)
  res_beta_pdf_gw_42_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grids[[4]],pdf_beta,range_beta)
  # res_beta_pdf_gw_42_method2[i] = L2_from_res(x,beta_grids[[4]],pdf_beta,range_beta)
  # res_beta_pdf_gw_42_method3[i] = L2_binnednp(x,beta_grids[[4]],pdf_beta,range_beta)
  # res_beta_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,beta_grids[[4]],pdf_beta,range_beta)
  
  
  # Beta Distribution - GW5
  
  # res_beta_pdf_gw_51_method1[i] = L1_Distance_calc(x,beta_grids[[5]],pdf_beta,range_beta)
  res_beta_pdf_gw_52_method1[i] = L2_Distance_calc_nonlogconc(x,beta_grids[[5]],pdf_beta,range_beta)
  # res_beta_pdf_gw_52_method2[i] = L2_from_res(x,beta_grids[[5]],pdf_beta,range_beta)
  # res_beta_pdf_gw_52_method3[i] = L2_binnednp(x,beta_grids[[5]],pdf_beta,range_beta)
  # res_beta_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,beta_grids[[5]],pdf_beta,range_beta)
  
  print(i)
  
}

# Method-1 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-GW/L1-PDF-1-Beta-gw-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11_method1,res_beta_pdf_gw_21_method1,res_beta_pdf_gw_31_method1,res_beta_pdf_gw_41_method1,res_beta_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method1,res_beta_pdf_gw_22_method1,res_beta_pdf_gw_32_method1,res_beta_pdf_gw_42_method1,res_beta_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method2,res_beta_pdf_gw_22_method2,res_beta_pdf_gw_32_method2,res_beta_pdf_gw_42_method2,res_beta_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method3,res_beta_pdf_gw_22_method3,res_beta_pdf_gw_32_method3,res_beta_pdf_gw_42_method3,res_beta_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method4,res_beta_pdf_gw_22_method4,res_beta_pdf_gw_32_method4,res_beta_pdf_gw_42_method4,res_beta_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()


# Gamma Distribution - Grid Width Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grids = list(seq(0,MM, by = gamma_deltas[5]),
                   seq(0,MM, by = gamma_deltas[4]),
                   seq(0,MM, by = gamma_deltas[3]),
                   seq(0,MM, by = gamma_deltas[2]),
                   seq(0,MM, by = gamma_deltas[1]))

# Method-1 # 

res_gamma_pdf_gw_11_method1 = rep(0,B)
res_gamma_pdf_gw_21_method1 = rep(0,B)
res_gamma_pdf_gw_31_method1 = rep(0,B)
res_gamma_pdf_gw_41_method1 = rep(0,B)
res_gamma_pdf_gw_51_method1 = rep(0,B)

res_gamma_pdf_gw_12_method1 = rep(0,B)
res_gamma_pdf_gw_22_method1 = rep(0,B)
res_gamma_pdf_gw_32_method1 = rep(0,B)
res_gamma_pdf_gw_42_method1 = rep(0,B)
res_gamma_pdf_gw_52_method1 = rep(0,B)

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
  
  # res_gamma_pdf_gw_11_method1[i] = L1_Distance_calc(x,gamma_grids[[1]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_12_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grids[[1]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_12_method2[i] = L2_from_res(x,gamma_grids[[1]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_12_method3[i] = L2_binnednp(x,gamma_grids[[1]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,gamma_grids[[1]],pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW2
  
  # res_gamma_pdf_gw_21_method1[i] = L1_Distance_calc(x,gamma_grids[[2]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_22_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grids[[2]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_22_method2[i] = L2_from_res(x,gamma_grids[[2]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_22_method3[i] = L2_binnednp(x,gamma_grids[[2]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,gamma_grids[[2]],pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW3
  
  # res_gamma_pdf_gw_31_method1[i] = L1_Distance_calc(x,gamma_grids[[3]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_32_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grids[[3]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_32_method2[i] = L2_from_res(x,gamma_grids[[3]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_32_method3[i] = L2_binnednp(x,gamma_grids[[3]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,gamma_grids[[3]],pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW4
  
  # res_gamma_pdf_gw_41_method1[i] = L1_Distance_calc(x,gamma_grids[[4]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_42_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grids[[4]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_42_method2[i] = L2_from_res(x,gamma_grids[[4]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_42_method3[i] = L2_binnednp(x,gamma_grids[[4]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,gamma_grids[[4]],pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW5
  
  # res_gamma_pdf_gw_51_method1[i] = L1_Distance_calc(x,gamma_grids[[5]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_52_method1[i] = L2_Distance_calc_nonlogconc(x,gamma_grids[[5]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_52_method2[i] = L2_from_res(x,gamma_grids[[5]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_52_method3[i] = L2_binnednp(x,gamma_grids[[5]],pdf_gamma,range_gamma)
  # res_gamma_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,gamma_grids[[5]],pdf_gamma,range_gamma)
  
  print(i)
  
}

# Method # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-GW/L1-PDF-1-Gamma-gw-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11_method1,res_gamma_pdf_gw_21_method1,res_gamma_pdf_gw_31_method1,res_gamma_pdf_gw_41_method1,res_gamma_pdf_gw_51_method1,
        names = c("1 (7)", "0.8 (8)", "0.6 (11)", "0.4 (16)", "0.2 (31)"),
        col = colors3,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method1,res_gamma_pdf_gw_22_method1,res_gamma_pdf_gw_32_method1,res_gamma_pdf_gw_42_method1,res_gamma_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method2,res_gamma_pdf_gw_22_method2,res_gamma_pdf_gw_32_method2,res_gamma_pdf_gw_42_method2,res_gamma_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method3,res_gamma_pdf_gw_22_method3,res_gamma_pdf_gw_32_method3,res_gamma_pdf_gw_42_method3,res_gamma_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method4,res_gamma_pdf_gw_22_method4,res_gamma_pdf_gw_32_method4,res_gamma_pdf_gw_42_method4,res_gamma_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.4))
dev.off()


# Logistic Distribution - Grid Width Simulation 

MM = qlogis(0.9999999,mu,sigma) + 3
logistic_grids = list(seq(-MM,MM, by = logistic_deltas[5]),
                      seq(-MM,MM, by = logistic_deltas[4]),
                      seq(-MM,MM, by = logistic_deltas[3]),
                      seq(-MM,MM, by = logistic_deltas[2]),
                      seq(-MM,MM, by = logistic_deltas[1]))


# Method-1 # 

res_logis_pdf_gw_11_method1 = rep(0,B)
res_logis_pdf_gw_21_method1 = rep(0,B)
res_logis_pdf_gw_31_method1 = rep(0,B)
res_logis_pdf_gw_41_method1 = rep(0,B)
res_logis_pdf_gw_51_method1 = rep(0,B)

res_logis_pdf_gw_12_method1 = rep(0,B)
res_logis_pdf_gw_22_method1 = rep(0,B)
res_logis_pdf_gw_32_method1 = rep(0,B)
res_logis_pdf_gw_42_method1 = rep(0,B)
res_logis_pdf_gw_52_method1 = rep(0,B)

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

set.seed(10)

for (i in 1:B) {
  
  x = sort(rlogis(n))
  
  # Logistic Distribution - GW1
  
  # res_logis_pdf_gw_11_method1[i] = L1_Distance_calc(x,logistic_grids[[1]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_12_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grids[[1]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_12_method2[i] = L2_from_res(x,logistic_grids[[1]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_12_method3[i] = L2_binnednp(x,logistic_grids[[1]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,logistic_grids[[1]],pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW2
  
  # res_logis_pdf_gw_21_method1[i] = L1_Distance_calc(x,logistic_grids[[2]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_22_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grids[[2]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_22_method2[i] = L2_from_res(x,logistic_grids[[2]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_22_method3[i] = L2_binnednp(x,logistic_grids[[2]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,logistic_grids[[2]],pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW3
  
  # res_logis_pdf_gw_31_method1[i] = L1_Distance_calc(x,logistic_grids[[3]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_32_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grids[[3]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_32_method2[i] = L2_from_res(x,logistic_grids[[3]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_32_method3[i] = L2_binnednp(x,logistic_grids[[3]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,logistic_grids[[3]],pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW4
  
  # res_logis_pdf_gw_41_method1[i] = L1_Distance_calc(x,logistic_grids[[4]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_42_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grids[[4]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_42_method2[i] = L2_from_res(x,logistic_grids[[4]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_42_method3[i] = L2_binnednp(x,logistic_grids[[4]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,logistic_grids[[4]],pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW5
  
  # res_logis_pdf_gw_51_method1[i] = L1_Distance_calc(x,logistic_grids[[5]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_52_method1[i] = L2_Distance_calc_nonlogconc(x,logistic_grids[[5]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_52_method2[i] = L2_from_res(x,logistic_grids[[5]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_52_method3[i] = L2_binnednp(x,logistic_grids[[5]],pdf_logistic,range_logistic)
  # res_logis_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,logistic_grids[[5]],pdf_logistic,range_logistic)
  
  print(i)
  
}

# Method # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-GW/L1-PDF-1-Logistic-gw-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11_method1,res_logis_pdf_gw_21_method1,res_logis_pdf_gw_31_method1,res_logis_pdf_gw_41_method1,res_logis_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method1,res_logis_pdf_gw_22_method1,res_logis_pdf_gw_32_method1,res_logis_pdf_gw_42_method1,res_logis_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method2,res_logis_pdf_gw_22_method2,res_logis_pdf_gw_32_method2,res_logis_pdf_gw_42_method2,res_logis_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method3,res_logis_pdf_gw_22_method3,res_logis_pdf_gw_32_method3,res_logis_pdf_gw_42_method3,res_logis_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method4,res_logis_pdf_gw_22_method4,res_logis_pdf_gw_32_method4,res_logis_pdf_gw_42_method4,res_logis_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# Student's t Distribution - Grid Width Simulation 

MM = qt(0.9999999,df) + 3
t_grids = list(seq(-MM,MM, by = t_deltas[5]),
               seq(-MM,MM, by = t_deltas[4]),
               seq(-MM,MM, by = t_deltas[3]),
               seq(-MM,MM, by = t_deltas[2]),
               seq(-MM,MM, by = t_deltas[1]))

# Method-1 # 

res_t_pdf_gw_11_method1 = rep(0,B)
res_t_pdf_gw_21_method1 = rep(0,B)
res_t_pdf_gw_31_method1 = rep(0,B)
res_t_pdf_gw_41_method1 = rep(0,B)
res_t_pdf_gw_51_method1 = rep(0,B)

res_t_pdf_gw_12_method1 = rep(0,B)
res_t_pdf_gw_22_method1 = rep(0,B)
res_t_pdf_gw_32_method1 = rep(0,B)
res_t_pdf_gw_42_method1 = rep(0,B)
res_t_pdf_gw_52_method1 = rep(0,B)

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

for (i in 1:B) {
  
  x = sort(rt(n,df = df))
  # t_grids[[1]] <- make_breaks_cover(x, t_grids[[1]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[2]] <- make_breaks_cover(x, t_grids[[2]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[3]] <- make_breaks_cover(x, t_grids[[3]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[4]] <- make_breaks_cover(x, t_grids[[4]])   # use the matching grid (t_grid1...t_grid5)
  # t_grids[[5]] <- make_breaks_cover(x, t_grids[[5]])   # use the matching grid (t_grid1...t_grid5)
  
  # t Distribution - GW1
  
  # res_t_pdf_gw_11_method1[i] = L1_Distance_calc(x,t_grids[[1]],pdf_t,range_t)
  res_t_pdf_gw_12_method1[i] = L2_Distance_calc_nonlogconc(x,t_grids[[1]],pdf_t,range_t)
  # res_t_pdf_gw_12_method2[i] = L2_from_res(x,t_grids[[1]],pdf_t,range_t)
  # res_t_pdf_gw_12_method3[i] = L2_binnednp(x,t_grids[[1]],pdf_t,range_t)
  # res_t_pdf_gw_12_method4[i] = L2_from_kernsmooth(x,t_grids[[1]],pdf_t,range_t)
  
  # t Distribution - GW2
  
  # res_t_pdf_gw_21_method1[i] = L1_Distance_calc(x,t_grids[[2]],pdf_t,range_t)
  res_t_pdf_gw_22_method1[i] = L2_Distance_calc_nonlogconc(x,t_grids[[2]],pdf_t,range_t)
  # res_t_pdf_gw_22_method2[i] = L2_from_res(x,t_grids[[2]],pdf_t,range_t)
  # res_t_pdf_gw_22_method3[i] = L2_binnednp(x,t_grids[[2]],pdf_t,range_t)
  # res_t_pdf_gw_22_method4[i] = L2_from_kernsmooth(x,t_grids[[2]],pdf_t,range_t)
  
  # t Distribution - GW3
  
  # res_t_pdf_gw_31_method1[i] = L1_Distance_calc(x,t_grids[[3]],pdf_t,range_t)
  res_t_pdf_gw_32_method1[i] = L2_Distance_calc_nonlogconc(x,t_grids[[3]],pdf_t,range_t)
  # res_t_pdf_gw_32_method2[i] = L2_from_res(x,t_grids[[3]],pdf_t,range_t)
  # res_t_pdf_gw_32_method3[i] = L2_binnednp(x,t_grids[[3]],pdf_t,range_t)
  # res_t_pdf_gw_32_method4[i] = L2_from_kernsmooth(x,t_grids[[3]],pdf_t,range_t)
  
  # t Distribution - GW4
  
  # res_t_pdf_gw_41_method1[i] = L1_Distance_calc(x,t_grids[[4]],pdf_t,range_t)
  res_t_pdf_gw_42_method1[i] = L2_Distance_calc_nonlogconc(x,t_grids[[4]],pdf_t,range_t)
  # res_t_pdf_gw_42_method2[i] = L2_from_res(x,t_grids[[4]],pdf_t,range_t)
  # res_t_pdf_gw_42_method3[i] = L2_binnednp(x,t_grids[[4]],pdf_t,range_t)
  # res_t_pdf_gw_42_method4[i] = L2_from_kernsmooth(x,t_grids[[4]],pdf_t,range_t)
  
  # t Distribution - GW5
  
  # res_t_pdf_gw_51_method1[i] = L1_Distance_calc(x,t_grids[[5]],pdf_t,range_t)
  res_t_pdf_gw_52_method1[i] = L2_Distance_calc_nonlogconc(x,t_grids[[5]],pdf_t,range_t)
  # res_t_pdf_gw_52_method2[i] = L2_from_res(x,t_grids[[5]],pdf_t,range_t)
  # res_t_pdf_gw_52_method3[i] = L2_binnednp(x,t_grids[[5]],pdf_t,range_t)
  # res_t_pdf_gw_52_method4[i] = L2_from_kernsmooth(x,t_grids[[5]],pdf_t,range_t)
  
  print(i)
  
}

# Method # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L1-GW/L1-PDF-1-t-gw-Method.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11_method1,res_t_pdf_gw_21_method1,res_t_pdf_gw_31_method1,res_t_pdf_gw_41_method1,res_t_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5),xlab('Grid Width'))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method1,res_t_pdf_gw_22_method1,res_t_pdf_gw_32_method1,res_t_pdf_gw_42_method1,res_t_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method2,res_t_pdf_gw_22_method2,res_t_pdf_gw_32_method2,res_t_pdf_gw_42_method2,res_t_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method3,res_t_pdf_gw_22_method3,res_t_pdf_gw_32_method3,res_t_pdf_gw_42_method3,res_t_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method4.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method4,res_t_pdf_gw_22_method4,res_t_pdf_gw_32_method4,res_t_pdf_gw_42_method4,res_t_pdf_gw_52_method4,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.5))
dev.off()

# Data Save

# Create the data frame
df_norm_n_method1 <- data.frame(
  `1e2` = res_norm_pdf_n_12_method1,
  `1e3` = res_norm_pdf_n_22_method1,
  `1e4` = res_norm_pdf_n_32_method1,
  `1e5` = res_norm_pdf_n_42_method1,
  `1e6` = res_norm_pdf_n_52_method1
)

# Create the data frame
df_beta_n_method1 <- data.frame(
  `1e2` = res_beta_pdf_n_12_method1,
  `1e3` = res_beta_pdf_n_22_method1,
  `1e4` = res_beta_pdf_n_32_method1,
  `1e5` = res_beta_pdf_n_42_method1,
  `1e6` = res_beta_pdf_n_52_method1
)


# Create the data frame
df_gamma_n_method1 <- data.frame(
  `1e2` = res_gamma_pdf_n_12_method1,
  `1e3` = res_gamma_pdf_n_22_method1,
  `1e4` = res_gamma_pdf_n_32_method1,
  `1e5` = res_gamma_pdf_n_42_method1,
  `1e6` = res_gamma_pdf_n_52_method1
)

# Create the data frame
df_logis_n_method1 <- data.frame(
  `1e2` = res_logis_pdf_n_12_method1,
  `1e3` = res_logis_pdf_n_22_method1,
  `1e4` = res_logis_pdf_n_32_method1,
  `1e5` = res_logis_pdf_n_42_method1,
  `1e6` = res_logis_pdf_n_52_method1
)


# Create the data frame
df_t_n_method1 <- data.frame(
  `1e2` = res_t_pdf_n_12_method1,
  `1e3` = res_t_pdf_n_22_method1,
  `1e4` = res_t_pdf_n_32_method1,
  `1e5` = res_t_pdf_n_42_method1,
  `1e6` = res_t_pdf_n_52_method1
)


write.csv(df_norm_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_beta_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_gamma_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_logis_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_t_n_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1_NS.csv", row.names = FALSE)


# Create the data frame
df_norm_gw_method1 <- data.frame(
  `2.5` = res_norm_pdf_gw_12_method1,
  `2.0` = res_norm_pdf_gw_22_method1,
  `1.5` = res_norm_pdf_gw_32_method1,
  `1.0` = res_norm_pdf_gw_42_method1,
  `0.5` = res_norm_pdf_gw_52_method1
)


# Create the data frame
df_beta_gw_method1 <- data.frame(
  `2.5` = res_beta_pdf_gw_12_method1,
  `2.0` = res_beta_pdf_gw_22_method1,
  `1.5` = res_beta_pdf_gw_32_method1,
  `1.0` = res_beta_pdf_gw_42_method1,
  `0.5` = res_beta_pdf_gw_52_method1
)


# Create the data frame
df_gamma_gw_method1 <- data.frame(
  `2.5` = res_gamma_pdf_gw_12_method1,
  `2.0` = res_gamma_pdf_gw_22_method1,
  `1.5` = res_gamma_pdf_gw_32_method1,
  `1.0` = res_gamma_pdf_gw_42_method1,
  `0.5` = res_gamma_pdf_gw_52_method1
)


# Create the data frame
df_logis_gw_method1 <- data.frame(
  `2.5` = res_logis_pdf_gw_12_method1,
  `2.0` = res_logis_pdf_gw_22_method1,
  `1.5` = res_logis_pdf_gw_32_method1,
  `1.0` = res_logis_pdf_gw_42_method1,
  `0.5` = res_logis_pdf_gw_52_method1
)

# Create the data frame
df_t_gw_method1 <- data.frame(
  `2.5` = res_t_pdf_gw_12_method1,
  `2.0` = res_t_pdf_gw_22_method1,
  `1.5` = res_t_pdf_gw_32_method1,
  `1.0` = res_t_pdf_gw_42_method1,
  `0.5` = res_t_pdf_gw_52_method1
)

write.csv(df_norm_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_beta_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_gamma_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_logis_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1_NS.csv", row.names = FALSE)
write.csv(df_t_gw_method1, file = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1_NS.csv", row.names = FALSE)

# Recreating the plots

colors1 <- rep("#92CAEB",5)
colors2 <- rep("#F3CF70",5)
colors3 <- rep("#66D7A5",5)
colors4 <- rep("#E6A4C6",5)
colors5 <- rep("#D55E00",5)

# METHOD 1 - NS #
df_norm_n_method1_NS  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_sample_sim_data_method1_NS.csv")
df_beta_n_method1_NS  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_sample_sim_data_method1_NS.csv")
df_gamma_n_method1_NS <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_sample_sim_data_method1_NS.csv")
df_logis_n_method1_NS <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_sample_sim_data_method1_NS.csv")
df_t_n_method1_NS     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_sample_sim_data_method1_NS.csv")

df_norm_gw_method1_NS  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/normal_grid_sim_data_method1_NS.csv")
df_beta_gw_method1_NS  <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/beta_grid_sim_data_method1_NS.csv")
df_gamma_gw_method1_NS <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/gamma_grid_sim_data_method1_NS.csv")
df_logis_gw_method1_NS <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/logistic_grid_sim_data_method1_NS.csv")
df_t_gw_method1_NS     <- read.csv("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/t_grid_sim_data_method1_NS.csv")


############################
# METHOD 1 - NS
############################

# Normal - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Normal-n-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_n_method1_NS$X1e2,
        df_norm_n_method1_NS$X1e3,
        df_norm_n_method1_NS$X1e4,
        df_norm_n_method1_NS$X1e5,
        df_norm_n_method1_NS$X1e6,
        names = c(expression(10^2),
                  expression(10^3),
                  expression(10^4),
                  expression(10^5),
                  expression(10^6)),
        xlab = '',
        col = colors1, ylim=c(0,0.2))
dev.off()

# Normal - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Normal-gw-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_norm_gw_method1_NS$X2.5, df_norm_gw_method1_NS$X2.0, df_norm_gw_method1_NS$X1.5,
        df_norm_gw_method1_NS$X1.0, df_norm_gw_method1_NS$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors1, ylim=c(0,0.4))
dev.off()

# Beta - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Beta-n-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_n_method1_NS$X1e2,
        df_beta_n_method1_NS$X1e3,
        df_beta_n_method1_NS$X1e4,
        df_beta_n_method1_NS$X1e5,
        df_beta_n_method1_NS$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors2, ylim=c(0,0.45))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Beta - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Beta-gw-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_beta_gw_method1_NS$X2.5, df_beta_gw_method1_NS$X2.0, df_beta_gw_method1_NS$X1.5,
        df_beta_gw_method1_NS$X1.0, df_beta_gw_method1_NS$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors2, ylim=c(0,1))
dev.off()

# Gamma - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Gamma-n-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_n_method1_NS$X1e2,
        df_gamma_n_method1_NS$X1e3,
        df_gamma_n_method1_NS$X1e4,
        df_gamma_n_method1_NS$X1e5,
        df_gamma_n_method1_NS$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors3, ylim=c(0,0.45))
dev.off()

# Gamma - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Gamma-gw-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_gamma_gw_method1_NS$X2.5, df_gamma_gw_method1_NS$X2.0, df_gamma_gw_method1_NS$X1.5,
        df_gamma_gw_method1_NS$X1.0, df_gamma_gw_method1_NS$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors3, ylim=c(0,1))
dev.off()

# Logistic - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-Logistic-n-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_n_method1_NS$X1e2,
        df_logis_n_method1_NS$X1e3,
        df_logis_n_method1_NS$X1e4,
        df_logis_n_method1_NS$X1e5,
        df_logis_n_method1_NS$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors4, ylim=c(0,0.15))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Logistic - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-Logistic-gw-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_logis_gw_method1_NS$X2.5, df_logis_gw_method1_NS$X2.0, df_logis_gw_method1_NS$X1.5,
        df_logis_gw_method1_NS$X1.0, df_logis_gw_method1_NS$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors4, ylim=c(0,0.4))
dev.off()

# t - Sample Size
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-n/L2-PDF-1-t-n-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_n_method1_NS$X1e2,
        df_t_n_method1_NS$X1e3,
        df_t_n_method1_NS$X1e4,
        df_t_n_method1_NS$X1e5,
        df_t_n_method1_NS$X1e6,
        xlab = '',
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
        col = colors5, ylim=c(0,0.2))
dev.off()

# t - Grid Width
grDevices::pdf("~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method/L2-GW/L2-PDF-1-t-gw-Method1_NS.pdf", width = 6, height = 6)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(df_t_gw_method1_NS$X2.5, df_t_gw_method1_NS$X2.0, df_t_gw_method1_NS$X1.5,
        df_t_gw_method1_NS$X1.0, df_t_gw_method1_NS$X0.5,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        xlab = '',
        col = colors5, ylim=c(0,0.5))
dev.off()



