library(jmuOutlier)
library(EnvStats)

####### NORMAL DISTRIBUTION #######

# N = 100
set.seed(10)
grid <- seq(-4, 4, by = 0.5)
x <- sort(rnorm(100, 0, 1))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Norm_100.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(-4,4), ylim = c(0,0.5))
plot(function(x) dnorm(x), xlim=c(-4,4), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)

abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

# N = 200

set.seed(10)
grid <- seq(-4, 4, by = 0.5)
x <- sort(rnorm(500, 0, 1))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Norm_200.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(-4,4), ylim = c(0,0.5))
plot(function(x) dnorm(x), xlim=c(-4,4), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()


# N = 1000

set.seed(10)
grid <- seq(-4, 4, by = 0.5)
x <- sort(rnorm(1000, 0, 1))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Norm_1000.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))

# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(-4,4), ylim = c(0,0.5))
plot(function(x) dnorm(x), xlim=c(-4,4), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()


####### GAMMA DISTRIBUTION #######

# N = 100 

alpha=2;beta=5;
set.seed(10)
grid <- seq(0, 3, by = 0.1414)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
x <- sort(rgamma(100, alpha, beta))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Gamma_100.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,2.5), ylim = c(0,2))
plot(function(x) dgamma(x,alpha,beta), xlim=c(0,2.5), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)


# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()


# N = 200

alpha=2;beta=5;
set.seed(10)
grid <- seq(0, 3, by = 0.1414)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
x <- sort(rgamma(200, alpha, beta))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Gamma_200.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,2.5), ylim = c(0,2))
plot(function(x) dgamma(x,alpha,beta), xlim=c(0,2.5), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)


# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

# N = 1000

alpha=2;beta=5;
set.seed(10)
grid <- seq(0, 3, by = 0.1414)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
x <- sort(rgamma(1000, alpha, beta))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Gamma_1000.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,2.5), ylim = c(0,2))
plot(function(x) dgamma(x,alpha,beta), xlim=c(0,2.5), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

####### PARETO DISTRIBUTION #######

# N = 100 
shape=2;location=2;
shape_p = 4;
set.seed(10)
grid <- seq(location, 20, by = 0.4714)
eval_grid <- seq(min(grid),max(grid), length.out=1001)

x <- sort(rpareto(100,location = location,shape = shape_p))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Pareto_100.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
plot(0, 0, type = 'n', xlim=c(2,7), ylim = c(0,2))
# Base: true density as a curve (continuous), no grid dependence
plot(function(x) dpareto(x,location,shape_p), xlim=c(2,7), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()


# N = 200

set.seed(11)
grid <- seq(location, 20, by = 0.4714)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
x <- sort(rpareto(200,location = location,shape = shape_p))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Pareto_200.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(2,7), ylim = c(0,2))
plot(function(x) dpareto(x,location,shape_p), xlim=c(2,7), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

# N = 1000

set.seed(10)
grid <- seq(location, 20, by = 0.4714)
eval_grid <- seq(min(grid),max(grid), length.out=1001)
x <- sort(rpareto(1000,location = location,shape = shape_p))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_Pareto_1000.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(2,7), ylim = c(0,2))
plot(function(x) dpareto(x,location,shape_p), xlim=c(2,7), 
     col="black", ylim=c(0,2),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

####### CHI-SQUARE DISTRIBUTION #######

# N = 100 
df=3
set.seed(10)
grid <- seq(0, 15, by = 1.2247)
eval_grid <- seq(min(grid),max(grid), length.out=1001)

x <- sort(rchisq(100, df))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_ChiSqr_100.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,10), ylim = c(0,0.5))
plot(function(x) dchisq(x,df), xlim=c(0,10), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

# N = 200

set.seed(10)
grid <- seq(0, 20, by = 1.2247)
eval_grid <- seq(min(grid),max(grid), length.out=1001)

x <- sort(rchisq(200, df))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_ChiSqr_200.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,19), ylim = c(0,0.5))
plot(function(x) dchisq(x,df), xlim=c(0,10), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()

# N = 1000

set.seed(10)
grid <- seq(0, 20, by = 1.2247)
eval_grid <- seq(min(grid),max(grid), length.out=1001)

x <- sort(rchisq(1000, df))

# Estimation 

# MALC 
fhat2 = get_fhatn(x,grid,log_conc = FALSE)
fhat2plot <- evaluateLogConDens(eval_grid, fhat2$fhatn, which=4)
fhat_malc <- fhat2plot[,5]

# Res 
hobj <- hist(x, breaks = grid, plot = FALSE)
res <- nonlinear_kde_binned_BK2002(
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  h = NULL,        # <- forces SJ-from-binned as suggested in the paper
  bw_B = 50,
  bw_method = "ste"
)
x_res <- eval_grid
fhat_res <- approx(res$x, res$fhat, xout = eval_grid, rule = 2)$y

# Binned 
fhat_ndpe <- binned_np.func_for_plot(x, grid,eval_grid) 

# KernelSmooth 
hobj <- hist(x, breaks=grid, plot=FALSE)
counts <- hobj$counts
centers <- (grid[-1] + grid[-length(grid)]) / 2

y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)
h <- dpik_with_fallback_jitter(
  y = y,
  counts = hobj$counts,
  breaks = hobj$breaks,
  grid = grid,
  seed = 1
)
fit <- bkde(y, bandwidth = h,
            kernel = "normal", canonical = FALSE,
            gridsize = length(grid), range.x = range(grid),
            truncate = TRUE)

x_ks <- eval_grid
fhat_ks <- approx(fit$x, fit$y, xout = eval_grid, rule = 2)$y

# Plot 

grDevices::pdf("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/MALC_ChiSqr_1000.pdf",width = 4, height = 4)
#par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",
#    font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
#    font.axis=1,cex.main=3)
par(mar=rep(0.5,4))
# Base: true density as a curve (continuous), no grid dependence
plot(0, 0, type = 'n', xlim=c(0,19), ylim = c(0,0.5))
plot(function(x) dchisq(x,df), xlim=c(0,10), 
     col="black", ylim=c(0,0.5),lwd = 1.5,add = TRUE)
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)

# Overlay estimated densities on their own x grids
lines(eval_grid,fhat_malc,col = "blue",lwd = 1.5)
lines(x_res,fhat_res,col = "green",lwd = 1.5)
lines(eval_grid,fhat_ndpe,col = "red",lwd = 1.5)
lines(x_ks,fhat_ks,col = "purple",lwd = 1.5)

legend("topright",
       legend = c("True Density", "Dens-OLog", "BK2002", "BinnedNP", "KernSmooth"),
       col    = c("black", "blue", "green", "red", "purple"),
       lwd    = 1.5,
       lty = 1,
       bty    = "n")

dev.off()
