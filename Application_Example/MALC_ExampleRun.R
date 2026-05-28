#----------------------------
# examples

set.seed(20180621)
n <- 100
grid <- seq(0,60, length.out=61)
x <- rgamma(n, shape=5, scale=1)

fhat05 <- get_fhatn(x, grid, alpha=0.5)
fhat1 <- get_fhatn(x, grid, alpha=1)
fhat2 <- get_fhatn(x, grid, alpha=2)
fhat_comp <- logConDens(x)

# this next step is where I change to the smoothing version

grid2 <- seq(min(grid),max(grid), length.out=1001)
fhat05plot <- evaluateLogConDens(grid2, fhat05$fhatn, which=4)
fhat05plot <- fhat05plot[,5]
fhat1plot <- evaluateLogConDens(grid2, fhat1$fhatn, which=4)
fhat1plot <- fhat1plot[,5]
fhat2plot <- evaluateLogConDens(grid2, fhat2$fhatn, which=4)
fhat2plot <- fhat2plot[,5]
fhatcompplot <- evaluateLogConDens(grid2, fhat_comp, which=4)
fhatcompplot <- fhatcompplot[,5]

#pdf(file = "gamma_eg_n200.pdf",   # The directory you want to save the file in
#    width = 4, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

par(mar=rep(0.5,4))
plot(0, 0, type = 'n', xlim = c(0,15), ylim = c(0, 0.25))
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)
#rug(res$x)
lines(grid2, dgamma(grid2, shape=5, scale=1), lwd=1.5)
lines(grid2, fhat05plot, col = "blue", lwd=1.5)
lines(grid2, fhat1plot, col = "green", lwd=1.5)
lines(grid2, fhat2plot, col = "red", lwd=1.5)
lines(grid2, fhatcompplot, col = "purple", lwd=1.5)
legend("topright", c("true density", "fhat with alpha=0.5", "fhat with alpha=1", 
                     "fhat with alpha=2", "complete data"), lty = 1, lwd=1.5, col = c("black","blue","green","red", "purple"), bty = "n")

#dev.off()

#pdf(file = "gamma_eg_n200_log.pdf",   # The directory you want to save the file in
#    width = 4, # The width of the plot in inches
#    height = 4) # The height of the plot in inches


par(mar=rep(0.5,4))
plot(0, 0, type = 'n', xlim = c(0,15), ylim = c(-10, 2))
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)
#rug(res$x)
lines(grid2, log(dgamma(grid2, shape=5, scale=1)), lwd=1.5)
lines(grid2, log(fhat05plot), col = "blue", lwd=1.5)
lines(grid2, log(fhat1plot), col = "green", lwd=1.5)
lines(grid2, log(fhat2plot), col = "red", lwd=1.5)
lines(grid2, log(fhatcompplot), col = "purple", lwd=1.5)
legend("topright", c("true density", "fhat with alpha=0.5", "fhat with alpha=1", 
                     "fhat with alpha=2", "complete data"), lty = 1, lwd=1.5, col = c("black","blue","green","red", "purple"), bty = "n")

#dev.off()

Fhat05plot <- evaluateLogConDens(grid2, fhat05$fhatn, which=5)
Fhat05plot <- Fhat05plot[,6]
Fhat1plot <- evaluateLogConDens(grid2, fhat1$fhatn, which=5)
Fhat1plot <- Fhat1plot[,6]
Fhat2plot <- evaluateLogConDens(grid2, fhat2$fhatn, which=5)
Fhat2plot <- Fhat2plot[,6]
Fhatcompplot <- evaluateLogConDens(grid2, fhat_comp, which=5)
Fhatcompplot <- Fhatcompplot[,6]


par(mar=rep(0.5,4))
plot(0, 0, type = 'n', xlim = c(0,15), ylim = c(0, 1))
abline(v=grid, col=gray(0.9), lty=2, lwd=0.5)
#rug(res$x)
lines(grid2, pgamma(grid2, shape=5, scale=1), lwd=1.5)
lines(grid2, Fhat05plot, col = "blue", lwd=1.5)
lines(grid2, Fhat1plot, col = "green", lwd=1.5)
lines(grid2, Fhat2plot, col = "red", lwd=1.5)
lines(grid2, Fhatcompplot, col = "purple", lwd=1.5)
legend("bottomright", c("true density", "fhat with alpha=0.5", "fhat with alpha=1", 
                        "fhat with alpha=2", "complete data"), lty = 1, lwd=1.5, col = c("black","blue","green","red", "purple"), bty = "n")

#dev.off()
