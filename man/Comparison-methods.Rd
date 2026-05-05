\name{Comparison-methods}
\alias{BK2002}
\alias{BinnedNP}
\alias{KernSmooth}
\alias{Plot_Kernel}
\title{Comparison density estimators}
\usage{
BK2002(counts, grid, h = NULL, bw_B = 50, bw_seed = 1,
  bw_method = c("ste", "dpi"), max_iter = 200, tol_L1 = 1e-6,
  verbose = FALSE)

BinnedNP(counts, grid, bandwidth = c("boot", "plugin"),
  seed_start = 1L, max_tries = 200L)

KernSmooth(counts, grid, seed = 1, scalest = "minim", level = 2L,
  kernel = "normal", canonical = FALSE, truncate = TRUE)

Plot_Kernel(object, eval_grid = NULL, add = FALSE, xlab = "x",
  ylab = "Density", col = "darkgreen", lwd = 2, ...)
}
\arguments{
\item{counts}{Integer observed frequencies for each bin.}
\item{grid}{Numeric bin breaks. Must have length \code{length(counts) + 1}.}
\item{h}{Optional bandwidth for BK2002.}
\item{bw_B, bw_seed, bw_method}{Bandwidth-selection controls for BK2002.}
\item{max_iter, tol_L1, verbose}{Fixed-point iteration controls for BK2002.}
\item{bandwidth}{BinnedNP bandwidth choice stored in the fitted object.}
\item{seed_start, max_tries}{Bandwidth retry controls for BinnedNP.}
\item{seed, scalest, level, kernel, canonical, truncate}{Controls passed to KernSmooth bandwidth and density routines.}
\item{object}{A fitted comparison-method object.}
\item{eval_grid}{Optional grid used for plotting the fitted density.}
\item{add}{If \code{TRUE}, add the density to the current plot.}
\item{xlab, ylab, col, lwd, ...}{Graphical controls passed to \code{plot()} or \code{lines()}.}
}
\description{
Fit the three baseline estimators used in the DensOLog simulation study:
BK2002, BinnedNP, and KernSmooth. All methods use observed frequencies and
grid breaks as inputs. \code{Plot_Kernel()} plots any fitted comparison-method
object.
}
\value{
Fitting functions return method-specific fitted objects. \code{Plot_Kernel()}
returns invisibly a data frame with the plotted grid and density values.
}
