\name{DensOLog}
\alias{DensOLog}
\alias{ddensolog}
\alias{pdensolog}
\alias{qdensolog}
\alias{mean_densolog}
\alias{evaluate_DensOLog}
\alias{Plot}
\alias{plot.DensOLog}
\title{Fit and evaluate DensOLog}
\usage{
DensOLog(counts, grid, smooth = FALSE, B = 10000, alpha = 2,
  seed = 20180621, log_conc = TRUE)

ddensolog(object, x, smooth = NULL)

pdensolog(object, q, smooth = NULL, length_grid = 1001)

qdensolog(object, p, smooth = NULL, length_grid = 1001)

mean_densolog(object)

evaluate_DensOLog(object, x, smooth = NULL)

Plot(object, ...)

\method{plot}{DensOLog}(x, smooth = NULL, eval_grid = NULL, add = FALSE,
  xlab = "x", ylab = "Density", col = "navy", lwd = 2, ...)
}
\arguments{
\item{counts}{Integer observed frequencies for each bin.}
\item{grid}{Numeric bin breaks. Must have length \code{length(counts) + 1}.}
\item{smooth}{Logical; use the smoothed DensOLog estimate when \code{TRUE}.}
\item{B}{Number of Monte Carlo draws used internally by DensOLog.}
\item{alpha}{Shape parameter used in the beta perturbation.}
\item{seed}{Random seed used internally.}
\item{log_conc}{Passed to \code{logcondens::logConDens}.}
\item{object}{A fitted \code{DensOLog} object.}
\item{x}{Evaluation points, or a fitted object for \code{plot}.}
\item{q}{Quantile locations for the CDF.}
\item{p}{Probabilities for the quantile function.}
\item{length_grid}{Number of points in the numerical evaluation grid.}
\item{eval_grid}{Optional grid used for plotting.}
\item{add}{Logical; add to an existing plot when \code{TRUE}.}
\item{xlab, ylab, col, lwd, ...}{Graphical arguments.}
}
\value{
\code{DensOLog()} returns a fitted \code{DensOLog} object. Evaluation functions
return numeric vectors. \code{mean_densolog()} returns the EM-step mean estimate
stored in the fitted object.
}
\description{
Fit DensOLog from observed bin frequencies and evaluate density, distribution,
quantile, EM-step mean, and plot methods.
}
