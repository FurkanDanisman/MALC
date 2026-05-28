\name{MALC}
\alias{MALC}
\alias{dmalc}
\alias{pmalc}
\alias{qmalc}
\alias{mean_malc}
\alias{Plot}
\alias{plot.MALC}
\title{Fit and evaluate MALC}
\usage{
MALC(counts, grid, smooth = FALSE, B = 10000, alpha = 2,
  seed = 20180621, log_conc = TRUE)

dmalc(object, x)

pmalc(object, q, length_grid = 1001)

qmalc(object, p, length_grid = 1001)

mean_malc(object)

Plot(object, ...)

\method{plot}{MALC}(x, eval_grid = NULL, add = FALSE,
  xlab = "x", ylab = "Density", col = "navy", lwd = 2, ...)
}
\arguments{
\item{counts}{Integer observed frequencies for each bin.}
\item{grid}{Numeric bin breaks. Must have length \code{length(counts) + 1}.}
\item{smooth}{Logical; store the fitted object as smoothed when \code{TRUE}.}
\item{B}{Number of Monte Carlo draws used internally by MALC.}
\item{alpha}{Shape parameter used in the beta perturbation.}
\item{seed}{Random seed used internally.}
\item{log_conc}{Passed to \code{logcondens::logConDens}.}
\item{object}{A fitted \code{MALC} object.}
\item{x}{Evaluation points, or a fitted object for \code{plot}.}
\item{q}{Quantile locations for the CDF.}
\item{p}{Probabilities for the quantile function.}
\item{length_grid}{Number of points in the numerical evaluation grid.}
\item{eval_grid}{Optional grid used for plotting.}
\item{add}{Logical; add to an existing plot when \code{TRUE}.}
\item{xlab, ylab, col, lwd, ...}{Graphical arguments.}
}
\value{
\code{MALC()} returns a fitted \code{MALC} object. Evaluation functions
return numeric vectors. \code{mean_malc()} returns the EM-step mean estimate
stored in the fitted object.
}
\description{
Fit MALC from observed bin frequencies and evaluate density, distribution,
quantile, EM-step mean, and plot methods.
}
