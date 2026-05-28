\name{Eval}
\alias{Eval}
\title{L2 error calculations}
\usage{
Eval(object, pdf, length_grid = 1001)
}
\arguments{
\item{object}{A fitted object returned by \code{MALC()}, \code{BK2002()},
\code{BinnedNP()}, or \code{KernSmooth()}.}
\item{pdf}{True density function used for comparison.}
\item{length_grid}{Number of evaluation grid points.}
}
\description{
Compute numerical L2 error on a dense evaluation grid, matching the simulation
study convention. The fitted object supplies the density estimator; \code{pdf}
supplies the target density.
}
\value{A numeric L2 error.}
