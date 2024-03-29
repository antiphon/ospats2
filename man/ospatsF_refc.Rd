% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ospatsF-refactored-c.R
\name{ospatsF_refc}
\alias{ospatsF_refc}
\title{Legacy Ospats Algorithm}
\usage{
ospatsF_refc(
  x,
  covmodel_range,
  nstrata,
  niter = 100,
  niter_outer = 100,
  verbose = 0,
  temperature = 1,
  coolingrate = 0.95,
  rsquared = 1,
  Cov
)
}
\arguments{
\item{x}{data, with variables \link{x y pred var}}

\item{covmodel_range}{Range parameter for assumed exponential correlation model}

\item{nstrata}{number of starta to consider}

\item{niter}{Number of iterations per one run of the optimisation algorithm}

\item{niter_outer}{Number of independent runs of the optimisation algorithm}

\item{verbose}{Print runtime diagnostics?}

\item{temperature}{Annealing factor, will accept slightly bad moves with prob ~exp(-abs(delta)/temperature)}

\item{coolingrate}{Change temperature each interation by this factor. Should be at most 1.}

\item{rsquared}{The R^2 in the paper (default: 1)}

\item{Cov}{Optional, overrides the covariance matrix calculation using exp-correlation. No checks with data variances.}
}
\description{
Refactored Ospats::ospatsF and using RCpp
}
\details{
The algorithm in \link{ospatsF} with new inputs and outputs. Only change is that
the covariance matrix can be given to avoid the exponential covariance assumption and the "mean" covariance formula
used in the original scripts.

See the vignette \code{vignette("legacy",)}
}
