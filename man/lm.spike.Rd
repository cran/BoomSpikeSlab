% Copyright 2010 Google Inc. All Rights Reserved.
% Author: stevescott@google.com (Steve Scott)

\name{lm.spike}
\Rdversion{1.1}
\alias{lm.spike}
\alias{spikeslab}
\title{
  Spike and slab regression
}
\description{
  MCMC algorithm for linear regression models with a 'spike-and-slab'
  prior that places some amount of posterior probability at zero for a
  subset of the regression coefficients.
}
\usage{
lm.spike(formula,
         niter,
         data,
         subset,
         prior = NULL,
         contrasts = NULL,
         drop.unused.levels = TRUE,
         method = c("SSVS", "DA"),
         ping = niter / 10,
         seed = NULL,
         ...)
}
\arguments{
  \item{formula}{
    formula for the maximal model (with all variables included), this is
    parsed the same way as a call to \code{lm}.
  }
  \item{niter}{
    The number of MCMC iterations to run.  Be sure to include enough so
    you can throw away a burn-in set.
  }
  \item{data}{
    an optional data frame, list or environment (or object coercible by
    'as.data.frame' to a data frame) containing the variables in the
    model.  If not found in 'data', the variables are taken from
    'environment(formula)', typically the environment from which
    'lm.spike' is called..
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
  }

  \item{prior}{ An optional list returned by
    \code{\link{SpikeSlabPrior}}.  If \code{prior} is missing then a
    default prior will be used.  See \code{\link{SpikeSlabPrior}}.
  }

  \item{contrasts}{ An optional list. See the \code{contrasts.arg}
    argument of \code{\link{model.matrix.default}}.  }

  \item{drop.unused.levels}{Logical indicating whether unobserved factor
    levels should be dropped from the model.}

  \item{method}{The MCMC method to use.  SSVS is the stochastic search
    variable selection algorithm from George and McCulloch (1998).  DA
    is the data augmentation method from Clyde and Ghosh (2011).}

  \item{ping}{The frequency with which to print status update messages
    to the screen.  For example, if \code{ping == 10} then an update
    will be printed every 10 MCMC iterations.}

  \item{seed}{An integer to use as the random seed for the underlying
    C++ code.  If \code{NULL} then the seed will be set using the
    clock.}

  \item{\dots}{
    Extra arguments to be passed to \code{\link{SpikeSlabPrior}} (if
    \code{method == "SSVS"}) or \code{\link{IndependentSpikeSlabPrior}}
    (if \code{method == "DA"}).
}

}

\details{

  Both methods ("DA" and "SSVS") draw each variable inclusion indicator
  given all others, in a Gibbs sampler.  The DA method includes an extra
  data augmentation step that renders each indicator conditionally
  independent of the others given the latent data.  There is residual
  dependence between successive MCMC steps introduced by the latent
  data, but the paper by Ghosh and Clyde suggested that on balance
  mixing should be improved.

  Regarding the overall compute time, the DA method decomposes the
  (potentially very large) model matrix one time, at the start of the
  algorithm.  But it then works with independent scalar updates.  The
  SSVS algorithm does not have the upfront cost, but it works with many
  small matrix decompositions each MCMC iteration.  The DA algorithm is
  very likely to be faster in terms of time per iteration.

  Finally, note that the two algorithms require slightly different
  priors.  The DA algorithm requires a priori independence, while the
  SSVS algorithm can work with arbitrary conjugate priors.

}

\value{
  Returns an object of class \code{lm.spike}, which is a list with the
  following elements
  \item{beta}{A \code{niter} by \code{ncol(x)} matrix of regression
    coefficients, many of which may be zero.  Each row corresponds to an
    MCMC iteration. }

  \item{sigma}{ A vector of length \code{niter} containing the MCMC
    draws of the residual standard deviation parameter. }

  \item{prior}{The prior used to fit the model.  If a \code{prior} was
    supplied as an argument it will be returned.  Otherwise this will be
    the automatically generated prior based on the other function
    arguments. }

}

\references{
  George and McCulloch (1997), "Approaches to Bayesian Variable
  Selection", \emph{Statistica Sinica}, \bold{7}, 339 -- 373.
  \url{http://www3.stat.sinica.edu.tw/statistica/oldpdf/A7n26.pdf}

  Ghosh and Clyde (2011) "Rao-Blackwellization for Bayesian variable
  selection and model averaging in linear and binary regression: A novel
  data augmentation approach", \emph{Journal of the American Statistical
  Association}, \bold{106} 1041-1052.
  \url{http://homepage.stat.uiowa.edu/~jghsh/ghosh_clyde_2011_jasa.pdf}
}

\author{
  Steven L. Scott
}

\seealso{
  \code{\link{SpikeSlabPrior}},
  \code{\link{plot.lm.spike}},
  \code{\link{summary.lm.spike}},
  \code{\link{predict.lm.spike}}.
}

\examples{
  n <- 100
  p <- 10
  ngood <- 3
  niter <- 1000
  sigma <- .8

  x <- cbind(1, matrix(rnorm(n * (p-1)), nrow=n))
  beta <- c(rnorm(ngood), rep(0, p - ngood))
  y <- rnorm(n, x \%*\% beta, sigma)
  x <- x[,-1]
  model <- lm.spike(y ~ x, niter=niter)
  plot.ts(model$beta)
  opar <-  par(ask=TRUE)
  on.exit(par(opar))
  hist(model$sigma)  ## should be near 8
  plot(model)
  summary(model)
}
\keyword{models}
\keyword{regression}
