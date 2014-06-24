% Copyright 2011 Google Inc. All Rights Reserved.
% Author: stevescott@google.com (Steve Scott)

\name{PlotModelSize}
\Rdversion{1.1}
\alias{PlotModelSize}
\title{
  Plot a distribution of model size
}
\description{
  Produces a histogram of number of nonzero coefficients in a
  spike-and-slab regression.
}

\usage{
  PlotModelSize(model, burn = 0, xlab= "Number of nonzero coefficients", ...)
}
\arguments{
  \item{model}{
    An object of class \code{lm.spike}.
  }

  \item{burn}{ The number of MCMC iterations to be discarded as burn-in.
    }

  \item{xlab}{ Label for the horizontal axis.}

  \item{\dots}{
    Additional arguments to be passed to \code{\link{hist}}
  }

}
\value{
  Invisibly returns the vector of MCMC draws of model sizes.
}

\author{
  Steven L. Scott
}

\seealso{
  \code{\link{lm.spike}}
  \code{\link{plot.lm.spike}}
}

\examples{
simulate.lm.spike <- function(n = 100, p = 10, ngood = 3, niter=1000, sigma = 8){
  x <- cbind(matrix(rnorm(n * (p-1)), nrow=n))
  beta <- c(rnorm(ngood), rep(0, p - ngood))
  y <- rnorm(n, beta[1] + x \%*\% beta[-1], sigma)
  draws <- lm.spike(y ~ x, niter=niter)
  return(invisible(draws))
}
model <- simulate.lm.spike(n = 1000, p = 50, sigma = .3)
PlotModelSize(model, burn = 10)
}
\keyword{models}
\keyword{regression}
