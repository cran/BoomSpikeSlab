% Copyright 2010 Google Inc. All Rights Reserved.
% Author: stevescott@google.com (Steve Scott)

\name{summary.logit.spike}
\Rdversion{1.1}
\alias{summary.logit.spike}
\title{
  Numerical summaries of the results from a spike and slab logistic regression.
}
\description{
  Produces a summary of the marginal distribution of model coefficients
  from a spike and slab logistic regression.
}

\usage{
\method{summary}{logit.spike}(object, burn = 0, order = TRUE, ...)
}
\arguments{
  \item{object}{
    An object of class \code{logit.spike}.
  }

  \item{burn}{
    The number of MCMC iterations in the ojbect to be discarded as
    burn-in.
  }

  \item{order}{ Logical.  If \code{TRUE} then the coefficients are
    presented in order of their posterior inclusion probabilities.
    Otherwise the order of the coefficients is the same as in
    \code{object}.  }

  \item{\dots}{
    Unused.  Present for compatibility with generic summary().
  }

}

\value{

 A five-column matrix with rows representing model coefficients.  The
 first two columns are the posterior mean and standard deviation of each
 coefficient, including the point mass at zero.  The next two columns
 are the posterior mean and standard deviations conditional on the
 coefficient being nonzero.  The last column is the probability of a
 nonzero coefficient.
}

\author{
  Steven L. Scott
}

\seealso{
  \code{\link{logit.spike}}
  \code{\link{SpikeSlabPrior}}
}

\examples{

  n <- 100
  p <- 10
  ngood <- 3
  niter <- 1000

  x <- cbind(1, matrix(rnorm(n * (p-1)), nrow=n))
  beta <- c(rnorm(ngood), rep(0, p - ngood))
  prob <- plogis(x \%*\% beta)
  y <- runif(n) < prob
  x <- x[,-1]
  model <- logit.spike(y ~ x, niter=niter)
  summary(model)
}
\keyword{ models, regression}
