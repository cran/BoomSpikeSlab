\name{spike.slab.prior.base}
\Rdversion{1.1}
\alias{SpikeSlabPriorBase}
\title{Base class for spike and slab priors}

\description{
  A base class for SpikeSlabPrior and SpikeSlabPriorBase to ensure that
  elements common to both classes are handled consistently.  Users
  will not normally interact with this function.
}

\usage{
SpikeSlabPriorBase(x = NULL,
                   y = NULL,
                   expected.r2 = .5,
                   prior.df = .01,
                   expected.model.size = 1,
                   optional.coefficient.estimate = NULL,
                   mean.y = mean(y, na.rm = TRUE),
                   sdy = sd(as.numeric(y), na.rm = TRUE),
                   prior.inclusion.probabilities = NULL,
                   number.of.observations = nrow(x),
                   number.of.variables = ncol(x) )
}

\arguments{
  \item{x}{
    The design matrix for the regression problem.  Missing data is not allowed.
  }

  \item{y}{
    The vector of responses for the regression.  Missing data is not allowed.
  }

  \item{expected.r2}{
    The expected R-square for the regression.  The spike and slab prior
    requires an inverse gamma prior on the residual variance of the
    regression.  The prior can be parameterized in terms of a guess at
    the residual variance, and a "degrees of freedom" representing the
    number of observations that the guess should weigh.  The guess at
    sigma^2 is set to \code{ (1-expected.r2) * var(y) }.
  }

  \item{prior.df}{
    A positive scalar representing the prior 'degrees of freedom' for
    estimating the residual variance.  This can be thought of as the
    amount of weight (expressed as an observation count) given to the
    \code{expected.r2} argument.
  }

  \item{expected.model.size}{
    A positive number less than \code{ncol(x)}, representing a guess at
    the number of significant predictor p variables.  Used to obtain the
    'spike' portion of the spike and slab prior.  The default value
    leads to \code{pi[i] = .5}, which is the uniform distribution over
    the space of models.
  }

  \item{optional.coefficient.estimate}{ If desired, an estimate of the
    regression coefficients can be supplied.  In most cases this will be
    a difficult parameter to specify.  If omitted then a prior mean of
    zero will be used for all coordinates except the intercept, which
    will be set to mean(y).  }

  \item{mean.y}{The mean of the response vector, for use in cases when
    specifying the response vector is undesirable.}

  \item{sdy}{The standard deviation of the response vector, for use in
    cases when specifying the response vector is undesirable.}

  \item{prior.inclusion.probabilities}{A vector giving the prior
    probability of inclusion for each variable.}

  \item{number.of.observations}{The number of rows in \code{x} and
    \code{y}.}

  \item{number.of.variables}{The number of columns in \code{x}.}

}

\author{
  Steven L. Scott
}
