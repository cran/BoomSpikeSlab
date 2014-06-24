SpikeSlabPriorBase <- function(x = NULL,
                               y = NULL,
                               expected.r2 = .5,
                               prior.df = .01,
                               expected.model.size = 1,
                               optional.coefficient.estimate = NULL,
                               mean.y = mean(y, na.rm = TRUE),
                               sdy = sd(as.numeric(y), na.rm = TRUE),
                               prior.inclusion.probabilities = NULL,
                               number.of.observations = nrow(x),
                               number.of.variables = ncol(x) ) {
  ## Computes information that is shared by the different
  ## implementation of spike and slab priors.  Currently, the only
  ## difference between the different priors is the prior variance on
  ## the regression coefficients.  When that changes, change this
  ## function accordingly, and change all the classes that inherit
  ## from it.

  if (!is.null(x)) {
    stopifnot(is.matrix(x))
  }
  if (!is.null(y) && length(y) != number.of.observations) {
    stop("x and y do not conform")
  }

  ## Compute prior.inclusion.probabilities, the prior inclusion probability
  if (is.null(prior.inclusion.probabilities)) {
    stopifnot(is.numeric(expected.model.size))
    stopifnot(length(expected.model.size) == 1)
    stopifnot(expected.model.size > 0)
    if (expected.model.size < number.of.variables) {
      prior.inclusion.probabilities <-
        rep(expected.model.size/number.of.variables, number.of.variables)
    } else {
      prior.inclusion.probabilities <- rep(1, number.of.variables)
    }
  }
  stopifnot(length(prior.inclusion.probabilities) == number.of.variables)
  stopifnot(all(prior.inclusion.probabilities >= 0))
  stopifnot(all(prior.inclusion.probabilities <= 1))

  ## Compute sigma.guess (guess at the residual standard deviation)
  stopifnot(is.numeric(expected.r2))
  stopifnot(length(expected.r2) == 1)
  stopifnot(expected.r2 < 1 && expected.r2 > 0)
  stopifnot(is.numeric(sdy) && length(sdy) == 1 && sdy > 0)
  sigma.guess <- sqrt((1 - expected.r2)) * sdy

  ## Compute mu, the conditional prior mean of beta (given nonzero)
  if (!is.null(optional.coefficient.estimate)) {
    mu <- optional.coefficient.estimate
    if (length(mu) != number.of.variables) {
      stop("The vector of estimated coefficients has length",
           length(mu),
           "but the design matrix has",
           number.of.variables,
           "columns.")
    }
  } else {
    ## The default behavior is to set the prior mean of all slopes to
    ## zero, and the prior mean of the intercept to ybar.
    mu <- rep(0, number.of.variables)
    mu[1] <- mean.y
  }

  stopifnot(is.numeric(prior.df) && length(prior.df) == 1 && prior.df >= 0)

  ans <- list(prior.inclusion.probabilities = prior.inclusion.probabilities,
              mu = mu,
              sigma.guess = sigma.guess,
              prior.df = prior.df)
  class(ans) <- "SpikeSlabPriorBase"
  return(ans)
}

SpikeSlabPrior <- function(x,
                           y = NULL,
                           expected.r2 = .5,
                           prior.df = .01,
                           expected.model.size = 1,
                           prior.information.weight = .01,
                           diagonal.shrinkage = .5,
                           optional.coefficient.estimate = NULL,
                           max.flips = -1,
                           mean.y = mean(y, na.rm = TRUE),
                           sdy = sd(as.numeric(y), na.rm = TRUE),
                           prior.inclusion.probabilities = NULL) {
  ## Produces a list containing the prior distribution needed for
  ## lm.spike.  This is Zellner's g-prior, where
  ##
  ##    1 / sigsq ~ Ga(prior.weight, prior.sumsq)
  ## beta | gamma ~ N(b, sigsq * V)
  ##        gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## with V^{-1} = prior.information.weight *
  ##             ((1 - diagonal.shrinkage) * xtx / n
  ##               + (diagonal.shrinkage) * diag(diag(xtx / n)))
  ## and prior.sumsq = prior.df * (1 - expected.r2) * sdy^2
  ##
  ## Args:
  ##   x:  The design matrix in the regression problem
  ##   y:  The vector of responses in the regression problem
  ##   expected.r2: A positive scalar less than 1.  The expected
  ##     R-squared statistic in the regression.  Used to obtain the
  ##     prior distribution for the residual variance.
  ##   prior.df: A positive scalar representing the prior 'degrees of
  ##     freedom' for estimating the residual variance.  This can be
  ##     thought of as the amount of weight (expressed as an
  ##     observation count) given to the expected.r2 argument.
  ##   expected.model.size: A positive number less than ncol(x),
  ##     representing a guess at the number of significant predictor
  ##     variables.  Used to obtain the 'spike' portion of the spike
  ##     and slab prior.  The default value leads to pi[i] = .5, which
  ##     is the uniform distribution over the space of models.  This
  ##     is not used if prior.inclusion.probabilities is supplied.
  ##   prior.information.weight: A positive scalar.  Number of
  ##     observations worth of weight that should be given to the
  ##     prior estimate of beta.
  ##   diagonal.shrinkage: The conditionally Gaussian prior for beta
  ##     (the "slab") starts with a precision matrix equal to the
  ##     information in a single observation.  However, this matrix
  ##     might not be full rank.  The matrix can be made full rank by
  ##     averaging with its diagonal.  diagonal.shrinkage is the
  ##     weight given to the diaonal in this average.  Setting this to
  ##     zero gives Zellner's g-prior.
  ##   optional.coefficient.estimate: If desired, an estimate of the
  ##     regression coefficients can be supplied.  In most cases this
  ##     will be a difficult parameter to specify.  If omitted then a
  ##     prior mean of zero will be used for all coordinates except
  ##     the intercept, which will be set to mean(y).
  ##   max.flips: The maximum number of variable inclusion indicators
  ##     the sampler will attempt to sample each iteration.  If negative
  ##     then all indicators will be sampled.
  ##   mean.y:  The mean of the response variable.
  ##   sdy:  The standard deviation of the response variable.
  ##   prior.inclusion.probabilities: A vector giving the prior
  ##     probability of inclusion for each variable.
  ## Returns:
  ##   A list with the values needed to run lm.spike
  ## Details:
  ## Let gamma be a logical vector with length ncol(x), where
  ##   gamma[i] == TRUE means beta[i] != 0.  The prior is:
  ##
  ## Pr(beta[i] != 0) = pi
  ## p(beta[gamma] | x, sigma)
  ##     = N(b[gamma], sigma^2 * Omega[gamma, gamma])
  ## p(1/sigma^2) = Gamma(prior.df / 2, prior.sum.of.squares / 2)
  ##
  ## The parameters are set as follows:
  ##   prior.inclusion.probabilities = expected.model.size / ncol(x)
  ##   b is a vector of zeros, except for the intercept (first term),
  ##     which is set to mean(y)
  ##   Omega^{-1} = [(1 - diagonal.srhinkage) * x^Tx +
  ##     diagonal.shrinkage * diag(x^Tx)]
  ##     * prior.information.weight / n
  ##   prior.sum.of.squares = var(y) * (1 - expected.r2) * prior.df

  ans <- SpikeSlabPriorBase(
      x = x,
      y = y,
      expected.r2 = expected.r2,
      prior.df = prior.df,
      expected.model.size = expected.model.size,
      optional.coefficient.estimate = optional.coefficient.estimate,
      mean.y = mean.y,
      sdy = sdy,
      prior.inclusion.probabilities = prior.inclusion.probabilities)

  ans$max.flips <- max.flips

  p <- length(ans$mu)
  n <- nrow(x)

  ## Compute siginv, solve(siginv) times the residual variance is the
  ## prior variance of beta, conditional on nonzero elements.
  w <- diagonal.shrinkage
  if (w < 0 | w > 1) {
    stop("Illegal value of diagonal shrinkage: ",
         w,
         "(must be between 0 and 1)")
  }
  xtx <- crossprod(x) / n
  if (p == 1) {
    d <- xtx
  } else {
    d <- diag(diag(xtx))
  }
  xtx <- w*d + (1 - w)*xtx
  xtx <- xtx * prior.information.weight

  ans$siginv <- xtx
  class(ans) <- c("SpikeSlabPrior", class(ans))
  return(ans)
}

IndependentSpikeSlabPrior <- function(x = NULL,
                                      y = NULL,
                                      expected.r2 = .5,
                                      prior.df = .01,
                                      expected.model.size = 1,
                                      prior.beta.sd = NULL,
                                      optional.coefficient.estimate = NULL,
                                      mean.y = mean(y, na.rm = TRUE),
                                      sdy = sd(as.numeric(y), na.rm = TRUE),
                                      sdx = apply(as.matrix(x), 2, sd,
                                        na.rm = TRUE),
                                      prior.inclusion.probabilities = NULL,
                                      number.of.observations = nrow(x),
                                      number.of.variables = ncol(x)) {
  ## This version of the spike and slab prior assumes
  ##
  ##       1.0 / sigsq ~ Ga(prior.weight, prior.sumsq),
  ##      beta | gamma ~ N(b, V),
  ##             gamma ~ Independent Bernoulli(prior.inclusion.probabilities)
  ##
  ## where V is a diagonal matrix independent of the residual
  ## variance sigsq.
  ##
  ## Args:
  ##   See SpikeSlabPriorBase.
  ##   Additional args:
  ##   prior.beta.sd: A vector of positive numbers giving the prior
  ##     standard deviation of each model coefficient, conditionl on
  ##     inclusion.  If NULL it will be set to 10 * the ratio of sdy /
  ##     sdx.
  ans <- SpikeSlabPriorBase(
             x = x,
             y = y,
             expected.r2 = expected.r2,
             prior.df = prior.df,
             expected.model.size = expected.model.size,
             optional.coefficient.estimate = optional.coefficient.estimate,
             mean.y = mean.y,
             sdy = sdy,
             prior.inclusion.probabilities = prior.inclusion.probabilities,
             number.of.observations = number.of.observations,
             number.of.variables = number.of.variables)

  ## If any columns of x have zero variance (e.g. the intercept) then
  ## don't normalize by their standard deviation.
  sdx[sdx == 0] <- 1
  stopifnot(length(sdy) == 1)
  stopifnot(is.numeric(sdy))
  stopifnot(sdy > 0)
  stopifnot(is.numeric(sdx))
  stopifnot(all(sdx > 0))

  ## Each 1 sd change in x could have up to a 10 sd change in y.  That's BIG
  if (is.null(prior.beta.sd)) {
    prior.beta.sd <- 10 * sdy / sdx
  }
  ans$prior.variance.diagonal <- prior.beta.sd^2
  class(ans) <- c("IndependentSpikeSlabPrior", class(ans))
  return(ans)
}
