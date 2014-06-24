# Copyright 2010 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

logit.spike <- function(formula, niter, data, subset, prior = NULL,
                        na.action = options("na.action"), contrasts = NULL,
                        drop.unused.levels = TRUE,
                        initial.value = NULL, ping = niter / 10, nthreads = 4,
                        clt.threshold = 2,
                        mh.chunk.size = 10, proposal.df = 3,
                        seed = NULL,
                        ...) {
  ## Uses Bayesian MCMC to fit a logistic regression model with a
  ## spike-and-slab prior.
  ##
  ## Args:
  ##   formula: model formula, as would be passed to 'glm', specifying
  ##     the maximal model (i.e. the model with all predictors
  ##     included).
  ##   niter:  desired number of MCMC iterations
  ##   ping: if positive, then print a status update every 'ping' MCMC
  ##     iterations.
  ##   nthreads:  The number of threads to use when imputing latent data.
  ##   data:  optional data.frame containing the data described in 'formula'
  ##   subset: an optional vector specifying a subset of observations
  ##     to be used in the fitting process.
  ##   prior: an optional list such as that returned from
  ##     SpikeSlabPrior.  If missing, SpikeSlabPrior
  ##     will be called with the remaining arguments.
  ##   na.action: a function which indicates what should happen when
  ##     the data contain ‘NA’s.  The default is set by the
  ##     ‘na.action’ setting of ‘options’, and is ‘na.fail’ if that is
  ##     unset.  The ‘factory-fresh’ default is ‘na.omit’.  Another
  ##     possible value is ‘NULL’, no action.  Value ‘na.exclude’ can
  ##     be useful.
  ##   contrasts: an optional list. See the ‘contrasts.arg’ of
  ##     ‘model.matrix.default’.  An optional list.
  ##   drop.unused.levels: should factor levels that are unobserved be
  ##     dropped from the model?
  ##   initial.value: Initial value of logistic regression
  ##     coefficients for the MCMC algorithm.  Can be given as a
  ##     numeric vector, a 'logit.spike' object, or a 'glm' object.
  ##     If a 'logit.spike' object is used for initialization, it is
  ##     assumed to be a previous MCMC run to which 'niter' futher
  ##     iterations should be added.  If a 'glm' object is supplied,
  ##     its coefficients will be used as the initial values in the
  ##     MCMC simulation.
  ##   clt.threshold: The smallest number of successes or failures
  ##     needed to do use asymptotic data augmentation.
  ##   mh.chunk.size: The largest number of parameters to draw
  ##     together in a single Metropolis-Hastings proposal.  A
  ##     non-positive number means to use a single chunk.
  ##   proposal.df: The degrees of freedom parameter for the
  ##     multivariate T proposal distribution used for
  ##     Metropolis-Hastings updates.  A nonpositive number means to
  ##     use a Gaussian proposal.
  ##   seed: Seed to use for the C++ random number generator.  NULL or
  ##     an int.  If NULL, then the seed will be taken from the global
  ##     .Random.seed object.
  ##   ... : parameters to be passed to SpikeSlabPrior
  ##
  ## Returns:
  ##   An object of class 'logit.spike', which is a list containing the
  ##   following values
  ##   beta: A 'niter' by 'ncol(X)' matrix of regression coefficients
  ##     many of which may be zero.  Each row corresponds to an MCMC
  ##     iteration.
  ##   prior:  The prior that was used to fit the model.
  ##  In addition, the returned object contains sufficient details for
  ##  the call to model.matrix in the predict.lm.spike method.

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")

  if (!is.null(dim(y)) && length(dim(y)) > 1) {
    stopifnot(length(dim(y)) == 2, ncol(y) == 2)
    ## If the user passed a formula like "cbind(successes, failures) ~
    ## x", then y will be a two column matrix
    ny <- y[, 1] + y[, 2]
    y <- y[, 1]
  } else {
    ## The following line admits y's which are TRUE/FALSE, 0/1 or 1/-1.
    y <- y > 0
    ny <- rep(1, length(y))
  }

  x <- model.matrix(mt, mf, contrasts)
  if (is.null(prior)) {
    prior <- SpikeSlabPrior(x, y, ...)
  }

  if (!is.null(initial.value)) {
    if (inherits(initial.value, "logit.spike")) {
      stopifnot(colnames(initial.value$beta) == colnames(x))
      beta0 <- as.numeric(tail(initial.value$beta, 1))
    } else if (inherits(initial.value, "glm")) {
      stopifnot(colnames(initial.value$beta) == colnames(x))
      beta0 <- coef(initial.value)
    } else if (is.numeric(initial.value)) {
      stopifnot(length(initial.value) == ncol(x))
      beta0 <- initial.value
    } else {
      stop("initial.value must be a 'logit.spike' object, a 'glm' object,",
           "or a numeric vector")
    }
  } else {
    ## No initial value was supplied
    beta0 <- rep(0, ncol(x))
    if (all(x[, 1] == 1)) {
      p.hat <- sum(y) / sum(ny)
      beta0[1] <- qlogis(p.hat)
    }
  }

  ans <- .logit.spike.fit(x, y, ny, prior, niter, ping, nthreads, beta0,
                          mh.chunk.size, clt.threshold, seed)

  ## The stuff below will be needed by predict.lm.spike.
  ans$contrasts <- attr(x, "contrasts")
  ans$xlevels <- .getXlevels(mt, mf)
  ans$call <- cl
  ans$terms <- mt

  if (!is.null(initial.value) && inherits(initial.value, "logit.spike")) {
    ans$beta <- rbind(initial.value$beta, ans$beta)
  }

  ## Make the answer a class, so that the right methods will be used.
  class(ans) <- c("logit.spike", "lm.spike")
  return(ans)
}

predict.logit.spike <- function(object, newdata, burn = 0,
                                type = c("prob", "logit", "link", "response"),
                                na.action = na.pass, ...) {
  ## Prediction method for logit.spike
  ## Args:
  ##   object: object of class "logit.spike" returned from the logit.spike
  ##     function
  ##   newdata: A data frame including variables with the same names
  ##     as the data frame used to fit 'object'.
  ##   burn: The number of MCMC iterations in 'object' that should be
  ##     discarded.  If burn < 0 then all iterations are kept.
  ##   type: The type of prediction desired.  If 'prob' then the
  ##     prediction is returned on the probability scale.  If 'logit'
  ##     then it is returned on the logit scale (i.e. the scale of the
  ##     linear predictor).  Also accepts 'link' and 'response' for
  ##     compatibility with predict.glm.
  ##   ...: unused, but present for compatibility with generic predict().
  ## Returns:
  ##   A matrix of predictions, with each row corresponding to a row
  ##   in newdata, and each column to an MCMC iteration.

  type <- match.arg(type)

  tt <- terms(object)
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
  X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  if (nrow(X) != nrow(newdata)) {
    warning("Some entries in newdata have missing values, and  will",
            "be omitted from the prediction.")
  }

  beta <- object$beta
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }

  eta <- X %*% t(beta)
  if (type == "logit" || type == "link") return(eta)
  if (type == "prob" || type == "response") return(plogis(eta))
}

.logit.spike.fit <- function(x, y, ny, prior, niter, ping, nthreads, beta0,
                             mh.chunk.size, clt.threshold, seed) {
  ## Driver function for lm.spike.
  ## Args:
  ##   x: design matrix with 'n' rows corresponding to observations and
  ##     'p' columns corresponding to predictor variables.
  ##   y: vector of integer responses (success counts) of length n
  ##   ny: vector of integer trial counts of length n
  ##   prior: a list structured like the return value from
  ##     SpikeSlabPrior
  ##   niter:  the number of desired MCMC iterations
  ##   ping:  frequency with which to print MCMC status updates
  ##   nthreads:  number of threads to use when imputing latent data
  ##   beta0:  The initial value in the MCMC simulation.
  ##   seed: Seed to use for the C++ random number generator.  NULL or
  ##     an int.
  ##
  ## Returns:
  ##   An object of class logit.spike, which is a list with the elements
  ##   described in the 'logit.spike' function.
  ##
  ## TODO(stevescott): move this to a namespace so it can't be called
  ## without logit.spike

  stopifnot(nrow(x) == length(y),
            length(prior$mu) == ncol(x),
            length(prior$prior.inclusion.probabilities) == ncol(x),
            all(ny >= y),
            all(y >= 0))

  nobs <- nrow(x)
  p <- ncol(x)
  beta.draws <- matrix(0, nrow=niter, ncol = p)

  if (is.null(prior$max.flips)) {
    prior$max.flips <- -1
  }

  if (is.null(seed)) {
    tmp <- rnorm(1)  ## ensure that .Random.seed exists
    seed <- .Random.seed[1]
  }

  ans <- .C(logit_spike_slab_wrapper,
            as.double(x),
            as.integer(y),
            as.integer(ny),
            as.double(prior$mu),
            as.double(prior$siginv),
            as.double(prior$prior.inclusion.probabilities),
            as.integer(prior$max.flips),
            as.integer(nobs),
            as.integer(p),
            as.integer(niter),
            as.integer(ping),
            as.integer(nthreads),
            as.double(beta0),
            as.integer(clt.threshold),
            as.integer(mh.chunk.size),
            as.integer(seed),
            beta = as.double(beta.draws))[c("beta")]
  ans$beta <- matrix(ans$beta, nrow = niter, byrow = FALSE)
  variable.names <- dimnames(x)[[2]]
  if (!is.null(variable.names)) {
    dimnames(ans$beta)[[2]] <- variable.names
  }
  ans$prior <- prior
  class(ans) <- c("logit.spike", "lm.spike")
  return(ans)
}

summary.logit.spike <- function(object, burn = 0, order = TRUE, ...) {
  ## Summary method for logit.spike coefficients
  ## Args:
  ##   object:  an object of class 'logit.spike'
  ##   burn: an integer giving the number of MCMC iterations to
  ##     discard as burn-in
  ##   order: Logical indicating whether the output should be ordered
  ##     according to posterior inclusion probabilities
  ## Returns:
  ## An object of class 'summary.logit.spike' that summarizes the
  ## model coefficients as in SummarizeSpikeSlabCoefficients.
  ans <- SummarizeSpikeSlabCoefficients(object$beta, burn, order)
  class(ans) <- "summary.logit.spike"
  return(ans)
}

print.summary.logit.spike <- function(x, ...) {
  ## print method for summary.logit.spike objects.
  print.default(signif(x, 3))
}
