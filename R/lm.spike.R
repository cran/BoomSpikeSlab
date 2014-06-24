# Copyright 2010 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

lm.spike <- function(formula,
                     niter,
                     data,
                     subset,
                     prior = NULL,
                     contrasts = NULL,
                     drop.unused.levels = TRUE,
                     method = c("SSVS", "DA"),
                     ping = niter / 10,
                     seed = NULL,
                     ...) {
  ## Uses Bayesian MCMC to fit a linear regression model with a
  ## spike-and-slab prior.
  ##
  ## Args:
  ##   formula: Model formula, as would be passed to 'lm', specifying
  ##     the maximal model (i.e. the model with all predictors
  ##     included).
  ##   niter:  Desired number of MCMC iterations
  ##   data:  Optional data.frame containing the data described in 'formula'
  ##   subset: An optional vector specifying a subset of observations
  ##     to be used in the fitting process.
  ##   prior: An object of class SpikeSlabPrior or
  ##     IndependentSpikeSlabPrior.  If missing, SpikeSlabPrior will
  ##     be called with the arguments passed as ....
  ##   contrasts: An optional list. See the 'contrasts.arg' of
  ##     ‘model.matrix.default’.
  ##   drop.unused.levels: logical.  Should unobserved factor levels
  ##     be dropped from the model?
  ##   method: The MCMC method to use for posterior sampling.  SSVS is
  ##     the original "stochastic search variable selection" from
  ##     George and McCulloch 1998.  DA is the method from Clyde and
  ##     Ghosh (2011).
  ##   ping: Write a status update to the console every 'ping'
  ##     iterations.
  ##   seed:  An integer to use for the C++ seed.
  ##   ... : Parameters to be passed to SpikeSlabPrior or
  ##     IndependentSpikeSlabPrior., if 'prior' is not specified
  ##     directly.
  ##
  ## Returns:
  ##   An object of class 'lm.spike', which is a list containing the
  ##   following values
  ##   beta: A 'niter' by 'ncol(X)' matrix of regression coefficients
  ##     many of which may be zero.  Each row corresponds to an MCMC
  ##     iteration.
  ##   sigma: A vector of length 'niter' containing the MCMC draws of
  ##     the residual standard deviation parameter.
  ##   prior:  The prior that was used to fit the model.
  ## In addition, the returned object contains sufficient details for
  ## the call to model.matrix in the predict.lm.spike method.

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")

  x <- model.matrix(mt, mf, contrasts)
  method <- match.arg(method)
  if (is.null(prior)) {
    if (method == "SSVS") {
      prior <- SpikeSlabPrior(x, y, ...)
    } else if (method == "DA") {
      prior <- IndependentSpikeSlabPrior(x, y, ...)
    }
  }
  stopifnot(inherits(prior, "SpikeSlabPriorBase"))

  stopifnot(is.numeric(ping))
  stopifnot(length(ping) == 1)

  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }

  ans <- .lm.spike.fit(x, y, prior, as.integer(niter), as.integer(ping), seed)

  ## The stuff below will be needed by predict.lm.spike.
  ans$contrasts <- attr(x, "contrasts")
  ans$xlevels <- .getXlevels(mt, mf)
  ans$call <- cl
  ans$terms <- mt
  ans$sample.sd <- sd(y)

  ## Make the answer a class, so that the right methods will be used.
  class(ans) <- "lm.spike"
  return(ans)
}

plot.lm.spike <- function(x,
                          burn = 0,
                          inclusion.threshold = 0,
                          unit.scale = TRUE,
                          number.of.variables = NULL,
                          ...) {
  ## S3 plotting method for lm.spike objects.  Produces a barplot of the
  ## marginal inclusion probabilities of each variable, sorted by the
  ## marginal inclusion probability, and shaded by the conditional
  ## probability that a coefficient is positive, given that it is
  ## nonzero.
  ##
  ## Args:
  ##   x: An object of class 'lm.spike' produced by the 'lm.spike'
  ##     function.
  ##   burn: Number of MCMC iterations to discard.  If burn <= 0 then
  ##     no iterations are discarded.
  ##   inclusion.threshold: Only plot coefficients with posterior
  ##     inclusion probabilites exceeding this value.
  ##   unit.scale: Logical indicating if the axis of the plot should
  ##     be forced to [0, 1].
  ##   number.of.variables: If non-NULL this specifies the number of
  ##     coefficients to plot, taking precedence over
  ##     inclusion.threshold.
  ##   ...:  Additional arguments to be passed to barplot.
  ## Returns:
  ## A list with the following elements:
  ##   barplot: The midpoints of each bar, which is useful for adding
  ##     to the plot.
  ##   inclusion.prob: The marginal inclusion probabilities of each
  ##     variable, ordered smallest to largest (the same ordering as
  ##     the plot).
  ##   positive.prob: The probability that each variable has a
  ##     positive coefficient, in the same order as inclusion.prob.
  ##   permutation: The permutation of beta that puts the coefficients
  ##     in the same order as positive.prob and inclusion.prob.  That
  ##     is: beta[, permutation] will have the most significant
  ##     coefficients in the right hand columns.

  beta <- x$beta
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }
  inclusion.prob <- colMeans(beta!=0)
  index <- order(inclusion.prob)
  beta <- beta[, index, drop = FALSE]
  inclusion.prob <- inclusion.prob[index]

  compute.positive.prob <- function(x) {
    ## Compute the probability that x is positive, given that it is
    ## nonzero.  If all(x == 0) then zero is returned.
    x <- x[x != 0]
    if (length(x) == 0) {
      return(0)
    }
    return(mean(x > 0))
  }
  positive.prob <- apply(beta, 2, compute.positive.prob)

  ## If there are too many variables floating around you won't be able
  ## to read the plot.  You can make a nicer looking plot by only
  ## showing the variables with inclusion probabilities above a
  ## threshold, or by specifying number.of.variables.
  if (!is.null(number.of.variables)) {
    stopifnot(number.of.variables > 0)
    show <- tail(seq(along = inclusion.prob), number.of.variables)
  } else {
    show <- inclusion.prob >= inclusion.threshold
  }

  bar.plot <- NULL
  if (sum(show) > 0) {
    omai <- par("mai")
    variable.names <- dimnames(beta)[[2]]
    omai[2] <- max(strwidth(variable.names[show], units = "inches")) + .5
    oldpar <- par(mai = omai)
    on.exit(par(oldpar))
    xlim <- c(0, max(inclusion.prob[show]))
    if (unit.scale) {
      xlim <- c(0, 1)
    }
    bar.plot <- barplot(inclusion.prob[show],
                        horiz = TRUE,
                        xlim = xlim,
                        col = gray(positive.prob[show]),
                        names.arg = variable.names[show],
                        las = 1,
                        xlab = "Inclusion Probability",
                        ...
                        )
  } else {
    bar.plot <-
      plot(1, 1,
           axes = F, type = "n",
           main = paste("No variables exceeded the inclusion ",
             "probability threshold of", inclusion.threshold))
  }
  ans <- list(barplot = bar.plot,
              inclusion.prob = inclusion.prob,
              positive.prob = positive.prob,
              permutation = index)
  return(invisible(ans))
}

PlotModelSize <- function(model, burn = 0,
                            xlab = "Number of nonzero coefficients", ...) {
  ## Plot a histogram of the number of nonzero coefficient in the model.
  ## Args:
  ##   model: A model object produced by lm.spike.
  ##   burn:  The number of MCMC iterations to discard as burn-in.
  ##   xlab:  Label for the x axis.
  ##   ...:   Extra arguments to be passed to 'hist'.
  ## Returns:
  ##   A vector giving the number of nonzero coefficients in each MCMC
  ##     iteration.

  stopifnot(inherits(model, "lm.spike"))
  beta <- model$beta
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }
  size <- rowSums(beta!=0)
  hist(size, xlab = xlab, ...)
  return(invisible(size))
}

SummarizeSpikeSlabCoefficients <- function(beta, burn = 0, order = TRUE) {
  ## Compute posterior means, posterior standard deviations, and
  ## posterior inclusion probabilities for the coefficients in the
  ## lm.spike object.
  ## Args:
  ##   beta: a matrix containing MCMC draws of regression coefficients.
  ##     Each row is a draw.  Each column is a coefficient.
  ##   burn:  The number of MCMC draws to be discarded as burn-in.
  ##   order: Logical.  If TRUE the output will be in decending order
  ##     according to the probability that each coefficient is
  ##     nonzero.
  ## Returns:
  ## A five-column matrix giving the posterior mean and standard
  ## deviation of each coefficient, the conditional mean and standard
  ## deviation of each coefficient given that it is nonzero, and the
  ## probability that the coefficient is nonzero.
  beta <- as.matrix(beta)
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
  }
  inclusion.prob <- colMeans(beta != 0)
  if (order) {
    index <- rev(order(inclusion.prob))
    beta <- beta[, index , drop = FALSE]
    inclusion.prob <- inclusion.prob[index]
  }

  ConditionalMean <- function(b) {
    ## Returns the mean of the nonzero elements of b
    b <- b[b != 0]
    if (length(b) > 0) return(mean(b))
    return(0)
  }
  ConditionalSd <- function(b) {
    ## Returns the standard deviation of the nonzero elements of b
    b <- b[b != 0]
    if (length(b) > 1) return(sd(b))
    return(0)
  }

  ## Short names are used here because they are displayed as column
  ## headers.  Names that are wider than the corresponding column of
  ## numbers make the result poorly spaced.
  coef <- cbind(mean = colMeans(beta),
                sd = apply(beta, 2, sd),
                mean.inc = apply(beta, 2, ConditionalMean),
                sd.inc = apply(beta, 2, ConditionalSd),
                inc.prob = inclusion.prob)

  return(coef)
}

summary.lm.spike <- function(object, burn = 0, order = TRUE, ...) {
  ## Summarize the coefficients and residual standard deviation from
  ## an lm.spike object.
  ## Args:
  ##   object: An object of class 'lm.spike', or an equivalent
  ##     list-like object containing a matrix of coefficient draws named
  ##     'beta'.
  ##   burn: An integer giving the number of draws to be discarded as
  ##     burn-in.
  ##   order: Logical.  If TRUE then the coefficients are presented in
  ##     order of their posterior inclusion probabilities.  Otherwise
  ##     the coefficients retain the order they had in 'object'.
  ## Returns:
  ## A list with three elements:
  ##   coefficients: A summary of the model coefficients produced by
  ##     SummarizeSpikeSlabCoefficients
  ##   residual.sd: A summary of the posterior distribution of the
  ##     residual standard deviation parameter ("sigma")
  ##   rsquare: A summary of the posterior distribution of the R^2
  ##     statistic: 1 - resdual.sd^2 / sample.sd^2, where sample.sd is
  ##     the sample standard deviation of the response variable.
  sigma <- object$sigma
  if (burn > 0) {
    sigma <- sigma[-(1:burn)]
  }
  ans <- list(coefficients =
              SummarizeSpikeSlabCoefficients(object$beta, burn, order),
              residual.sd = summary(sigma),
              rsquare = summary(1 - sigma^2 / object$sample.sd^2))
  class(ans) <- "summary.lm.spike"
  return(ans)
}

print.summary.lm.spike <- function(x, ...) {
  ## Print method for summary.lm.spike objects.  Only print 3
  ## significant digits.
  cat("coefficients:\n")
  print.default(signif(x$coefficients, 3), ...)
  cat("\nresidual.sd = \n")
  print(x$residual.sd)
  cat("\nr-square    = \n")
  print(x$rsquare)
  cat("\n")
}

predict.lm.spike <- function(object, newdata, burn = 0,
                             na.action = na.pass, ...) {
  ## Prediction method for lm.spike
  ## Args:
  ##   object: object of class "lm.spike" returned from the lm.spike
  ##     function
  ##   newdata: A data frame, matrix, or vector containing the
  ##     predictors needed to make the prediction.  If 'newdata' is a
  ##     data.frame it must contain variables with the same names as
  ##     the data frame used to fit 'object'.  If it is a matrix, it
  ##     must have the same number of columns as object$beta.  (An
  ##     intercept term will be implicitly added if the number of
  ##     columns is one too small.)  If the dimension of object$beta
  ##     is 1 or 2, then
  ##   burn: The number of MCMC iterations in 'object' that should be
  ##     discarded.  If burn < 0 then all iterations are kept.
  ##   na.action: what to do about NA's.
  ##   ...:  unused, but present for compatibility with predict().
  ## Returns:
  ##   A matrix of predictions, with each row corresponding to a row
  ##   in newdata, and each column to an MCMC iteration.
  beta.dimension <- ncol(object$beta)
  if (is.data.frame(newdata)) {
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

    if (nrow(X) != nrow(newdata)) {
      warning("Some entries in newdata have missing values, and will",
                   "be omitted from the prediction.")
    }
  } else if (is.matrix(newdata)) {
    X <- newdata
    if (ncol(X) == beta.dimension - 1) {
      if (attributes(object$terms)$intercept) {
        X <- cbind(1, X)
        warning("Implicit intercept added to newdata.")
      }
    }
  } else if (is.vector(newdata) & beta.dimension == 2) {
    if (attributes(object$terms)$intercept) {
      X <- cbind(1, newdata)
    }
  } else if (is.vector(newdata & beta.dimension == 1)) {
    X <- matrix(newdata, ncol=1)
  } else {
    stop("Error in predict.lm.spike:  newdata must be a matrix,",
         "or a data.frame,",
         "unless dim(beta) <= 2, in which case it can be a vector")
  }

  if (ncol(X) != beta.dimension) {
    stop("The number of coefficients does not match the number",
         "of predictors in lm.spike")
  }

  beta <- object$beta
  sigma <- object$sigma
  if (burn > 0) {
    beta <- beta[-(1:burn), , drop = FALSE]
    sigma <- sigma[-(1:burn)]
  }

  eta <- X %*% t(beta)

  ## eta is the [n x niter] matrix of predicted values
  y.new <- matrix(rnorm(length(eta),
                        eta,
                        rep(sigma, each = nrow(X))),
                  nrow = nrow(X))
  return(y.new)
}

.lm.spike.fit <- function(x, y, prior, niter, ping, seed) {
  ## Driver function for lm.spike.
  ## Args:
  ##   x: design matrix with 'n' rows corresponding to observations and
  ##     'p' columns corresponding to predictor variables.
  ##   y: vector of responses of length n
  ##   prior: a list structured like the return value from SpikeSlabPrior.
  ##   niter:  the number of desired MCMC iterations
  ## Returns:
  ##   An object of class lm.spike, which is a list with the elements
  ##   described in the 'lm.spike' function.
  stopifnot(nrow(x) == length(y))
  stopifnot(length(prior$mu) == ncol(x))
  stopifnot(length(prior$prior.inclusion.probabilities) == ncol(x))

  ## The following is for backward compatibility.  A max.flips
  ## argument was added to SpikeSlabPrior.
  if (is.null(prior$max.flips)) {
    prior$max.flips <- -1
  }

  ans <- .Call(do_spike_slab,
               x,
               y,
               prior,
               niter,
               ping,
               seed)
  variable.names <- dimnames(x)[[2]]
  if (!is.null(variable.names)) {
    colnames(ans$beta) <- variable.names
  }
  ans$prior <- prior
  class(ans) <- "lm.spike"
  return(ans)
}
