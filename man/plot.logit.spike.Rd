\name{plot.logit.spike}
\Rdversion{1.1}
\alias{plot.logit.spike}
\title{
  Plot a \code{\link{logit.spike}} object.
}
\description{


  The default plot is a barplot of the marginal inclusion probabilities
  for each variable, as obtained by
  \code{\link{PlotMarginalInclusionProbabilities}}.

}

\usage{
  \method{plot}{logit.spike}(x,
     y = c("coefficients", "fit", "residuals", "size"),
     ...)
}
\arguments{
  \item{x}{
    An object of class \code{logit.spike}.
  }

  \item{y}{
    The type of plot desired.
  }

  \item{\dots}{ Additional arguments passed to the specific functions
    that do the plotting.  }

}

\author{
  Steven L. Scott
}

\seealso{
  \code{\link{PlotMarginalInclusionProbabilities}}
  \code{\link{PlotModelSize}}
  \code{\link{PlotLogitSpikeFitSummary}}
  \code{\link{PlotLogitSpikeResiduals}}
}

\examples{
## See the examples in ?logit.spike
}
\keyword{models}
\keyword{regression}
