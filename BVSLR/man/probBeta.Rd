% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BVSLR.R
\name{probBeta}
\alias{probBeta}
\title{probBeta}
\usage{
probBeta(betat1, At1, gt1, sig2t1, data.X, data.Y, n.size, p.size, log.p = F)
}
\arguments{
\item{betat1}{\eqn{\beta} at time t + 1}

\item{At1}{random set, A at time t + 1}

\item{gt1}{g-prior, g at time t + 1}

\item{sig2t1}{\eqn{\sigma^2} at time t + 1}

\item{data.X}{predictors}

\item{data.Y}{responses}

\item{n.size}{sample size}

\item{p.size}{number of predictor}

\item{log.p}{logical. If TRUE, probability p is given as log(p)}
}
\value{
proposal density of \eqn{\beta}
}
\description{
probBeta
}
