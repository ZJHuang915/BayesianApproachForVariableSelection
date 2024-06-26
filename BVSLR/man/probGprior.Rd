% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BVSLR.R
\name{probGprior}
\alias{probGprior}
\title{probGprior}
\usage{
probGprior(gt1, gt, eps_g, log.p = F)
}
\arguments{
\item{gt1}{g-prior, g at time t + 1}

\item{gt}{g-prior, g at time t}

\item{eps_g}{parameter for proposal function of g}

\item{log.p}{logical. If TRUE, probability p is given as log(p)}
}
\value{
proposal density of g
}
\description{
probGprior
}
