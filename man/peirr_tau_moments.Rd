% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peirr_tau_moments.R
\name{peirr_tau_moments}
\alias{peirr_tau_moments}
\title{Pair-based tau estimator of common infection and removal rates}
\usage{
peirr_tau_moments(r, i, N)
}
\arguments{
\item{r}{numeric vector: removal times}

\item{i}{numeric vector: infection times}

\item{N}{integer: population size}
}
\value{
numeric list (infection.rate, removal.rate, R0, tau.sum)
}
\description{
Estimate infection and removal rates with tau-based expectation-maximization.
}
