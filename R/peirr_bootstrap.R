#' Bootstrap pair-based estimators for epidemic parameters
#'
#' Perform a parametric bootstrap procedure to estimate variability in
#' infection and removal rate estimates of an epidemic model.
#'
#' This function repeatedly simulates epidemic data under the specified
#' model parameters, applies a user-specified estimator function (typically
#' \code{peirr}), and stores the resulting infection and removal rate
#' estimates for each bootstrap iteration.
#'
#' @param num.bootstrap Integer. Number of bootstrap replicates to perform.
#' @param beta Numeric. Infection rate parameter used for simulation.
#' @param gamma Numeric. Removal rate parameter used for simulation.
#' @param N Integer. Population size.
#' @param sample.size Integer. Number of individuals sampled in each
#'   bootstrap replicate.
#' @param p Numeric. Expected proportion of complete pairs observed.
#' @param q Numeric. Probability that an infection time is missing.
#' @param peirr Function (default \code{peirr_tau}). A function that estimates infection and removal
#'   rates given infection and removal time data (e.g., \code{peirr()}).
#' @param m Integer (default = 1). Positive shape parameter for the infection period.
#' @param e Numeric (default = 0). Fixed exposure period duration.
#' @param within Numeric (default = 0.1). Fractional range around
#'   \code{sample.size} used to generate random sample sizes in each
#'   bootstrap replicate.
#' @param etc List (default = NULL). Additional arguments passed to
#'   \code{peirr()} via \code{do.call()}.
#'
#' @details
#' For each bootstrap iteration, this function:
#' \enumerate{
#'   \item Simulates an epidemic using \code{simulator(beta, gamma, N, m, e, p, q, min.sample.size, max.sample.size)}.
#'   \item Extracts infection and removal times from the simulated data.
#'   \item Applies \code{peirr()} (or a user-supplied estimator) to obtain
#'         infection and removal rate estimates.
#'   \item Stores the resulting parameter estimates in a matrix.
#' }
#'
#' The first row of the returned matrix stores the true values of
#' \code{beta} and \code{gamma}.
#'
#' @return
#' A numeric matrix with \code{num.bootstrap + 1} rows and 2 columns:
#' \itemize{
#'   \item Column 1: Estimated infection rates (\code{beta}).
#'   \item Column 2: Estimated removal rates (\code{gamma}).
#' }
#' The first row contains the true values used for simulation.
#'
#' @examples
#' \dontrun{
#' results <- peirr_bootstrap(
#'   num.bootstrap = 100,
#'   beta = 2,
#'   gamma = 1,
#'   peirr = peirr_tau,
#'   N = 500,
#'   sample.size = 100,
#'   p = 0.5,
#'   q = 0.5
#' )
#' }
#'
#' @export
peirr_bootstrap <- function(num.bootstrap,
                            beta,
                            gamma,
                            N,
                            sample.size,
                            p,
                            q,
                            peirr=peirr_tau,
                            m=1,
                            e=0,
                            within=0.1,
                            etc=NULL
                            ){

  # override
  if(identical(body(peirr),body(peirr_pbla_infection_rate))){
    q=1
  }
  if(identical(body(peirr),body(peirr_pbla_both_rates))){
    p=0
    q=1
  }

  storage <- matrix(0,nrow=num.bootstrap+1,ncol=2)
  storage[1,1] = beta
  storage[1,2] = gamma
  min.sample.size = (1 - within) * sample.size
  max.sample.size = (1 + within) * sample.size
  for(l in 2:(num.bootstrap+1)){
    out = simulator(beta,gamma,N,m,e,p,q,min.sample.size,max.sample.size)
    X = out$matrix.time
    r = X[,2]
    i = X[,1]
    bth.estimate = do.call(peirr, c(list(r=r,i=i,N=N), etc))
    storage[l,1] = bth.estimate$infection.rate
    storage[l,2] = bth.estimate$removal.rate
  }
  return(storage)
}
