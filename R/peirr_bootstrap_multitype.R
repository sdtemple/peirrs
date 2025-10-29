#' Bootstrap pair-based estimators for multitype epidemic model
#'
#' Perform parametric bootstrapping to estimate infection and removal rates
#' in a multitype stochastic epidemic model. Each bootstrap replicate simulates
#' a new epidemic under the given parameters, re-estimates rates using
#' \code{peirr_tau_multitype}, and stores the results.
#'
#' @param num.bootstrap Integer: number of bootstrap replicates.
#' @param betas Numeric vector: infection rates for each infection type.
#' @param gammas Numeric vector: removal rates for each removal type.
#' @param beta.sizes Integer vector: population sizes or group sizes corresponding to \code{betas}.
#' @param gamma.sizes Integer vector: population sizes or group sizes corresponding to \code{gammas}.
#' @param sample.size Numeric: target epidemic size to condition on.
#' @param p Numeric: expected proportion of complete pairs observed.
#' @param q Numeric: probability that infection time is missing (conditional on missingness).
#' @param m Integer: positive shape parameter for infectious period distribution (default 1).
#' @param e Numeric: fixed exposure period (default 0).
#' @param within Numeric: acceptable proportional deviation from the target sample size (default 0.1).
#' @param etc List or NULL: additional arguments passed to \code{peirr_tau_multitype}.
#'
#' @return A numeric matrix of dimension \code{(num.bootstrap + 1) × (length(betas) + length(gammas))}.
#' The first row contains the original \code{betas} and \code{gammas}, and subsequent rows
#' contain bootstrap estimates of infection and removal rates.
#'
#' @details
#' Each bootstrap replicate:
#' \enumerate{
#'   \item Simulates a multitype epidemic using \code{simulator_multitype}.
#'   \item Filters and extracts infection and removal times by type.
#'   \item Re-estimates rates with \code{peirr_tau_multitype}.
#' }
#' The resulting matrix can be used to assess estimation uncertainty.
#'
#' @examples
#' \dontrun{
#' results <- peirr_bootstrap_multitype(
#'   num.bootstrap = 100,
#'   betas = c(2,2),
#'   gammas = c(1,2),
#'   beta.sizes = c(110,90)
#'   gamma.sizes = c(90,110),
#'   sample.size = 100,
#'   p = 0.5,
#'   q = 0.5
#' )
#' }
#'
#' @export
peirr_bootstrap_multitype <- function(num.bootstrap,
                            betas,
                            gammas,
                            beta.sizes,
                            gamma.sizes,
                            sample.size,
                            p,
                            q,
                            m=1,
                            e=0,
                            within=0.1,
                            etc=NULL
){

  # for matrix setup and filling
  num.beta = length(betas)
  num.gamma = length(gammas)
  num.col = num.beta + num.gamma
  
  # initialize matrix
  storage <- matrix(0,nrow=num.bootstrap+1,ncol=num.col)
  storage[1,1:num.beta] = betas
  storage[1,(num.beta+1):num.col] = gammas
  
  # condition on epidemic being of a certain size
  min.sample.size = (1 - within) * sample.size
  max.sample.size = (1 + within) * sample.size
  
  # bootstrapping
  for(l in 2:(num.bootstrap+1)){
    out = simulator_multitype(betas,
                              gammas,
                              beta.sizes,
                              gamma.sizes,
                              m,
                              e,
                              p,
                              q,
                              min.sample.size,
                              max.sample.size)
    X = out$matrix.time
    r = X[,2]
    cr = X[,5]
    i = X[,1]
    ci = X[,3]
    bth.estimate = do.call(peirr_tau_multitype, 
                           c(list(r=r,i=i,cr=cr,ci=ci,Ns=beta.sizes), etc)
                           )
    # print(bth.estimate)
    storage[l,1:num.beta] = bth.estimate$infection.rates
    storage[l,(num.beta+1):num.col] = bth.estimate$removal.rates
  }
  return(storage)
}