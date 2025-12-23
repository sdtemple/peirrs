#' Posterior parameters for infection and removal rates given complete data
#'
#' Parameters for independent Gibbs sampling of the infection and removal rates.
#' You can sample from \code{rgamma()} with the posterior parameters to get the posterior distribution.
#' For the infection rate, you should scale the result of \code{rgamma()} by the population size.
#' \code{infection.rate.samples} and \code{removal.rate.samples} contain the posterior samples.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param infection.rate.rate.prior numeric
#' @param infection.rate.shape.prior numeric
#' @param removal.rate.rate.prior numeric
#' @param removal.rate.shape.prior numeric
#' @param num.posterior.samples numeric
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list of exact posterior and prior parameters
#'
#' @export
bayes_complete_data <- function(r,i,N,
                                infection.rate.rate.prior = 1e-5,
                                infection.rate.shape.prior = 1e-3,
                                removal.rate.rate.prior = 1e-3,
                                removal.rate.shape.prior = 1e-3,
                                num.posterior.samples = 1e4,
                                lag=0
                                ){
  n = length(r)
  tau = matrix(0, nrow = n, ncol = N)
  for(j in 1:n){
    tau[j,1:n] = sapply(i - lag, min, r[j]) - sapply(i - lag, min, i[j])
  }
  tau[,(n+1):N] = r - i
  Ai = sum(tau)
  Ci = sum(r-i)
  beta.samples = rgamma(num.posterior.samples,
                        rate = infection.rate.rate.prior+Ai,
                        shape = infection.rate.shape.prior+n-1
                        )
  gamma.samples = rgamma(num.posterior.samples,
                         rate = removal.rate.rate.prior+Ci,
                         shape = removal.rate.shape.prior+n
                         )
  return(list(infection.rate.samples=beta.samples * N,
              removal.rate.samples=gamma.samples,
              infection.rate.rate=infection.rate.rate.prior+Ai,
              infection.rate.shape=infection.rate.shape.prior+n-1,
              infection.rate.rate.prior=infection.rate.rate.prior,
              infection.rate.shape.prior=infection.rate.shape.prior,
              removal.rate.rate=removal.rate.rate.prior+Ci,
              removal.rate.shape=removal.rate.shape.prior+n,
              removal.rate.rate.prior=removal.rate.rate.prior,
              removal.rate.shape.prior=removal.rate.shape.prior
  ))
}
