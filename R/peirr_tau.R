#' Pair-based tau estimator of common infection and removal rates
#' 
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging. 
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' 
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'  
#' @export 
peirr_tau <- function(r, i, N){
  
  # estimate of removal rate
  gamma.estim <- mle_removal_rate(r,i)
  
  # number of infected
  n <- sum(!is.na(r) | !is.na(i))
  
  ### first infected ###
  ralpha <- which.min(r)
  ialpha <- which.min(i)
  if(i[ialpha] < r[ralpha]){
    alpha <- ialpha
  } else{
    alpha <- ralpha
  }
  
  # sum up tau terms
  tau <- 0
  for(j in (1:n)[-alpha]){
    rj <- r[j]
    ij <- i[j]
    for(k in (1:n)[-j]){
      rk <- r[k]
      ik <- i[k]
      tau <- tau + tau_moment(rk, rj, ik, ij, gamma.estim, gamma.estim)
    }
  }

  # these are where we have full infectious period
  full.r <- r[(!is.na(r)) & (!is.na(i))]
  full.i <- i[(!is.na(r)) & (!is.na(i))]
  
  # maximizes give conditional expectations
  beta.estim <- (n - 1) / (tau + (N - n) / gamma.estim * n)
  
  return(list(infection.rate=beta.estim*N, 
              removal.rate=gamma.estim, 
              R0=beta.estim*N/gamma.estim, 
              tau.sum=tau))
}
