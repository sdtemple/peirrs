#' Pair-based tau estimator of common infection and removal rates
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param med bool: use median imputation if true
#'
#' @return numeric list (infection.rate, removal.rate, R0)
#'
#' @export
peirr_tau <- function(r, i, N, med=TRUE){

  # make sure one of the other is finite
  or.finite <- is.finite(r)|is.finite(i)
  i <- i[or.finite]
  r <- r[or.finite]

  # estimate of removal rate
  gamma.estim <- mle_removal_rate(r,i)

  # number of infected
  n <- sum(!is.na(r) | !is.na(i))
  if(sum(is.na(r))>=n){
    stop("There are no complete case periods to estimate the removal rate")
  }
  if(sum(is.na(i))>=n){
    stop("There are no complete case periods to estimate the removal rate")
  }

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
      tm <- tau_moment(rk, rj, ik, ij, gamma.estim, gamma.estim, med)
      if(is.na(tm)){print(c(rk,rj,ik,ij))}
      tau <- tau + tm
    }
  }

  # these are where we have full infectious period
  full.r <- r[(!is.na(r)) & (!is.na(i))]
  full.i <- i[(!is.na(r)) & (!is.na(i))]

  # only take expectation when we don't have the full peroid
  num.not.full <- length(r) - length(full.r)
  median.scalar <- 1
  if(med){median.scalar <- log(2)}
  ri.sum <- num.not.full / gamma.estim * median.scalar + sum(full.r - full.i)

  # maximizes give conditional expectations
  beta.estim <- (n - 1) / (tau + (N - n) * ri.sum)

  return(list(infection.rate=beta.estim*N,
              removal.rate=gamma.estim,
              R0=beta.estim*N/gamma.estim,
              tau.sum=tau,
              full.ri.sum=ri.sum,
              num.not.infected=N-n
              ))
}
