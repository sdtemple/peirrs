#' Pair-based tau estimator of common infection and removal rates with spatial effect
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param lag numeric: fixed exposure period
#' @param tau.med bool: use median imputation for tau if TRUE
#' @param gamma.med bool: TRUE for median, and FALSE for mean in estimating the removal rate
#' @param h function: symmetric function of distance
#' @param D numeric: two-dimensional distance matrix
#'
#' @return numeric list (infection.rate, removal.rate)
#'
#' @export
peirr_tau_spatial <- function(r, i, N, h, D,
                      lag = 0,
                      tau.med = FALSE,
                      gamma.med = FALSE
                      ) {

  # make sure one of the other is finite
  or.finite <- is.finite(r) | is.finite(i)
  i <- i[or.finite]
  r <- r[or.finite]

  # estimate of removal rate
  gamma.estim <- peirr_removal_rate(r, i, gamma.med)

  # number of infected
  n <- sum(!is.na(r) | !is.na(i))
  if (sum(is.na(r)) >= n) {
    stop("There are no complete case periods to estimate the removal rate")
  }
  if (sum(is.na(i)) >= n) {
    stop("There are no complete case periods to estimate the removal rate")
  }

  ### first infected ###
  ralpha <- which.min(r)
  ialpha <- which.min(i)
  if (i[ialpha] < r[ralpha]) {
    alpha <- ialpha
  } else {
    alpha <- ralpha
  }

  # sum up tau terms
  tau <- 0
  for (j in (1:n) [-alpha]) {
    rj <- r[j]
    ij <- i[j]
    for (k in (1:n)[-j]) {
      rk <- r[k]
      ik <- i[k]
      tm <- tau_moment(rk, rj, ik, ij, gamma.estim, gamma.estim, lag, tau.med) * h(D[k,j])
      if (is.na(tm)) {
        print(c(rk, rj, ik, ij, gamma.estim, gamma.estim))
        }
      tau <- tau + tm
    }
  }

  # these are where we have full infectious period
  full.r <- r[(!is.na(r)) & (!is.na(i))]
  full.i <- i[(!is.na(r)) & (!is.na(i))]

  # have to address this component with spatial function
  # only take expectation when we don't have the full period
  num.not.full <- length(r) - length(full.r)
  median.scalar <- 1
  if (tau.med) {median.scalar <- log(2)}
  if (n == N) {
    ri.sum <- 0
  } else{
    ri.sum <- 0
    for (j in 1:n){
      rjij <- r[j] - i[j]
      if (is.na(rjij)){ rjij <- 1 / gamma.estim * median.scalar }
      for (k in (n+1):N) {
        ri.sum <- ri.sum + rjij * h(D[j,k])
      }
    }
  }


  # maximizes give conditional expectations
  beta.estim <- (n - 1) / (tau + ri.sum)

  return(list(infection.rate = beta.estim * N,
              removal.rate = gamma.estim,
              tau.sum = tau,
              full.ri.sum = ri.sum,
              num.not.infected = N-n
              ))
}
