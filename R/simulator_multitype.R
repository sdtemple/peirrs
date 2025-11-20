#' Simulate multi-type stochastic epidemic model with formatting
#'
#' Draw infectious periods for stochastic epidemic with different classes.
#'
#' @param betas numeric vector: infection rates
#' @param gamma numeric: removal rates
#' @param beta.sizes integers: subpopulation sizes for infection rates
#' @param gamma.sizes integers: subpopulation sizes for removal rates
#' @param m integer: positive shape
#' @param lag numeric: fixed exposure period
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' @param min.sample.size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulator_multitype <- function(betas,
                      gammas,
                      beta.sizes,
                      gamma.sizes,
                      m = 1,
                      lag = 0,
                      p = 0,
                      q = 1,
                      min.sample.size = 10,
                      max.sample.size = Inf
                      ){
  sample.size <- 0
  gamma.estim <- NA
  while((sample.size <= min.sample.size) || (sample.size >= max.sample.size) || is.na(gamma.estim) ){
    # main simulation
    epi = simulate_sem_multitype(betas, gammas, beta.sizes, gamma.sizes, m, lag)
    epi$matrix.time = filter_sem(epi$matrix.time)
    epi$matrix.time = decomplete_sem(epi$matrix.time, p=p, q=q)
    epi$matrix.time = sort_sem(epi$matrix.time)
    # calculate the sample size
    sample.size = dim(epi$matrix.time)[1]
    # ensure there are r and i enough to estimate all gammas
    gamma.estim <- 0
    for(n in 1:length(gamma.sizes)){
      X <- epi$matrix.time
      r <- X[, 2]
      i <- X[, 1]
      cr <- X[, 5]
      rg <- r[cr==n]
      ig <- i[cr==n]
      rig <- rg - ig
      rig <- rig[!is.na(rig)]
      gamma.estim <- gamma.estim + 1 / mean(rig)
    }
  }
  return(epi)
}
