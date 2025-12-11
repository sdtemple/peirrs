#' Simulate general stochastic epidemic model with formatting
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param N integer: population size
#' @param m integer: positive shape
#' @param lag numeric: fixed exposure period
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' @param min.sample.size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulator <- function(beta,
                      gamma,
                      N,
                      m = 1,
                      lag = 0,
                      p = 0.5,
                      q = 1,
                      min.sample.size = 10,
                      max.sample.size = Inf) {
  sample.size <- 0
  gamma.estim <- NA
  while ((sample.size <= min.sample.size) || (sample.size >= max.sample.size) || is.na(gamma.estim)) {
    # main simulation
    epi <- simulate_sem(beta, gamma, N, m, lag)
    epi$matrix.time <- filter_sem(epi$matrix.time)
    epi$matrix.time <- decomplete_sem(epi$matrix.time, p = p, q = q)
    epi$matrix.time <- sort_sem(epi$matrix.time)
    # calculate the sample size
    sample.size <- dim(epi$matrix.time)[1]
    # ensure there are r and i enough to estimate gamma
    X <- epi$matrix.time
    r <- X[,2]
    i <- X[,1]
    gamma.estim <- peirr_removal_rate(r, i)
  }
  return(epi)
}
