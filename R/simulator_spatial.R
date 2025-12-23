#' Simulate spatial stochastic epidemic model with formatting
#'
#' Draw infectious periods for spatial stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param N integer: population size
#' @param h function: symmetric function of distance
#' @param D numeric: two-dimensional distance matrix
#' @param m integer: positive shape
#' @param lag numeric: fixed exposure period
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' @param min.sample.size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time), matrix of n by N distances
#'
#' @export
simulator_spatial <- function(beta,
                      gamma,
                      N,
                      h,
                      D,
                      m = 1,
                      lag = 0,
                      p = 0.5,
                      q = 1,
                      min.sample.size = 10,
                      max.sample.size = Inf) {
  sample.size <- 0
  gamma.estim <- NA
  if (p <= 0) {
    stop("p <= 0 error. Must have some complete infectious periods.")
  }
  while ((sample.size <= min.sample.size) || (sample.size >= max.sample.size) || is.na(gamma.estim)) {
    # main simulation
    epi <- simulate_sem_spatial(beta, gamma, N, h, D, m, lag)
    filter_indices <- is.finite(epi$matrix.time[,1]) & is.finite(epi$matrix.time[,2])
    epi$matrix.time <- filter_sem(epi$matrix.time)
    epi$matrix.distance <- epi$matrix.distance[filter_indices, ]
    epi$matrix.time <- decomplete_sem(epi$matrix.time, p = p, q = q)
    sort_indices <- order(epi$matrix.time[,2]) # sort the distance matrix
    epi$matrix.time <- sort_sem(epi$matrix.time)
    # calculate the sample size
    sample.size <- dim(epi$matrix.time)[1]
    if (sample.size > 1){
      # sorting runs into error when sample.size <= 1
      epi$matrix.distance <- epi$matrix.distance[sort_indices, ]
    }

    # ensure there are r and i enough to estimate gamma
    X <- epi$matrix.time
    r <- X[,2]
    i <- X[,1]
    gamma.estim <- peirr_removal_rate(r, i)
    }
    return(epi)
  }
