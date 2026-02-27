#' Simulate stochastic epidemic model with spatial distance effect
#'
#' Draw infectious periods for SEIR model with spatial effect.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param population_size integer: population size
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed incubation period
#' @param kernel_spatial function: symmetric function of distance
#' @param matrix_distance numeric: two-dimensional distance matrix
#'
#' @details
#' This function implements a spatial stochastic SIR epidemic model using an
#' event-driven (Gillespie) algorithm. The key difference from non-spatial models
#' is that the infection rate between individuals depends on their spatial distance
#' via the `kernel_spatial` function.
#'
#' The infection rate at time t is computed as:
#' beta/N * sum over (infectious, susceptible) pairs of kernel_spatial(distance).
#'
#' The `kernel_spatial` function should be a symmetric, non-negative function of
#' distance (e.g., exponential decay, power law, or step function). Common choices
#' include `function(d) exp(-lambda * d)` or `function(d) 1 / (1 + d^2)`.
#'
#' The `matrix_distance` is an N x N symmetric distance matrix where element (i, j)
#' gives the distance between individuals i and j. Typically constructed via
#' \code{as.matrix(dist(coordinates))} for coordinates in Euclidean space.
#'
#' @return A list with three elements:
#' \itemize{
#'   \item `matrix_time`: an N x 2 matrix with columns "infection" and "removal" containing
#'     exposure and removal times for each individual (Inf indicates never infected or still infectious)
#'   \item `matrix_record`: a T x 5 matrix with time-indexed columns St, Et, It, Rt, Time
#'     recording the susceptible, exposed, infectious, and removed counts plus elapsed time
#'   \item `matrix_distance`: the input distance matrix, returned for reference
#' }
#'
#' @examples
#' # Construct 2D spatial coordinates and distance matrix
#' set.seed(1)
#' n <- 50
#' coords <- cbind(runif(n), runif(n))
#' D <- as.matrix(dist(coords))
#'
#' # Define an exponential decay kernel
#' kernel_spatial <- function(d) exp(-2 * d)
#'
#' # Simulate spatial epidemic
#' epi1 <- simulate_sem_spatial(beta = 1.5, gamma = 1.0, population_size = n,
#'                              kernel_spatial = kernel_spatial,
#'                              matrix_distance = D)
#' head(epi1$matrix_time)
#'
#' # Simulation with lag and multiple renewals
#' set.seed(2)
#' kernel_spatial2 <- function(d) 1 / (1 + d^2)  # Power-law kernel
#' epi2 <- simulate_sem_spatial(beta = 2.0, gamma = 0.8, population_size = n,
#'                              kernel_spatial = kernel_spatial2,
#'                              matrix_distance = D,
#'                              lag = 0.5, num_renewals = 2)
#' hist(epi2$matrix_time[, "removal"] - epi2$matrix_time[, "infection"],
#'      main = "Spatial epidemic: infectious period distribution")
#'
#' @keywords internal
simulate_sem_spatial <- function(beta,
                                  gamma,
                                  population_size,
                                  kernel_spatial,
                                  matrix_distance,
                                  num_renewals = 1,
                                  lag = 0
                                  ) {

  # initialize vectors
  t <- 0
  betaN <- beta / population_size
  infections <- rep(Inf, population_size)
  removals <- rep(Inf, population_size)
  renewals <- rep(0, population_size)
  alpha <- sample(population_size, 1)
  infections[alpha] <- t

  # simulate epidemic
  St <- sum(is.infinite(infections))
  It <- sum(is.finite(infections)) - sum(is.finite(removals))
  Et <- 0
  Rt <- 0

  # recording the evolution
  susceptible_recording <- c(St)
  infectious_recording <- c(It)
  exposed_recording <- c(Et)
  removed_recording <- c(Rt)
  time_recording <- c(0)
  ctr <- 1

  while ( (It > 0) || (Et > 0) ) {

    # closest infectious time after exposure
    min_time <- min(
      infections[is.infinite(removals) & is.finite(infections) & (infections > t)],
      Inf
    )

    if (It == 0) {
      # no infecteds but there are exposeds
      # the closest exposure wait
      t <- min_time + .Machine$double.eps
    } else {
      # simulate time
      infection_rate <- betaN * sum(kernel_spatial(matrix_distance[is.infinite(removals) & is.finite(infections) & (infections <= t), is.infinite(infections)]))
      removal_rate <- gamma * It
      t <- t + rexp(1, rate = infection_rate + removal_rate)

      if (t > min_time) {
        # update time to make an exposed infectious
        t <- min_time + .Machine$double.eps
      } else {
        # there is infection or removal before
        # simulate transition
        x <- rbinom(1, size = 1, prob = removal_rate / (infection_rate + removal_rate))
        x <- (x + 1) %% 2
        if (x) {
          # infect a susceptible
          if(St > 1){
            argx <- sample(which(is.infinite(infections) & is.infinite(removals), arr.ind = TRUE), 1)
          } else{
            argx <- which(is.infinite(infections) & is.infinite(removals))
          }
          infections[argx] <- t + lag # fixed exposure period
        } else {

          # remove an infected
          if (It > 1) {
            argx <- sample(which(is.infinite(removals) & (infections <= t), arr.ind = TRUE), size=1)
            # & (i <= t) means can't be removed before infectious when exposed
          } else {
            # when epidemic is winding down
            # no more infecteds
            argx <- which(is.infinite(removals) & (infections <= t))
            # & (i <= t) means can't be removed before infectious when exposed
          }
          renewals[argx] <- renewals[argx] + 1
          if (renewals[argx] == num_renewals) { # after m renewals
            removals[argx] <- t
          }
        }
      }
    }

    # update (S,I) counts
    St = sum(is.infinite(infections))
    It = sum(is.finite(infections) & (infections <= t)) - sum(is.finite(removals))
    Rt = sum(is.finite(infections) & is.finite(removals))
    Et = sum(is.finite(infections) & (infections > t))
    if(St + Rt + Et + It != population_size){
      stop("S(t) + I(t) + E(t) + R(t) do not equal N")
    }
    # & (i <= t) delays the infectious period after exposure

    susceptible_recording <- c(susceptible_recording, St)
    infectious_recording <- c(infectious_recording, It)
    exposed_recording <- c(exposed_recording, Et)
    removed_recording <- c(removed_recording, Rt)
    time_recording <- c(time_recording, t)
    ctr <- ctr + 1
  }

  # there should be no negatives
  # and ignore the Inf values
  non_negative_check <- removals - infections
  non_negative_check <- non_negative_check[is.finite(non_negative_check)]
  if( any( (non_negative_check) < 0)){
    stop("At least one r - i value is negative")
  }

  # formatting
  output <- matrix(c(infections, removals),
                  nrow = population_size,
                  ncol = 2,
                  byrow = FALSE)
  colnames(output) <- c("infection", "removal")

  recording <- matrix(c(susceptible_recording, exposed_recording, infectious_recording, removed_recording, time_recording),
                     nrow = ctr,
                     ncol = 5,
                     byrow = FALSE
                     )
  colnames(recording) <- c("St", "Et", "It", "Rt", "Time")

  return(list(matrix_time = output,
              matrix_record = recording,
              matrix_distance = matrix_distance
              ))
}
