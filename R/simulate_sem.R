#' Simulate general stochastic epidemic model
#'
#' Draw infectious periods for SEIR model.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param population_size integer: population size
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed incubation period
#'
#' @details
#' This function implements a stochastic SIR (Susceptible-Infected-Removed) epidemic
#' model using an event-driven (Gillespie) algorithm. The simulation proceeds by
#' sampling infection and removal events at random times according to their rates.
#'
#' The infection rate is `beta`, scaled by the number of susceptible-infectious pairs
#' (beta*S*I/N). The removal rate is `gamma` per infected individual.
#'
#' The `lag` parameter represents a fixed incubation period: when a susceptible is infected
#' at time t, they become infectious at time t + lag.
#'
#' The `num_renewals` parameter allows for multiple stages of infection before removal,
#' implementing a Negative Binomial infectious period distribution (mean 1/gamma per stage).
#'
#' The function returns both the infection/removal times for each individual and a
#' time-indexed recording of the S, E, I, R counts throughout the epidemic.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `matrix_time`: an N x 2 matrix with columns "infection" and "removal" containing
#'     exposure and removal times for each individual (Inf indicates never infected or still infectious)
#'   \item `matrix_record`: a T x 5 matrix with time-indexed columns St, Et, It, Rt, Time
#'     recording the susceptible, exposed, infectious, and removed counts plus elapsed time
#' }
#'
#' @examples
#' # Basic SIR epidemic simulation
#' set.seed(1)
#' epi1 <- simulate_sem(beta = 1.5, gamma = 1.0, population_size = 100)
#' head(epi1$matrix_time)
#' tail(epi1$matrix_record)
#'
#' # Simulation with exposure lag and multiple infection renewals
#' set.seed(2)
#' epi2 <- simulate_sem(beta = 2.0, gamma = 0.8, population_size = 80,
#'                      lag = 0.5, num_renewals = 2)
#' hist(epi2$matrix_time[, "removal"] - epi2$matrix_time[, "infection"],
#'      main = "Infectious period distribution")
#'
#' @keywords internal
simulate_sem <- function(beta,
                          gamma,
                          population_size,
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
  infection_recording <- c(It)
  exposed_recording <- c(Et)
  removal_recording <- c(Rt)
  time_recording <- c(0)
  ctr <- 1

  while ( (It > 0) || (Et > 0) ) {

    # closest infectious time after exposure
    min.time <- min(
      infections[is.infinite(removals) &
        is.finite(infections) &
        (infections > t)],
      Inf
    )

    if (It == 0) {
      # no infecteds but there are exposeds
      # the closest exposure wait
      t <- min.time + .Machine$double.eps

    } else {
      # simulate time
      irate <- betaN * It * St
      rrate <- gamma * It
      t <- t + rexp(1, rate = irate + rrate)

      if (t > min.time) {
        # update time to make an exposed infectious
        t <- min.time + .Machine$double.eps

      } else {
        # there is infection or removal before
        # simulate transition
        x <- rbinom(1, size = 1, prob = rrate / (irate + rrate))
        x <- (x + 1) %% 2
        if (x) {
          # infect a susceptible
          if (St > 1) {
            argx <- sample(which( is.infinite(infections) & is.infinite(removals), arr.ind = TRUE), 1)
          } else {
            argx <- which( is.infinite(infections) & is.infinite(removals) )
          }
          infections[argx] <- t + lag # fixed exposure period

        } else {

          # remove an infected
          if (It > 1) {
            argx <- sample(which( is.infinite(removals) & (infections <= t), arr.ind = TRUE), size=1)
            # & (i <= t) means can't be removed before infectious when exposed
          } else {
            # when epidemic is winding down
            # no more infecteds
            argx <- which( is.infinite(removals) & (infections <= t) )
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
    infection_recording <- c(infection_recording, It)
    exposed_recording <- c(exposed_recording, Et)
    removal_recording <- c(removal_recording, Rt)
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
                    byrow = FALSE
                    )
  colnames(output) <- c("infection", "removal")

  recording <- matrix(c(susceptible_recording,
                        exposed_recording,
                        infection_recording,
                        removal_recording,
                        time_recording),
                      nrow = ctr,
                      ncol = 5,
                      byrow = FALSE
                      )
  colnames(recording) <- c("St", "Et", "It", "Rt", "Time")

  return(list(matrix_time = output,
              matrix_record = recording
              ))
}
