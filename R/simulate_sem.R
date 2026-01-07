#' Simulate general stochastic epidemic model
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param population_size integer: population size
#' @param num_renewals integer: positive shape
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
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
