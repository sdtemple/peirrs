#' Pair-based tau estimator of multiple infection and removal rates
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param removals numeric: removal times
#' @param infections numeric: infection times
#' @param removal_classes numeric: removal time classes
#' @param infection_classes numeric: infection time classes
#' @param infection_class_sizes integer: population sizes for infection classes (sorted)
#' @param lag numeric: fixed exposure period
#' @param median_tau bool: use median imputation for tau if TRUE
#' @param median_gamma bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection_rates, removal_rates, etc)
#'
#' @export
peirr_tau_multitype <- function(removals, 
                                infections,
                                removal_classes, 
                                infection_classes,
                                infection_class_sizes,
                                lag = 0,
                                median_tau = FALSE,
                                median_gamma = FALSE
                                ) {

  # make sure one or the other is finite
  or.finite <- is.finite(removals) | is.finite(infections)
  removal_classes <- removal_classes[or.finite]
  infection_classes <- infection_classes[or.finite]
  infections <- infections[or.finite]
  removals <- removals[or.finite]

  median_scalar <- 1
  if (median_tau) {median_scalar <- log(2)}

  # number of infected
  epidemic_size <- sum(!is.na(removals) | !is.na(infections))
  if (sum(is.na(removals))>=epidemic_size) {
    stop("There are no complete case periods to estimate removal rates")
  }
  if (sum(is.na(infections))>=epidemic_size) {
    stop("There are no complete case periods to estimate removal rates")
  }

  ### first infected ###
  alpha_r <- which.min(removals)
  alpha_i <- which.min(infections)
  if (infections[alpha_i] < removals[alpha_r]) {
    alpha <- alpha_i
  } else{
    alpha <- alpha_r
  }

  # initialize class specific rates
  removal_classes_unique <- sort(unique(removal_classes, na.rm=TRUE))
  infection_classes_unique <- sort(unique(infection_classes, na.rm=TRUE))
  removal_rates <- c()
  infection_rates <- c()
  removal_sizes <- c()
  infection_sizes <- c()
  removal_full_sizes <- c()
  removal_partial_sums <- c()

  # compute class-specific removal rates
  for (removal_class in removal_classes_unique) {

    # find those in same class
    indicators <- which(removal_classes == removal_class, arr.ind=TRUE)
    removals_kept <- removals[indicators]
    infections_kept <- infections[indicators]
    # find those that are complete
    filter_by <- (!is.na(removals_kept)) & (!is.na(infections_kept))
    filter_by <- which(filter_by == 1, arr.ind=TRUE)
    removals_kept_v2 <- removals_kept[filter_by]
    infections_kept_v2 <- infections_kept[filter_by]

    # estimate with complete obs
    removal_full_sizes <- c(removal_full_sizes, length(removals_kept_v2))
    if (!median_gamma) {
      rate_estim <- length(removals_kept_v2) / sum(removals_kept_v2 - infections_kept_v2)
    } else {
      rate_estim <- 1 / median(removals_kept_v2 - infections_kept_v2) * log(2)
    }
    removal_rates <- c(removal_rates, rate_estim)
    # compute expected value for incomplete obs
    num_not_complete <- length(removals_kept) - length(filter_by)
    removal_partial_sum <- sum(removals_kept_v2 - infections_kept_v2) + num_not_complete / rate_estim * median_scalar
    removal_partial_sums <- c(removal_partial_sums, removal_partial_sum)

  }

  # these sum with some expected values for observed infections
  complete_period_sum <- sum(removal_partial_sums)

  # sum up tau terms
  tau_sums <- c()
  num_not_infecteds <- c()
  for (l in 1:length(infection_classes_unique)) {
    tau_sum <- 0
    infection_class <- infection_classes_unique[l]
    num_infected <- length(which(infection_classes == infection_class))
    num_not_infected <- infection_class_sizes[l] - num_infected
    num_not_infecteds <- c(num_not_infecteds, num_not_infected)

    if (infection_classes[alpha] == infection_class) {
      num_infected <- num_infected - 1
    }

    for (j in (1:epidemic_size)[-alpha]) {
      infection_class_j <- infection_classes[j]
      removal_class_j <- removal_classes[j]

      if (infection_class_j == infection_class) {
        removal_j <- removals[j]
        infection_j <- infections[j]
        # should be a scalar
        rate_j <- removal_rates[which(removal_classes_unique == removal_class_j)]

        for (k in (1:epidemic_size)[-j]) {
          removal_class_k <- removal_classes[k]
          removal_k <- removals[k]
          infection_k <- infections[k]
          # should be a scalar
          rate_k <- removal_rates[which(removal_classes_unique == removal_class_k)]
          tau_kj <- tau_moment(removal_k, removal_j, infection_k, infection_j, rate_k, rate_j, lag, median_tau)
          if (is.na(tau_kj)) {print(c(removal_k, removal_j, infection_k, infection_j, rate_k, rate_j))}
          tau_sum <- tau_sum + tau_kj
        }
      }
    }
    tau_sums <- c(tau_sums, tau_sum)
    rate_estim <- num_infected / (tau_sum + num_not_infected * complete_period_sum)
    infection_rates <- c(infection_rates, rate_estim)
  }

  return(list(infection_rate = infection_rates * sum(infection_class_sizes),
              removal_rate = removal_rates,
              tau_sum = tau_sums,
              not_infected_sum = removal_partial_sums,
              num_not_infected = num_not_infecteds,
              num_complete = removal_full_sizes
              ))
}
