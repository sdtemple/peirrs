#' Pair-based tau estimator of group-specific infection and removal rates
#'
#' Estimate group-specific infection and removal rates with tau-based expectation-maximization.
#'
#' @param removals numeric: removal times
#' @param infections numeric: infection times
#' @param removal_classes numeric: removal time classes
#' @param infection_classes numeric: infection time classes
#' @param infection_class_sizes integer: population sizes for infection classes (sorted)
#' @param lag numeric: fixed incubation period
#'
#' @details
#' This function extends \code{peirr_tau()} to estimate class-specific infection and
#' removal rates in multitype epidemic models. Individuals belong to both an infection
#' class (determining their infectiousness) and a removal class (determining their
#' recovery rate).
#'
#' Step 1: Estimate class-specific removal rates (gamma_l) using maximum likelihood
#' on complete infection-removal pairs within each removal class. For each class l:
#' \deqn{\hat{\gamma}_l = \frac{n_l}{\sum_{i \in class l} (r_i - i_i)}}
#'
#' Step 2: Estimate class-specific infection rates (beta_m) using expectation-maximization
#' with pairwise transmission indicators (tau). For each infection class m, the function
#' computes expected tau values for all pairs where the infector belongs to class m,
#' using the class-specific removal rates estimated in Step 1.
#'
#' Class vectors should be sorted and aligned with the input vectors. The
#' `infection_class_sizes` vector specifies the total population size for each
#' infection class (including both infected and susceptible individuals).
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `infection_rate`: vector of estimated class-specific beta values
#'   \item `removal_rate`: vector of estimated class-specific gamma values
#' }
#'
#' @examples
#' # Simulate multitype epidemic with 2 infection classes and 2 removal classes
#' set.seed(1)
#' epi1 <- simulator_multitype(beta = c(1.5, 2.0), gamma = c(0.8, 1.2),
#'                             infection_class_sizes = c(50, 50),
#'                             removal_class_sizes = c(50, 50),
#'                             prop_complete = 1)
#' 
#' fit1 <- peirr_tau_multitype(removals = epi1$matrix_time[, "removal"],
#'                             infections = epi1$matrix_time[, "infection"],
#'                             removal_classes = epi1$matrix_time[, "removal_class"],
#'                             infection_classes = epi1$matrix_time[, "infection_class"],
#'                             infection_class_sizes = c(50, 50),
#'                             lag = 0)
#' fit1$infection_rate  # Should be near c(1.5, 2.0) * 100
#' fit1$removal_rate    # Should be near c(0.8, 1.2)
#'
#' # Multitype epidemic with missing data
#' set.seed(2)
#' epi2 <- simulator_multitype(beta = c(1.2, 2.5), gamma = c(0.7, 1.5),
#'                             infection_class_sizes = c(60, 40),
#'                             removal_class_sizes = c(50, 50),
#'                             prop_complete = 0.75, prop_infection_missing = 0.6)
#' 
#' fit2 <- peirr_tau_multitype(removals = epi2$matrix_time[, "removal"],
#'                             infections = epi2$matrix_time[, "infection"],
#'                             removal_classes = epi2$matrix_time[, "removal_class"],
#'                             infection_classes = epi2$matrix_time[, "infection_class"],
#'                             infection_class_sizes = c(60, 40),
#'                             lag = 0)
#' fit2$infection_rate
#' fit2$removal_rate
#'
#' @export
peirr_tau_multitype <- function(removals, 
                                infections,
                                removal_classes, 
                                infection_classes,
                                infection_class_sizes,
                                lag = 0
                                ) {

  # make sure one or the other is finite
  or.finite <- is.finite(removals) | is.finite(infections)
  removal_classes <- removal_classes[or.finite]
  infection_classes <- infection_classes[or.finite]
  infections <- infections[or.finite]
  removals <- removals[or.finite]


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
    rate_estim <- length(removals_kept_v2) / sum(removals_kept_v2 - infections_kept_v2)

    removal_rates <- c(removal_rates, rate_estim)
    # compute expected value for incomplete obs
    num_not_complete <- length(removals_kept) - length(filter_by)
    removal_partial_sum <- sum(removals_kept_v2 - infections_kept_v2) + num_not_complete / rate_estim
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
          tau_kj <- tau_moment(removal_k, removal_j, infection_k, infection_j, rate_k, rate_j, lag)
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
              removal_rate = removal_rates
              ))
}
