#' Simulate stochastic epidemic model with different infection and removal rates
#'
#' Draw infectious periods for the multi-type SEIR model.
#'
#' @param beta numeric vector: infection rates
#' @param gamma numeric: removal rates
#' @param infection_class_sizes integer: subpopulation sizes for infection rates
#' @param removal_class_sizes integer: subpopulation sizes for removal rates
#' @param num_renewals integer: positive erlang shape
#' @param lag numeric: fixed incubation period
#'
#' @details
#' This function implements a multitype SIR epidemic model using an event-driven
#' (Gillespie) algorithm. The key feature is that individuals have separate
#' infection-class and removal-class assignments, each drawn from a stratified
#' population structure.
#'
#' The `beta` vector gives class-specific infection rates; `gamma` vector gives
#' class-specific removal rates. When a susceptible is infected, they are assigned
#' an infection class (determining their infectiousness) and independently a
#' removal class (determining how quickly they recover).
#'
#' Population structure:
#' \itemize{
#'   \item `infection_class_sizes`: sizes of infection-rate subpopulations
#'   \item `removal_class_sizes`: sizes of removal-rate subpopulations
#'   \item Both must sum to the same total (the epidemic population size)
#' }
#'
#' During the epidemic, individuals are assigned infection and removal classes
#' with probability proportional to class sizes. The output matrix includes
#' infection_class, infection_rate, removal_class, and removal_rate in addition
#' to infection and removal times.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `matrix_time`: an N x 6 matrix with columns infection, removal,
#'     infection_class, infection_rate, removal_class, removal_rate
#'   \item `matrix_record`: a T x 5 matrix with time-indexed columns St, Et, It, Rt, Time
#'     recording the susceptible, exposed, infectious, and removed counts plus elapsed time
#' }
#'
#' @examples
#' # Multitype epidemic with 2 infection classes and 2 removal classes
#' set.seed(1)
#' epi1 <- simulate_sem_multitype(beta = c(1.5, 2.0), gamma = c(0.8, 1.2),
#'                                infection_class_sizes = c(50, 50),
#'                                removal_class_sizes = c(50, 50))
#' head(epi1$matrix_time)
#'
#' # Multitype epidemic with unequal class sizes
#' set.seed(2)
#' epi2 <- simulate_sem_multitype(beta = c(1.2, 2.5), gamma = c(0.7, 1.5),
#'                                infection_class_sizes = c(70, 30),
#'                                removal_class_sizes = c(60, 40),
#'                                lag = 0.5, num_renewals = 2)
#' # Compare infectious period distributions by removal class
#' class1_periods <- epi2$matrix_time[epi2$matrix_time[, "removal_class"] == 1, "removal"] -
#'                   epi2$matrix_time[epi2$matrix_time[, "removal_class"] == 1, "infection"]
#' class2_periods <- epi2$matrix_time[epi2$matrix_time[, "removal_class"] == 2, "removal"] -
#'                   epi2$matrix_time[epi2$matrix_time[, "removal_class"] == 2, "infection"]
#' boxplot(list(class1 = class1_periods[is.finite(class1_periods)],
#'              class2 = class2_periods[is.finite(class2_periods)]))
#'
#' @keywords internal
simulate_sem_multitype <- function(beta,
                                   gamma,
                                   infection_class_sizes,
                                   removal_class_sizes,
                                   num_renewals = 1,
                                   lag = 0
                                   ) {

  if (sum(removal_class_sizes) != sum(infection_class_sizes)) {
    stop("Infection and removal rate classes do not add up")
  }

  # initialize vectors
  t <- 0
  population_size <- sum(infection_class_sizes)
  betaN <- beta / population_size
  infections <- rep(Inf, population_size)
  removals <- rep(Inf, population_size)
  infection_classes <- rep(NA, population_size)
  renewals <- rep(0, population_size)
  infection_weights <- infection_class_sizes / population_size
  zero_class_infection <- sample(1:length(infection_class_sizes), size=1, prob=infection_weights)
  infections[1] <- t
  infection_classes[1] <- zero_class_infection
  infection_class_sizes[zero_class_infection] <- infection_class_sizes[zero_class_infection] - 1
  itr <- 1

  infection_rates <- rep(NA, population_size)
  removal_rates <- rep(NA, population_size)
  removal_classes <- rep(NA, population_size)
  removals_current <- rep(0, length(removal_class_sizes))
  removal_weights <- removal_class_sizes / population_size
  zero_class_removal <- sample(1:length(removal_class_sizes), size=1, prob=removal_weights)
  removal_classes[1] <- zero_class_removal
  removal_class_sizes[zero_class_removal] <- removal_class_sizes[zero_class_removal] - 1
  infection_rates[1] <- beta[zero_class_infection]
  removal_rates[1] <- gamma[zero_class_removal]
  removals_current[zero_class_removal] <- removals_current[zero_class_removal] + 1

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
    min_time <- min(
      infections[is.infinite(removals) & is.finite(infections) & (infections > t)],
      Inf
    )
    # for updating the number recovering
    valid_indices <- which(is.infinite(removals) & is.finite(infections) & (infections > t))
    index_of_min_in_subset <- which.min(infections[valid_indices])
    arg_min_time <- valid_indices[index_of_min_in_subset]

    if (It == 0) {
      # no infecteds but there are exposeds
      # the closest exposure wait
      t <- min_time + .Machine$double.eps
      # and update the number recovering
      sampled_class_removal <- removal_classes[arg_min_time]
      removals_current[sampled_class_removal] <- removals_current[sampled_class_removal] + 1
    } else {
      # simulate time
      infection_rate <- It * sum(infection_class_sizes * betaN)
      removal_rate <- sum(removals_current * gamma)
      t <- t + rexp(1, rate = infection_rate + removal_rate)

      if (t > min_time) {
        # update time to make an exposed infectious
        t <- min_time + .Machine$double.eps
        # and update the number recovering
        sampled_class_removal <- removal_classes[arg_min_time]
        removals_current[sampled_class_removal] <- removals_current[sampled_class_removal] + 1
      } else {
        # there is infection or removal before
        # simulate transition
        x <- rbinom(1, size = 1, prob = removal_rate / (infection_rate + removal_rate))
        x = (x + 1) %% 2
        if (x) {
          # infect a susceptible
          infection_weights <- infection_class_sizes * betaN / (sum(infection_class_sizes * betaN))
          sampled_class_infection <- sample(1:length(infection_class_sizes), size=1, prob=infection_weights)
          infection_class_sizes[sampled_class_infection] <- infection_class_sizes[sampled_class_infection] - 1
          itr <- itr + 1
          infection_classes[itr] <- sampled_class_infection
          infections[itr] <- t + lag # fixed exposure period
          infection_rates[itr] <- beta[sampled_class_infection]

          # give the infected a removal rate
          weights_uniform <- removal_class_sizes / sum(removal_class_sizes)
          sampled_class_removal <- sample(1:length(removal_class_sizes), size=1, prob=weights_uniform)
          removal_class_sizes[sampled_class_removal] <- removal_class_sizes[sampled_class_removal] - 1
          removal_classes[itr] <- sampled_class_removal
          removal_rates[itr] <- gamma[sampled_class_removal]
          if (lag == 0) {
            # don't automatically update this
            # if there is an exposure delay
            removals_current[sampled_class_removal] <- removals_current[sampled_class_removal] + 1
          }

        } else {
          # remove an infected
          if (It > 1) {

            # sample based on removal rates
            removal_bool = which(is.infinite(removals) & is.finite(infections) & (infections <= t), arr.ind=T)
            removal_weights = removal_rates[removal_bool]
            removal_weights = removal_weights / sum(removal_weights)
            argx = sample(removal_bool, size=1, prob=removal_weights)

            # & (i <= t) means can't be removed before infectious when exposed
          } else {
            # when epidemic is winding down
            # no more infecteds
            argx = which(is.infinite(removals) & is.finite(infections) & (infections <= t))
            # & (i <= t) means can't be removed before infectious when exposed
          }
          renewals[argx] = renewals[argx] + 1
          if (renewals[argx] == num_renewals) { # after m renewals
            removals[argx] = t
            argx_class = removal_classes[argx]
            removals_current[argx_class] = removals_current[argx_class] - 1
          }
        }
      }
    }

    # update (S,I) counts
    St = sum(is.infinite(infections))
    It = sum(is.finite(infections) & (infections <= t)) - sum(is.finite(removals))
    Rt = sum(is.finite(infections) & is.finite(removals))
    Et = sum(is.finite(infections) & (infections > t))
    if (St + Rt + Et + It != population_size) {
      stop("S(t) + I(t) + E(t) + R(t) do not equal N")
    }
    # & (i <= t) delays the infectious period after exposure

    susceptible_recording = c(susceptible_recording, St)
    infection_recording = c(infection_recording, It)
    exposed_recording = c(exposed_recording, Et)
    removal_recording = c(removal_recording, Rt)
    time_recording = c(time_recording, t)
    ctr = ctr + 1

  }

  # assign the remaining non-infecteds to classes
  # infection classes
  if (itr < population_size) {
    remain_infection_classes <- c()
    remain_infection_rates <- c()
    for (l in 1:length(infection_class_sizes)) {
      if (infection_class_sizes[l] > 0) {
        remain_infection_classes <- c(remain_infection_classes, rep(l, infection_class_sizes[l]))
        remain_infection_rates <- c(remain_infection_rates, rep(betaN[l] * population_size, infection_class_sizes[l]))
      }
    }
    shuffled_indices <- sample(seq_along(remain_infection_classes))
    infection_classes[(itr+1):population_size] <- remain_infection_classes[shuffled_indices]
    infection_rates[(itr+1):population_size] <- remain_infection_rates[shuffled_indices]

    # removal classes
    remain_removal_classes <- c()
    remain_removal_rates <- c()
    for (l in 1:length(removal_class_sizes)) {
      if (removal_class_sizes[l] > 0) {
        remain_removal_classes <- c(remain_removal_classes, rep(l, removal_class_sizes[l]))
        remain_removal_rates <- c(remain_removal_rates, rep(gamma[l], removal_class_sizes[l]))
      }
    }
    shuffled_indices <- sample(seq_along(remain_removal_classes))
    removal_classes[(itr+1):population_size] <- remain_removal_classes[shuffled_indices]
    removal_rates[(itr+1):population_size] <- remain_removal_rates[shuffled_indices]
  }

  # there should be no negatives
  # and ignore the Inf values
  non_negative_check <- removals - infections
  non_negative_check <- non_negative_check[is.finite(non_negative_check)]
  if( any( (non_negative_check) < 0)){
    stop("At least one r - i value is negative")
  }

  # formatting
  output = matrix(c(infections,removals,infection_classes,infection_rates,removal_classes,removal_rates),
                  nrow = population_size,
                  ncol = 6,
                  byrow = FALSE)
  colnames(output) = c('infection',
                       'removal',
                       'infection_class',
                       'infection_rate',
                       'removal_class',
                       'removal_rate')

  recording = matrix(c(susceptible_recording, exposed_recording, infection_recording, removal_recording, time_recording),
                     nrow = ctr,
                     ncol = 5,
                     byrow = FALSE
  )
  colnames(recording) = c('St','Et','It', 'Rt', 'Time')

  return(list(matrix_time = output,
              matrix_record = recording
  )
  )
}
