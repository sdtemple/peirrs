#' PEIRR likelihood estimator of infection rates and MLE removal rates
#'
#' Estimate multiple removal rates with duration date and multiple infection rates with PBLA
#'
#' @param removals numeric: removal times (as large a vector as population size, with NAs for uninfected)
#' @param infections numeric: infection times (as large a vector as population size, with NAs for uninfected)
#' @param removal_classes numeric: removal time classes (as large a vector as population size, with NAs for uninfected)
#' @param infection_classes numeric: infection time classes (as large a vector as population size, with NAs for uninfected)
#' @param num_renewals integer: erlang shape
#' @param num_patient_zeros integer: patient zeros
#' @param lag numeric: fixed lag
#' @param known_gamma numeric: removal rates
#' @param median_gamma bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection_rate, removal_rate)
#'
#' @export
peirr_pbla_infection_rate_multitype <- function(removals,
                                      infections,
                                      removal_classes,
                                      infection_classes,
                                      num_renewals = 1,
                                      num_patient_zeros = 1,
                                      lag = 0,
                                      known_gamma = NULL,
                                      median_gamma = FALSE
                                      ) {

  # PBLA function with fixed removal rate
  pb <- function(beta_estims, pbla, gamma_estims, beta_map, gamma_map, removals, num_renewals, num_patient_zeros, lag) {
    # define the betamap
    beta_mapped <- pblas::multitypes(beta_estims, beta_map)
    betas <- beta_mapped / ncol(beta_mapped)
    # define the gamma map
    gammas <- gamma_estims[gamma_map]
    return(pbla(removals, betas, gammas, num_renewals, num_patient_zeros, lag))
  }

  # determine unique categories
  unique_gammas <- sort(unique(removal_classes, na.rm=TRUE))
  unique_betas <- sort(unique(infection_classes, na.rm=TRUE))

  # estimate of removal rate
  if (is.null(known_gamma)) {
    removal_rates <- c()
    # compute class-specific removal rates
    for (removal_class_unique in unique_gammas) {
      # find those in same class
      indicators <- which(removal_classes == removal_class_unique)
      removals_kept <- removals[indicators]
      infections_kept <- infections[indicators]
      # find those that are complete
      filter_by <- (!is.na(removals_kept)) & 
        (!is.na(infections_kept)) & 
        (is.finite(removals_kept)) & 
        (is.finite(infections_kept)
        )
      filters <- which(filter_by == 1)
      removals_kept_v2 <- removals_kept[filters]
      infections_kept_v2 <- infections_kept[filters]
      # estimate with complete obs
      if(!median_gamma){
        rate_estim <- length(removals_kept_v2) / sum(removals_kept_v2 - infections_kept_v2)
      } else{
        rate_estim <- 1 / median(removals_kept_v2 - infections_kept_v2) * log(2)
      }
      removal_rates <- c(removal_rates, rate_estim)
    }
    gamma_estims <- removal_rates
  } else {
    gamma_estims <- known_gamma
    if (length(known_gamma) != length(unique_gammas)) {
      stop("Not enough known removal rates")
    }
    if (any(known_gamma <= 0)) {
      stop("At least one removal rate is not positive")
    }
  }

  # recode the categories if not 1,2,3,...
  for (g in 1:length(unique_gammas)) {
    removal_classes[removal_classes == unique_gammas[g]] <- g
  }
  for (b in 1:length(unique_betas)) {
    infection_classes[infection_classes == unique_betas[b]] <- b
  }

  epidemic_size <- sum(is.finite(removals))
  population_size <- length(removals)
  beta_map <- matrix(0, nrow=epidemic_size, ncol=population_size)
  for (b in 1:population_size) {
    beta_map[,b] <- infection_classes[b]
  }
  gamma_map <- removal_classes[is.finite(removals)]
  removals <- removals[is.finite(removals)]

  # maximizes give conditional expectations
  beta_estims <-   nlm(pb,
                      rep(1,length(unique_betas)),
                      pbla=pblas::pbla_std,
                      gamma_estims=gamma_estims,
                      beta_map=beta_map,
                      gamma_map=gamma_map,
                      removals=removals,
                      num_renewals=num_renewals,
                      num_patient_zeros=num_patient_zeros,
                      lag=lag)$estimate

  return(list(infection_rate=beta_estims,
              removal_rate=gamma_estims))
}
