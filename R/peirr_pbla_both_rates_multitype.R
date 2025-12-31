#' PBLA estimators of infection rates and removal rates
#'
#' Estimate multiple infection and removal rates with PBLA
#'
#' @param removals numeric: removal times (as large a vector as population size, with NAs for uninfected)
#' @param removal_classes numeric: removal time classes (as large a vector as population size, with NAs for uninfected)
#' @param infection_classes numeric: infection time classes (as large a vector as population size, with NAs for uninfected)
#' @param num_renewals integer: erlang shape
#' @param num_patient_zeros integer: patient zeros
#' @param lag numeric: fixed lag
#'
#' @return numeric list (infection_rate, removal_rate_)
#'
#' @export
peirr_pbla_both_rates_multitype <- function(removals,
                                  removal_classes,
                                  infection_classes,
                                  num_renewals=1,
                                  num_patient_zeros=1,
                                  lag=0){

  # jointly optimize the likelihood
  num_rates <- length(unique(removal_classes,na.rm=TRUE)) + length(unique(infection_classes,na.rm=TRUE))
  num_beta_rates <- length(unique(infection_classes,na.rm=TRUE))

  epidemic_size <- sum(is.finite(removals))
  population_size <- length(removals)
  beta_map <- matrix(0, nrow=epidemic_size, ncol=population_size)
  for(b in 1:population_size){
    beta_map[,b] = infection_classes[b]
  }
  gamma_map = removal_classes[is.finite(removals)]
  removals = removals[is.finite(removals)]

  etc = list(m=num_renewals, A=num_patient_zeros, lag=lag)
  pbla.estimates <- nlm(pblas::pbla_multi,
                        rep(1, num_rates),
                        R=num_beta_rates,
                        pbla=pblas::pbla_std,
                        betamap=beta_map,
                        gammamap=gamma_map,
                        r=removals,
                        etc=etc
                        )

  # estimate of removal rate
  gamma_estims <- pbla.estimates$estimate[(num_beta_rates+1):num_rates]

  # estimate of infection rate
  beta_estims <- pbla.estimates$estimate[1:num_beta_rates]

  return(list(infection.rate=beta_estims,
              removal.rate=gamma_estims))
}
