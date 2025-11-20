#' PBLA estimator of infection rate and removal rate
#' 
#' Estimate infection and removal rates with PBLA 
#' 
#' @param r numeric vector: removal times
#' @param N integer: population size
#' @param m positive integer shape
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#' 
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'  
#' @export 
peirr_pbla_both_rates <- function(r,
                                  N,
                                  m = 1,
                                  A = 1,
                                  lag = 0,
                                  i = NA
                                  ) {
  if (any(!is.na(i))) {
    stop("i is not all NAs. This (hidden) named parameter only exists for bootstrapping.")
  }
  # jointly optimize the likelihood
  etc <- list(m = m, A = A, lag = lag)
  pbla.estimates <- nlm(pblas::pbla_gsem,
                        c(1, 1),
                        pbla = pblas::pbla_std_gsem,
                        r = r,
                        N = N,
                        etc
                        )
  # estimate of removal rate
  gamma.estim <- pbla.estimates$estimate[2]
  # estimate of infection rate
  beta.estim <- pbla.estimates$estimate[1]
  return(list(infection.rate = beta.estim,
              removal.rate = gamma.estim, 
              R0 = beta.estim / gamma.estim))
}
