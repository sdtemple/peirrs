#' PBLA estimator of common infection rate and estimated removal rate
#' 
#' Estimate infection and removal rates with PBLA 
#' 
#' @param r numeric vector: removal times
#' @param N integer: population size
#' 
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'  
#' @export 
peirr_naive_pbla <- function(r, N){
  
  # jointly optimize the likelihood
  pbla.estimates <- nlm(pbla_gsem,
                        c(1,1),
                        pbla=pblas::pbla_std_gsem,
                        r=r,
                        N=N)
  
  # estimate of removal rate
  gamma.estim <- pbla.estimates$estimate[2]
  
  # estimate of infection rate
  beta.estim <- pbla.estimates$estimate[1]
  
  return(list(infection.rate=beta.estim, 
              removal.rate=gamma.estim, 
              R0=beta.estim/gamma.estim))
}
