#' PEIRR likelihood estimator of common infection rate and estimated removal rate
#' 
#' Estimate removal rate with duration date and infection rate with PBLA 
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' 
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'  
#' @export 
peirr_pbla_infection_rate <- function(r, i, N){
  
  # PBLA function with fixed removal rate
  pb <- function(beta.estim,pbla,gamma.estim,r,N){
    pbla(r,beta.estim,gamma.estim,N)
  }
  
  # estimate of removal rate
  gamma.estim <- mle.removal.rate(r,i)
  
  # maximizes give conditional expectations
  beta.estim <-   nlm(pb,
                      1,
                      pbla=pblas::pbla_std_gsem,
                      gamma.estim=gamma.estim,
                      r=r,
                      N=N)$estimate  
  
  return(list(infection.rate=beta.estim, 
              removal.rate=gamma.estim, 
              R0=beta.estim/gamma.estim))
}
