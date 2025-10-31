#' PEIRR likelihood estimator of infection rate with MLE removal rate
#'
#' Estimate removal rate with duration date and infection rate with PBLA
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param m positive integer shape
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#' @param known.gamma numeric: removal rate
#' @param gamma.med bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'
#' @export
peirr_pbla_infection_rate <- function(r,
                                      i,
                                      N,
                                      m=1,
                                      A=1,
                                      lag=0,
                                      known.gamma = NULL,
                                      gamma.med = FALSE
                                      ){

  # PBLA function with fixed removal rate
  pb <- function(beta.estim,pbla,gamma.estim,r,N,m,A,lag){
    return(pbla(r,beta.estim,gamma.estim,N,m,A,lag))
  }

  # estimate of removal rate
  if(is.null(known.gamma)){
    gamma.estim <- peirr_removal_rate(r,i,gamma.med)
  } else{
    gamma.estim <- known.gamma
    if(length(known.gamma)>1){
      stop("More than 1 removal rate")
    }
    if(known.gamma <= 0){
      stop("Removal rate is not positive")
    }
  }

  # maximizes give conditional expectations
  beta.estim <-   nlm(pb,
                      1,
                      pbla=pblas::pbla_std_gsem,
                      gamma.estim=gamma.estim,
                      r=r,
                      N=N,
                      m=m,
                      A=A,
                      lag=lag)$estimate

  return(list(infection.rate=beta.estim,
              removal.rate=gamma.estim,
              R0=beta.estim/gamma.estim))
}
