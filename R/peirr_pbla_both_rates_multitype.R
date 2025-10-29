#' PBLA estimators of infection rates and removal rates
#'
#' Estimate multiple infection and removal rates with PBLA
#'
#' @param r numeric vector: removal times
#' @param N integer: population size
#' @param cr numeric vector: removal time classes
#' @param ci numeric vector: infection time classes
#' @param m positive integer shape
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return numeric list (infection.rates, removal.rates)
#'
#' @export
peirr_pbla_both_rates_multitype <- function(r,
                                  cr,
                                  ci,
                                  m=1,
                                  A=1,
                                  lag=0){

  # jointly optimize the likelihood
  num.rates <- length(unique(cr,na.rm=TRUE)) + length(unique(ci,na.rm=TRUE))
  num.beta.rates <- length(unique(ci,na.rm=TRUE))

  n <- sum(is.finite(r))
  N <- length(r)
  betamap <- matrix(0, nrow=n, ncol=N)
  for(b in 1:N){
    betamap[,b] = ci[b]
  }
  gammamap = cr[is.finite(r)]
  r = r[is.finite(r)]

  etc = list(m=m,A=A,lag=lag)
  pbla.estimates <- nlm(pblas::pbla_multi,
                        rep(1, num.rates),
                        R=num.beta.rates,
                        pbla=pblas::pbla_std,
                        betamap=betamap,
                        gammamap=gammamap,
                        r=r,
                        etc=etc
                        )

  # estimate of removal rate
  gamma.estims <- pbla.estimates$estimate[(num.beta.rates+1):num.rates]

  # estimate of infection rate
  beta.estims <- pbla.estimates$estimate[1:num.beta.rates]

  return(list(infection.rates=beta.estims,
              removal.rates=gamma.estims))
}
