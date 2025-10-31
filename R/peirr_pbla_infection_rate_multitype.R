#' PEIRR likelihood estimator of infection rates and MLE removal rates
#'
#' Estimate multiple removal rates with duration date and multiple infection rates with PBLA
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param cr numeric vector: removal time classes
#' @param ci numeric vector: infection time classes
#' @param m positive integer shape
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#' @param known.gammas numeric: removal rates
#' @param gamma.med bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection.rates, removal.rates)
#'
#' @export
peirr_pbla_infection_rate_multitype <- function(r,
                                      i,
                                      cr,
                                      ci,
                                      m=1,
                                      A=1,
                                      lag=0,
                                      known.gammas = NULL,
                                      gamma.med = FALSE
                                      ){

  # PBLA function with fixed removal rate
  pb <- function(beta.estims,pbla,gamma.estims,betamap,gammamap,r,m,A,lag){
    # define the betamap
    betamapped <- pblas::multitypes(beta.estims, betamap)
    betas <- betamapped / ncol(betamapped)
    # define the gamma map
    gammas <- gamma.estims[gammamap]
    return(pbla(r,betas,gammas,m,A,lag))
  }

  # determine unique categories
  unique.gammas = sort(unique(cr, na.rm=TRUE))
  unique.betas = sort(unique(ci, na.rm=TRUE))

  # estimate of removal rate
  if(is.null(known.gammas)){
    removal.classes <- unique.gammas
    removal.rates <- c()
    # compute class-specific removal rates
    for(cl in removal.classes){
      # find those in same class
      indicators <- which(cr == cl)
      rg <- r[indicators]
      ig <- i[indicators]
      # find those that are complete
      filters <- (!is.na(rg)) & (!is.na(ig)) & (is.finite(rg)) & (is.finite(ig))
      filters <- which(filters == 1)
      rg2 <- rg[filters]
      ig2 <- ig[filters]
      # estimate with complete obs
      if(!gamma.med){
        rate.estim <- length(rg2) / sum(rg2 - ig2)
      } else{
        rate.estim <- 1 / median(rg2 - ig2)
      }
      removal.rates <- c(removal.rates, rate.estim)
    }
    gamma.estims <- removal.rates
  } else{
    gamma.estims <- known.gammas
    if(length(known.gammas)!=length(unique.gammas)){
      stop("Not enough known removal rates")
    }
    if(any(known.gammas <= 0)){
      stop("At least one removal rate is not positive")
    }
  }

  # recode the categories if not 1,2,3,...
  for(g in 1:length(unique.gammas)){
    cr[cr == unique.gammas[g]] = g
  }
  for(b in 1:length(unique.betas)){
    ci[ci == unique.betas[b]] = b
  }

  n <- sum(is.finite(r))
  N <- length(r)
  betamap <- matrix(0, nrow=n, ncol=N)
  for(b in 1:N){
    betamap[,b] = ci[b]
  }
  gammamap = cr[is.finite(r)]
  r = r[is.finite(r)]


  # maximizes give conditional expectations
  beta.estims <-   nlm(pb,
                      rep(1,length(unique.betas)),
                      pbla=pblas::pbla_std,
                      gamma.estims=gamma.estims,
                      betamap=betamap,
                      gammamap=gammamap,
                      r=r,
                      m=m,
                      A=A,
                      lag=lag)$estimate

  return(list(infection.rates=beta.estims,
              removal.rates=gamma.estims))
}
