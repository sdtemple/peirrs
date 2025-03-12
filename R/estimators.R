# source("incompleteSEM/R/psiFormulas.R") # is this necessary after package built

source("incompleteSEM/R/tauFormulas.R") # is this necessary after package built
library(pblas)

#' MLE for Removal Rate
#' 
#' Compute the maximum likelihood estimate for the (global) removal rate.
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' 
#' @return float: Removal rate estimate
#'  
#' @export 
mle_removal_rate <- function(r, i){
  ind <- (!is.na(r)) * (!is.na(i))
  ind <- which(ind == 1, arr.ind=T)
  r <- r[ind]
  i <- i[ind]
  ri <- r - i
  return(length(ri) / sum(ri))
}

#' MLE for complete stochastic epidemic model
#'
#' Estimate infection and removal rates with complete epidemic observations.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#'
#' @return MLEs for (infection.rate, removal.rate, R0)
#'
#' @export
mle_complete = function(r, i, N){
  n = length(r)
  t = 0
  for(j in 1:n){
    t = t + sum(sapply(i, min, r[j]) - sapply(i, min, i[j]))
  }
  ri = sum(r - i)
  g = n / ri
  b = (n - 1) / (t + (N - n) * ri) * N
  return(list(infection.rate=b,removal.rate=g,R0=b/g,tausum=t))
}


#' Pair-based tau estimator of common infection and removal rates
#' 
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' 
#' @return numeric list (infection.rate, removal.rate, R0, tau.sum)
#'  
#' @export 
peirr_tau_moments <- function(r, i, N){

    # estimate of removal rate
    gamma.estim <- mle_removal_rate(r,i)

    # number of infected
    n <- sum(!is.na(r) | !is.na(i))

    ### first infected ###
    ralpha <- which.min(r)
    ialpha <- which.min(i)
    if(i[ialpha] < r[ralpha]){
    alpha <- ialpha
    } else{
    alpha <- ralpha
    }

    # sum up tau terms
    tau <- 0
    for(j in (1:n)[-alpha]){
        rj <- r[j]
        ij <- i[j]
        for(k in (1:n)[-j]){
            rk <- r[k]
            ik <- i[k]
            tau <- tau + E.tau(rk, rj, ik, ij, gamma.estim, gamma.estim)
            }
        }
    }

    # these are where we have full infectious period
    full.r <- r[(!is.na(r)) & (!is.na(i))]
    full.i <- i[(!is.na(r)) & (!is.na(i))]

    # maximizes give conditional expectations
    beta.estim <- (n - 1) / (tau + (N - n) / gamma.estim * n)

    return(list(infection.rate=beta.estim*N, 
    removal.rate=gamma.estim, 
    R0=beta.estim*N/gamma.estim, 
    tau.sum=tau))
}

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
peirr_informed_pbla <- function(r, i, N){
  
  # PBLA function with fixed removal rate
  pb <- function(beta.estim,pbla,gamma.estim,r,N){
    pbla_std_gsem(r,beta.estim,gamma.estim,N)
  }
  
  # estimate of removal rate
  gamma.estim <- mle.removal.rate(r,i)
  
  # maximizes give conditional expectations
  beta.estim <-   nlm(pb,1,pbla=pbla_std_gsem,
                      gamma.estim=gamma.estim,
                      r=r,N=N)$estimate  
  
  return(list(infection.rate=beta.estim, 
              removal.rate=gamma.estim, 
              R0=beta.estim/gamma.estim))
}

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
  pbla.estimates <- nlm(pbla_gsem,c(1,1),pbla=pbla_std_gsem,r=r,N=N)
  
  # estimate of removal rate
  gamma.estim <- pbla.estimates$estimate[2]
  
  # estimate of infection rate
  beta.estim <- pbla.estimates$estimate[1]
  
  return(list(infection.rate=beta.estim, 
              removal.rate=gamma.estim, 
              R0=beta.estim/gamma.estim))
}


# fill this in
peirr_bayes <- function(r,i,N){
  return(list(infection.rate=0,
              removal.rate=0,
              R0=0
              ))
}


# multitype implementations


