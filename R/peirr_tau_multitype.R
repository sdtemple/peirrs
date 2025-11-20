#' Pair-based tau estimator of multiple infection and removal rates
#'
#' Estimate infection and removal rates with tau-based expectation-maximization.
#' The output value \code{tau.sum} is useful for debugging.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param cr numeric vector: removal time classes
#' @param ci numeric vector: infection time classes
#' @param Ns integer vector: population sizes for infection classes (sorted)
#' @param lag numeric: fixed exposure period
#' @param tau.med bool: use median imputation for tau if TRUE
#' @param gamma.med bool: TRUE for median, and FALSE for mean in estimating the removal rate
#'
#' @return numeric list (infection.rates, removal.rates, removal.full.sizes)
#'
#' @export
peirr_tau_multitype <- function(r, i,
                                cr, ci,
                                Ns,
                                lag=0,
                                tau.med=TRUE,
                                gamma.med=FALSE
                                ){

  # make sure one or the other is finite
  or.finite <- is.finite(r)|is.finite(i)
  cr <- cr[or.finite]
  ci <- ci[or.finite]
  i <- i[or.finite]
  r <- r[or.finite]

  median.scalar <- 1
  if(tau.med){median.scalar <- log(2)}

  # number of infected
  n <- sum(!is.na(r) | !is.na(i))
  if(sum(is.na(r))>=n){
    stop("There are no complete case periods to estimate removal rates")
  }
  if(sum(is.na(i))>=n){
    stop("There are no complete case periods to estimate removal rates")
  }

  ### first infected ###
  ralpha <- which.min(r)
  ialpha <- which.min(i)
  if(i[ialpha] < r[ralpha]){
    alpha <- ialpha
  } else{
    alpha <- ralpha
  }

  # initialize class specific rates
  removal.classes <- sort(unique(cr, na.rm=TRUE))
  infection.classes <- sort(unique(ci, na.rm=TRUE))
  removal.rates <- c()
  infection.rates <- c()
  removal.sizes <- c()
  infection.sizes <- c()
  removal.full.sizes <- c()
  removal.partial.sums <- c()

  # compute class-specific removal rates
  for(cl in removal.classes){

    # find those in same class
    indicators <- which(cr == cl, arr.ind=TRUE)
    rg <- r[indicators]
    ig <- i[indicators]
    # find those that are complete
    filters <- (!is.na(rg)) & (!is.na(ig))
    filters <- which(filters == 1, arr.ind=TRUE)
    rg2 <- rg[filters]
    ig2 <- ig[filters]

    # estimate with complete obs
    removal.full.sizes <- c(removal.full.sizes, length(rg2))
    if(!gamma.med){
      rate.estim <- length(rg2) / sum(rg2 - ig2)
    } else{
      rate.estim <- 1 / median(rg2 - ig2)
    }
    removal.rates <- c(removal.rates, rate.estim)
    # compute expected value for incomplete obs
    num.not.complete <- length(rg) - length(filters)
    removal.partial.sum <- sum(rg2 - ig2) + num.not.complete / rate.estim * median.scalar
    removal.partial.sums <- c(removal.partial.sums, removal.partial.sum)

  }

  # these sum with some expected values for observed infections
  ri.sum <- sum(removal.partial.sums)

  # sum up tau terms
  taus <- c()
  for(l in 1:length(infection.classes)){
    tau <- 0
    cl <- infection.classes[l]
    numer <- length(which(ci == cl))
    not.denom <- Ns[l] - numer
    if(ci[alpha] == cl){
      numer <- numer - 1
    }
    for(j in (1:n)[-alpha]){
      cj <- ci[j]
      if(cj == cl){
        rj <- r[j]
        ij <- i[j]
        # should be a scalar
        vj <- removal.rates[which(removal.classes == cj)]
        for(k in (1:n)[-j]){
          ck <- ci[k]
          rk <- r[k]
          ik <- i[k]
          # should be a scalar
          vk <- removal.rates[which(removal.classes == ck)]
          tm <- tau_moment(rk, rj, ik, ij, vk, vj, lag, tau.med)
          if(is.na(tm)){print(c(rk,rj,ik,ij,vk,vj))}
          tau <- tau + tm
        }
      }
    }
    taus <- c(taus, tau)
    rate.estim <- numer / (tau + not.denom * ri.sum)
    infection.rates <- c(infection.rates, rate.estim)
  }

  return(list(infection.rates=infection.rates*sum(Ns),
              removal.rates=removal.rates,
              removal.full.sizes=removal.full.sizes,
              tau.sums=taus,
              full.ri.sum=ri.sum,
              num.not.infected=sum(Ns)-n
              ))
}
