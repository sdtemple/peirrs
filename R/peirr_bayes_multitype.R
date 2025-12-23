
#' Bayesian Inference for Epidemic Parameters using MCMC
#'
#' Performs Bayesian inference for the transmission rate (beta) and recovery rate (gamma)
#' parameters of an epidemic model using Markov Chain Monte Carlo (MCMC) sampling.
#' Implements data augmentation for missing infection and removal times via Metropolis-Hastings
#' updates, with a Gibbs step for the beta parameter.
#'
#' @param r A numeric vector of removal/recovery times. NA values indicate unobserved times.
#' @param i A numeric vector of infection times. NA values indicate unobserved times.
#' @param cr numeric vector: removal time classes
#' @param ci numeric vector: infection time classes
#' @param Ns integer vector: population sizes for infection classes (sorted)
#' @param m A numeric shape parameter for the gamma distribution used in data augmentation.
#'          Default is 1.
#' @param binit A numeric initial value for the beta (transmission rate) parameter.
#'              Default is 1.
#' @param ginit A numeric initial value for the gamma (recovery rate) parameter.
#'              Default is 1.
#' @param bshape A numeric shape parameter for the gamma prior on beta. Default is 1.
#' @param brate A numeric rate parameter for the gamma prior on beta. Default is 1e-4.
#' @param gshape A numeric shape parameter for the gamma prior on gamma. Default is 1.
#' @param grate A numeric rate parameter for the gamma prior on gamma. Default is 1e-4.
#' @param num.time.update An integer specifying the number of Metropolis-Hastings updates
#'                        for infection and removal times per iteration. Default is 10.
#' @param num.iter An integer specifying the total number of MCMC iterations. Default is 500.
#' @param num.print An integer specifying the print frequency for iteration progress.
#'                  Default is 100.
#'
#' @return A numeric matrix with 2 rows and num.iter columns containing posterior samples.
#'         Row 1 contains samples for beta (transmission rate).
#'         Row 2 contains samples for gamma (recovery rate).
#'
#' @details
#' The function implements a data augmentation MCMC algorithm for epidemic models.
#' It alternates between: (1) updating gamma via Gibbs sampling,
#' (2) updating missing infection times via Metropolis-Hastings,
#' (3) updating missing removal times via Metropolis-Hastings, and
#' (4) updating beta via Gibbs sampling. All epidemic configurations are validated
#' to ensure consistency with epidemic dynamics.
#' @export
peirr_bayes_multitype <- function(r,
                        i,
                        cr,
                        ci,
                        Ns,
                        m=1,
                        binit=1,
                        ginit=1,
                        bshape=1,
                        gshape=1,
                        num.time.update=10,
                        num.iter=500,
                        num.print=100
){

  ### set up initialization and prior ###

  # unique classes
  ur <- unique(cr)
  ui <- unique(ci)
  vr <- length(ur)
  vi <- length(ui)

  # cardinality of the infected removal classes
  cs <- c()
  for (u in ur) {
    cs <- c(cs, length(cr[cr == u]))
  }
  # cardinality of the infected infected classes
  ds <- c()
  for (u in ui) {
    ds <- c(ds, length(ci[ci == u]))
  }
  # cardinality of the not infected infected classes
  dstars <- Ns - ds

  # priors on beta and gamma come from initial estimate
  brate <- bshape / binit
  grate <- gshape / ginit
  b <- brate / Ns
  g <- rgamma(vr, shape=gshape, rate=grate)

  ### utility function local to ###

  # indicates data consistent with epidemic
  is_epidemic = function(r, i){
    n = length(r)
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < i[j]) * (r > i[j])
    }
    x = apply(ind, 1, sum)
    x = x[x == 0]
    if(length(x) > 1){
      return(FALSE)
    } else{
      return(TRUE)
    }
  }

  # utlity for infection time metropolis hastings step
  iupdate = function(r, i, ip, bshape, brate){

    # initialize
    n = length(r)
    N = length(i)

    # compute tau matrices
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
    }
    taup = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      taup[j,] = sapply(ip, min, r[j]) - sapply(ip, min, ip[j])
    }

    # compute
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < i[j]) * (r > i[j])
    }
    Blong = apply(ind, 1, sum)
    Blong = Blong[Blong > 0]
    indp = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      indp[j,] = (ip[1:n] < ip[j]) * (r > ip[j])
    }
    Bplong = apply(indp, 1, sum)
    Bplong = Bplong[Bplong > 0]

    ellratio = sum(log(Bplong)) - sum(log(Blong))
    ellratio = ellratio + (bshape + n - 1) *
      (log(brate + sum(tau)) - log(brate + sum(taup)))
    return(ellratio)
  }

  # utlity for infection time metropolis hastings step
  rupdate = function(r, i, rp, bshape, brate){

    # initialize
    n = length(r)
    N = length(i)

    # compute tau matrices
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i, min, r[j]) - sapply(i, min, i[j])
    }
    taup = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      taup[j,] = sapply(i, min, rp[j]) - sapply(i, min, i[j])
    }

    # compute
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < i[j]) * (r > i[j])
    }
    Blong = apply(ind, 1, sum)
    Blong = Blong[Blong > 0]
    indp = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      indp[j,] = (i[1:n] < i[j]) * (rp > i[j])
    }
    Bplong = apply(indp, 1, sum)
    Bplong = Bplong[Bplong > 0]

    ellratio = sum(log(Bplong)) - sum(log(Blong))
    ellratio = ellratio + (bshape + n - 1) *
      (log(brate + sum(tau)) - log(brate + sum(taup)))
    return(ellratio)
  }

  # start the sampler

  K = num.iter
  J = num.time.update
  U = num.print

  # initialize
  n = sum((!is.na(i)) | (!is.na(r)))
  storage = array(NA, dim = c(2, K))
  tau = matrix(0, nrow = n, ncol = N)

  # first data augmentation
  ni = sum(is.na(i))
  nr = sum(is.na(r))
  ii = i
  ri = r
  ii[is.na(i)] <- r[is.na(i)] - (rgamma(ni, shape = m, rate = 1) / g)
  ri[is.na(r)] <- i[is.na(r)] + (rgamma(nr, shape = m, rate = 1) / g)
  while(!is_epidemic(ri, ii)){ # must be epidemic
    ii[is.na(i)] <- r[is.na(i)] - (rgamma(ni, shape = m, rate = 1) / g)
    ri[is.na(r)] <- i[is.na(r)] + (rgamma(nr, shape = m, rate = 1) / g)
  }
  ii = c(ii, rep(Inf, N - n))
  ri = c(ri, rep(Inf, N - n))

  # sampling
  for(k in 1:K){

    # gamma metropolis hastings step
    g = rgamma(1, gshape + n, grate + sum((ri - ii)[is.finite(ri)]))
    storage[2,k] = g

    # infection times metropolis hastings step
    for(j in 1:J){
      l = sample((1:n)[is.na(i)], 1)
      il = r[l] - (rgamma(1, shape = m, rate = 1) / g)
      ip = ii
      ip[l] = il
      if(is_epidemic(ri, ip)){ # must be epidemic
        a = min(1, iupdate(ri, ii, ip, gshape, grate))
        if(runif(1) < a){
          ii[l] = il
        }
      }
    }

    # removal times metropolis hastings step
    for(j in 1:J){
      l = sample((1:n)[is.na(i)], 1)
      rl = i[l] + (rgamma(1, shape = m, rate = 1) / g)
      rp = ri
      rp[l] = rl
      if(is_epidemic(rp, ii)){ # must be epidemic
        a = min(1, rupdate(ri, ii, rp, gshape, grate))
        if(runif(1) < a){
          ri[l] = rl
        }
      }
    }

    # beta gibbs step
    for(j in 1:n){
      tau[j,] = sapply(ii, min, ri[j]) - sapply(ii, min, ii[j])
    }
    b = rgamma(1, shape = bshape + n - 1, rate = brate + sum(tau))
    storage[1,k] = b * N

    # iteration update
    if(!(k %% U)){
      print(k)
    }
  }

  return(storage)

}


