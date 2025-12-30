
#' Bayesian Inference for Epidemic Parameters using MCMC
#'
#' Performs Bayesian inference for the transmission rate (beta) and recovery rate (gamma)
#' parameters of an epidemic model using Markov Chain Monte Carlo (MCMC) sampling.
#' Implements data augmentation for missing infection and removal times via Metropolis-Hastings
#' updates, with a Gibbs step for the beta parameter.
#'
#' @param r A numeric vector of removal/recovery times. NA values indicate unobserved times.
#' @param i A numeric vector of infection times. NA values indicate unobserved times.
#' @param N An integer specifying the total population size.
#' @param m A numeric shape parameter for the gamma distribution used in data augmentation.
#'          Default is 1.
#' @param binit A numeric initial estimate for the beta (transmission rate) parameter.
#'              Default is 1.
#' @param ginit A numeric initial estimate for the gamma (recovery rate) parameter.
#'              Default is 1.
#' @param bshape A numeric shape parameter for the gamma prior on beta. Default is 1.
#' @param gshape A numeric shape parameter for the gamma prior on gamma. Default is 1.
#' @param num.time.update An integer specifying the number of Metropolis-Hastings updates
#'                        for infection and removal times per iteration. Default is 10.
#' @param num.iter An integer specifying the total number of MCMC iterations. Default is 500.
#' @param num.print An integer specifying the print frequency for iteration progress.
#'                  Default is 100.
#' @param num.tries An integer specifying the total number of draws
#'                  to check if proposal is consistent with an epidemic.
#'                  Default is 20.
#' @param update.gamma bool: TRUE to update removal rate estimate from initial estimate
#'                  Default is FALSE.
#' @param lag numeric: fixed exposure period
#'
#' @return A numeric matrix with 2 rows and num.iter columns containing posterior samples.
#'         Row 1 contains samples for beta (transmission rate).
#'         Row 2 contains samples for gamma (recovery rate).
#'         Row 3 contains proportion of infection times augmented/updated.
#'         Row 4 contains proportion of removal times augmented/updated.
#'
#' @details
#' The function implements a data augmentation MCMC algorithm for epidemic models.
#' It alternates between: (1) updating gamma via Gibbs sampling,
#' (2) updating missing infection times via Metropolis-Hastings,
#' (3) updating missing removal times via Metropolis-Hastings, and
#' (4) updating beta via Gibbs sampling. All epidemic configurations are validated
#' to ensure consistency with epidemic dynamics.
#' @export
peirr_bayes <- function(r,
                        i,
                        N,
                        m=1,
                        binit=1,
                        ginit=1,
                        bshape=1,
                        gshape=1,
                        num.time.update=10,
                        num.iter=500,
                        num.print=100,
                        num.tries=20,
                        update.gamma=FALSE,
                        lag=0
){

  ### set up initialization and prior ###

  # priors on beta and gamma come from initial estimate
  brate <- bshape / binit
  grate <- gshape / ginit
  b <- brate / N
  if (update.gamma) {
    g <- rgamma(1, shape=gshape, rate=grate)
  } else {
    g <- ginit
  }

  ### utility function local to ###

  # indicates data consistent with epidemic
  is_epidemic = function(r, i, lag){
    n = length(r)
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < (i[j] - lag)) * (r > (i[j] - lag))
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
  iupdate = function(r, i, ip, bshape, brate, lag){

    # initialize
    n = length(r)
    N = length(i)

    if(length(ip) != length(i)){
      stop("Observed and proposed infection time vectors must have the same length.")
    }

    # compute tau matrices
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i - lag, min, r[j]) - sapply(i - lag, min, i[j])
    }
    taup = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      taup[j,] = sapply(ip - lag, min, r[j]) - sapply(ip - lag, min, ip[j])
    }

    # compute
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < (i[j] - lag)) * (r > (i[j] - lag))
    }
    # L1 comes from page 25 of the stockdale 2019 thesis
    # Integrated out formula is page 32 of my prelim
    # There are typos in my exam, though
    # Assuming gamma density function for nice proposal
    L1.stockdale19.long = apply(ind, 1, sum)
    L1.stockdale19.long = L1.stockdale19.long[L1.stockdale19.long > 0]

    indp = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      indp[j,] = (ip[1:n] < (ip[j] - lag)) * (r > (ip[j] - lag))
    }
    L1.stockdale19.long.p = apply(indp, 1, sum)
    L1.stockdale19.long.p = L1.stockdale19.long.p[L1.stockdale19.long.p > 0]

    ellratio = sum(log(L1.stockdale19.long.p)) - sum(log(L1.stockdale19.long))
    ellratio = ellratio + (bshape + n - 1) *
      (log(brate + sum(tau)) - log(brate + sum(taup)))
    return(ellratio)
  }

  # utlity for infection time metropolis hastings step
  rupdate = function(r, i, rp, bshape, brate, lag){

    # initialize
    n = length(r)
    N = length(i)

    # check that rp has the same length as r
    if(length(rp) != length(r)){
      stop("Observed and proposed removal time vectors must have the same length.")
    }

    # compute tau matrices
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(i - lag, min, r[j]) - sapply(i - lag, min, i[j])
    }
    taup = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      taup[j,] = sapply(i - lag, min, rp[j]) - sapply(i - lag, min, i[j])
    }

    # compute
    ind = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      ind[j,] = (i[1:n] < (i[j] - lag)) * (r > (i[j] - lag))
    }
    # L1 comes from page 25 of the stockdale 2019 thesis
    # Integrated out formula is page 32 of my prelim
    # There are typos in my exam, though
    # Assuming gamma density function for nice proposal
    L1.stockdale19.long = apply(ind, 1, sum)
    L1.stockdale19.long = L1.stockdale19.long[L1.stockdale19.long > 0]

    indp = matrix(0, nrow = n, ncol = n)
    for(j in 1:n){
      indp[j,] = (i[1:n] < (i[j] - lag)) * (rp > (i[j] - lag))
    }
    L1.stockdale19.long.p = apply(indp, 1, sum)
    L1.stockdale19.long.p = L1.stockdale19.long.p[L1.stockdale19.long.p > 0]

    ellratio = sum(log(L1.stockdale19.long.p)) - sum(log(L1.stockdale19.long))
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
  storage = array(NA, dim = c(4, K))
  tau = matrix(0, nrow = n, ncol = N)

  if(sum(!is.na(i)) == n && sum(!is.na(r)) == n){
    # premature exit
    # because complete data
    out <- bayes_complete_data(r,
                               i,
                               N,
                               infection.rate.rate.prior = b,
                               infection.rate.shape.prior = bshape,
                               removal.rate.rate.prior = grate,
                               removal.rate.shape.prior = gshape,
                               num.posterior.samples = num.iter,
                               lag=lag
                               )
    storage[1, ] <- out$infection.rate.samples
    storage[2, ] <- out$removal.rate.samples
    storage[3, ] <- rep(NA, K)
    storage[4, ] <- rep(NA, K)
    return(storage)
  }

  # first data augmentation
  ni = sum(is.na(i))
  nr = sum(is.na(r))
  J1 = min(ni, J)
  J2 = min(nr, J)
  ii = i
  ri = r
  ii[is.na(i)] <- r[is.na(i)] - (rgamma(ni, shape = m, rate = 1) / g)
  ri[is.na(r)] <- i[is.na(r)] + (rgamma(nr, shape = m, rate = 1) / g)
  while(!is_epidemic(ri, ii, lag)){ # must be epidemic
    # can get hung for poorly drawn gamma
    if (update.gamma) {
      g <- rgamma(1, shape=gshape, rate=grate)
    } else {
      g <- ginit
    }
    ii[is.na(i)] <- r[is.na(i)] - (rgamma(ni, shape = m, rate = 1) / g)
    ri[is.na(r)] <- i[is.na(r)] + (rgamma(nr, shape = m, rate = 1) / g)
  }
  ii = c(ii, rep(Inf, N - n))

  # sampling
  for(k in 1:K){

    # gamma metropolis hastings step
    if (update.gamma) {
      g = rgamma(1, gshape + n, grate + sum((ri - ii[1:n])))
    }
    storage[2,k] = g

    # infection times metropolis hastings step
    successes <- 0
    if (ni > 0) {
      for(j in 1:J1){
        ctr = 1
        if (sum(is.na(i))==1) {
          l = (1:n)[is.na(i)]
        } else {
          l = sample((1:n)[is.na(i)], 1)
        }
        il = r[l] - (rgamma(1, shape = m, rate = 1) / g)
        ip = ii
        ip[l] = il
        while ((!is_epidemic(ri, ip[1:n], lag)) && (ctr <= num.tries)){
          ctr = ctr + 1
          if (sum(is.na(i))==1) {
            l = (1:n)[is.na(i)]
          } else {
            l = sample((1:n)[is.na(i)], 1)
          }
          il = r[l] - (rgamma(1, shape = m, rate = 1) / g)
          ip = ii
          ip[l] = il
        }
        if (is_epidemic(ri, ip[1:n], lag)) { # must be epidemic
          a = min(1, exp(iupdate(ri, ii, ip, bshape, brate, lag=lag)))
          if(runif(1) < a){
            ii[l] = il
            successes <- successes + 1
          }
        }
      }
    }
    storage[3, k] <- successes / J1

    # removal times metropolis hastings step
    successes <- 0
    if (nr > 0) {
      for(j in 1:J2){
        ctr = 1
        if (sum(is.na(r))==1) {
          l = (1:n)[is.na(r)]
        } else {
          l = sample((1:n)[is.na(r)], 1)
        }
        rl = i[l] + (rgamma(1, shape = m, rate = 1) / g)
        rp = ri
        rp[l] = rl
        while ((!is_epidemic(rp, ii[1:n], lag)) && (ctr <= num.tries)) {
          ctr = ctr + 1
          if (sum(is.na(r))==1) {
            l = (1:n)[is.na(r)]
          } else {
            l = sample((1:n)[is.na(r)], 1)
          }
          rl = i[l] + (rgamma(1, shape = m, rate = 1) / g)
          rp = ri
          rp[l] = rl
        }
        if (is_epidemic(rp, ii[1:n], lag)) { # must be epidemic
          a = min(1, exp(rupdate(ri, ii, rp, bshape, brate, lag=lag)))
          if (runif(1) < a) {
            ri[l] = rl
            successes <- successes + 1
          }
        }
      }
    }
    storage[4, k] <- successes / J2

    # beta gibbs step
    tau = matrix(0, nrow = n, ncol = N)
    for(j in 1:n){
      tau[j,] = sapply(ii - lag, min, ri[j]) - sapply(ii - lag, min, ii[j])
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


