
#' Bayesian Inference for Multitype Epidemic Parameters using MCMC
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
#' @param binit Numeric initial estimates for the beta (transmission rates) parameters.
#'              Default is 1.
#' @param ginit Numeric initial estimates for the gamma (recovery rates) parameters.
#'              Default is 1.
#' @param bshape Numeric shape parameters for the gamma prior on beta.
#' @param gshape Numeric shape parameters for the gamma prior on gamma.
#' @param m A numeric shape parameter for the gamma distribution used in data augmentation.
#'          Default is 1.
#' @param num.time.update An integer specifying the number of Metropolis-Hastings updates
#'                        for infection and removal times per iteration. Default is 10.
#' @param num.iter An integer specifying the total number of MCMC iterations. Default is 500.
#' @param num.print An integer specifying the print frequency for iteration progress.
#'                  Default is 100.
#' @param num.tries An integer specifying the total number of draws
#'                  to check if proposal is consistent with an epidemic.
#'                  Default is 20.
#' @param update.gamma bool: TRUE to update removal rate estimates from initial estimate
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
peirr_bayes_multitype <- function(r,
                        i,
                        cr,
                        ci,
                        Ns,
                        binit,
                        ginit,
                        bshape,
                        gshape,
                        m=1,
                        num.time.update=5,
                        num.iter=500,
                        num.print=100,
                        num.tries=20,
                        update.gamma=FALSE,
                        lag=0
){

  ### some input checks ###

  removal.classes <- sort(unique(cr, na.rm=TRUE))
  infection.classes <- sort(unique(ci, na.rm=TRUE))
  if (length(removal.classes) != length(ginit)) {
    stop("Incorrect vector size of initial removal rates")
  }
  if (length(removal.classes) != length(gshape)) {
    stop("Incorrect vector size of removal rate shape parameters")
  }
  if (length(infection.classes) != length(binit)) {
    stop("Incorrect vector size of initial infection rates")
  }
  if (length(infection.classes) != length(bshape)) {
    stop("Incorrect vector size of infection rate shape parameters")
  }
  if (length(infection.classes) != length(Ns)) {
    stop("Incorrect vector size of infection rate class sizes")
  }
  num.removal.classes <- length(removal.classes)
  num.infection.classes <- length(infection.classes)

  ### set up initialization and prior ###

  # priors on beta and gamma come from initial estimate
  brate <- bshape / binit
  grate <- gshape / ginit
  b <- brate / Ns
  if (update.gamma) {
    g <- rgamma(num.removal.classes, shape=gshape, rate=grate)
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
  iupdate = function(r, i, ip, ci, Ns, bshape, brate, lag){

    infection.classes <- sort(unique(ci, na.rm=TRUE))
    num.infection.classes <- length(infection.classes)

    # initialize
    n = length(r)
    sum.ri <- r - i[1:n]
    sum.ri.p <- r - ip[1:n]

    if(length(ip) != length(i)){
      stop("Observed and proposed infection time vectors must have the same length.")
    }

    # compute tau matrices
    log.sum.tau <- 0
    for (u in 1:num.infection.classes) {
      cl = infection.classes[u]
      bool.cl = ci == cl
      num.type = sum(bool.cl)
      tau = matrix(0, nrow = n, ncol = n)
      for(j in (1:n)){
        tau[j,] = sapply(i[1:n][bool.cl] - lag, min, r[j]) - sapply(i[1:n][bool.cl] - lag, min, i[j])
      }
      log.sum.tau <- log.sum.tau + log(brate[u] + sum(tau) + (Ns[u] - num.type) * sum.ri) * (bshape[i] + num.type)
    }


    # compute tau matrices
    log.sum.tau.p <- 0
    for (u in 1:num.infection.classes) {
      cl = infection.classes[u]
      bool.cl = ci == cl
      num.type = sum(bool.cl)
      tau = matrix(0, nrow = n, ncol = n)
      for(j in (1:n)){
        tau[j,] = sapply(ip[1:n][bool.cl] - lag, min, rp[j]) - sapply(ip[1:n][bool.cl] - lag, min, ip[j])
      }
      log.sum.tau.p <- log.sum.tau.p + log(brate[u] + num.type + sum(tau) + (Ns[u] - num.type) * sum.ri.p)  * (bshape[i] + num.type)
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
    ellratio = ellratio + log.sum.tau - log.sum.tau.p
    return(ellratio)
  }

  # utlity for infection time metropolis hastings step
  rupdate = function(r, i, rp, ci, Ns, bshape, brate, lag){

    infection.classes <- sort(unique(ci, na.rm=TRUE))
    num.infection.classes <- length(infection.classes)

    # initialize
    n = length(r)
    sum.ri <- r - i[1:n]
    sum.ri.p <- rp - i[1:n]

    # check that rp has the same length as r
    if(length(rp) != length(r)){
      stop("Observed and proposed removal time vectors must have the same length.")
    }

    # compute tau matrices
    log.sum.tau <- 0
    for (u in 1:num.infection.classes) {
      cl = infection.classes[u]
      bool.cl = ci == cl
      num.type = sum(bool.cl)
      tau = matrix(0, nrow = n, ncol = n)
      for(j in (1:n)){
        tau[j,] = sapply(i[1:n][bool.cl] - lag, min, r[j]) - sapply(i[1:n][bool.cl] - lag, min, i[j])
      }
      log.sum.tau <- log.sum.tau + log(brate[u] + sum(tau) + (Ns[u] - num.type) * sum.ri) * (bshape[i] + num.type)
    }

    # compute tau matrices
    log.sum.tau.p <- 0
    for (u in 1:num.infection.classes) {
      cl = infection.classes[u]
      bool.cl = ci == cl
      num.type = sum(bool.cl)
      tau = matrix(0, nrow = n, ncol = n)
      for(j in (1:n)){
        tau[j,] = sapply(i[1:n][bool.cl] - lag, min, rp[j]) - sapply(i[1:n][bool.cl] - lag, min, i[j])
      }
      log.sum.tau.p <- log.sum.tau.p + log(brate[u] + num.type + sum(tau) + (Ns[u] - num.type) * sum.ri.p)  * (bshape[i] + num.type)
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
    ellratio = ellratio + log.sum.tau - log.sum.tau.p
    return(ellratio)
  }

  # start the sampler

  K = num.iter
  J = num.time.update
  U = num.print

  # initialize
  n = sum((!is.na(i)) | (!is.na(r)))
  storage = array(NA, dim = c((num.infection.classes + num.removal.classes + 2), K))
  tau = matrix(0, nrow = n, ncol = N)

  if(sum(!is.na(i)) == n && sum(!is.na(r)) == n){
    # premature exit
    # because complete data
    out <- bayes_complete_data_multitype(r,
                               i,
                               cr,
                               ci,
                               Ns,
                               infection.rate.rate.prior = b,
                               infection.rate.shape.prior = bshape,
                               removal.rate.rate.prior = grate,
                               removal.rate.shape.prior = gshape,
                               num.posterior.samples = num.iter,
                               lag=lag
                               )
    # work on this
    storage[1, ] <- out$infection.rate.samples
    storage[2, ] <- out$removal.rate.samples
    storage[3, ] <- rep(NA, K)
    storage[4, ] <- rep(NA, K)
    return(storage)
  }

  # first data augmentation
  nis <- rep(0, num.removal.classes)
  nrs <- rep(0, num.removal.classes)
  for (u in 1:num.removal.classes) {
    cl = removal.classes[u]
    nrc = sum(is.na(r) & (cr == cl))
    nic = sum(is.na(i) & (cr == cl))
    nrs[u] = nrc
    nis[u] = nic
  }
  ni = sum(nis)
  nr = sum(nrs)
  J1s = pmin(nis, J)
  J2s = pmin(nrs, J)
  ii = i
  ri = r
  for (u in 1:num.removal.classes) {
    cl = removal.classes[u]
    ii[is.na(i) & (cr == cl)] <- r[is.na(i) & (cr == cl)] - (rgamma(nis[u], shape = m, rate = 1) / g[u])
    ri[is.na(r) & (cr == cl)] <- i[is.na(r) & (cr == cl)] + (rgamma(nrs[u], shape = m, rate = 1) / g[u])
  }

  while(!is_epidemic(ri, ii, lag)){ # must be epidemic
    # can get hung for poorly drawn gamma
    if (update.gamma) {
      g <- rgamma(num.removal.classes, shape=gshape, rate=grate)
    } else {
      g <- ginit
    }
    for (u in 1:num.removal.classes) {
      cl = removal.classes[u]
      ii[is.na(i) & (cr == cl)] <- r[is.na(i) & (cr == cl)] - (rgamma(nis[u], shape = m, rate = 1) / g[u])
      ri[is.na(r) & (cr == cl)] <- i[is.na(r) & (cr == cl)] + (rgamma(nrs[u], shape = m, rate = 1) / g[u])
    }
  }
  ii = c(ii, rep(Inf, N - n))

  # sampling
  for(k in 1:K){

    # gamma metropolis hastings step
    if (update.gamma) {
      for (u in 1:num.removal.classes) {
        cl <- removal.classes[u]
        g[u] <- rgamma(1, gshape[u] + sum(cr==cl), grate[u] + sum((ri[cr==cl] - ii[1:n][cr==cl])))
      }
    }
    storage[(num.infection.classes+1):(num.infection.classes+num.removal.classes),k] = g

    # infection times metropolis hastings step
    successes <- 0
    if (ni > 0) { # can draw an infection time

      for (u in 1:num.removal.classes) {

        if (nis[u] > 0) { # can draw this type of infection time
          cl = removal.classes[u]
          J1 = J1s[u]

          for (j in 1:J1) {
            ctr = 1

            if (sum(is.na(i) & (cr == cl))==1) {
              l = (1:n)[is.na(i) & (cr == cl)]
            } else {
              l = sample((1:n)[is.na(i) & (cr == cl)], 1)
            }

            il = r[l] - (rgamma(1, shape = m, rate = 1) / g[u])
            ip = ii
            ip[l] = il

            while ((!is_epidemic(ri, ip[1:n], lag)) && (ctr <= num.tries)){
              ctr = ctr + 1

              if (sum(is.na(i) & (cr == cl))==1) {
                l = (1:n)[is.na(i) & (cr == cl)]
              } else {
                l = sample((1:n)[is.na(i) & (cr == cl)], 1)
              }

              il = r[l] - (rgamma(1, shape = m, rate = 1) / g[u])
              ip = ii
              ip[l] = il
            }

            if (is_epidemic(ri, ip[1:n], lag)) { # must be epidemic
              a = min(1, exp(iupdate(ri, ii, ip, ci, Ns, bshape, brate, lag=lag)))
              if(runif(1) < a){
                ii[l] = il
                successes <- successes + 1
              }
            }
          }
        }
      }
    }
    storage[(num.removal.classes + num.infection.classes + 1), k] <- successes / sum(J1s)

    # removal times metropolis hastings step
    successes <- 0
    if (nr > 0) {
      for (u in 1:num.removal.classes) {

        if (nrs[u] > 0) {
          cl = removal.classes[u]
          J2 = J2s[u]

          for(j in 1:J2){
            ctr = 1

            if (sum(is.na(r) & (cr == cl)) == 1) {
              l = (1:n)[is.na(r) & (cr == cl)]
            } else {
              l = sample((1:n)[is.na(r) & (cr == cl)], 1)
            }

            rl = i[l] + (rgamma(1, shape = m, rate = 1) / g[u])
            rp = ri
            rp[l] = rl

            while ((!is_epidemic(rp, ii[1:n], lag)) && (ctr <= num.tries)) {
              ctr = ctr + 1

              if (sum(is.na(r) & (cr == cl)) == 1) {
                l = (1:n)[is.na(r) & (cr == cl)]
              } else {
                l = sample((1:n)[is.na(r) & (cr == cl)], 1)
              }

              rl = i[l] + (rgamma(1, shape = m, rate = 1) / g[u])
              rp = ri
              rp[l] = rl
            }

            if (is_epidemic(rp, ii[1:n], lag)) { # must be epidemic
              a = min(1, exp(rupdate(ri, ii, rp, ci, Ns, bshape, brate, lag=lag)))
              if (runif(1) < a) {
                ri[l] = rl
                successes <- successes + 1
              }
            }
          }
        }
      }
    }
    storage[(num.removal.classes + num.infection.classes + 2), k] <- successes / sum(J2s)

    # beta gibbs step
    sum.ri <- ri - ii[1:n]
    for (u in 1:num.infection.classes) {
      cl = infection.classes[u]
      bool.cl = ci == cl
      num.type = sum(bool.cl)
      tau = matrix(0, nrow = n, ncol = n)
      for(j in (1:n)){
        tau[j,] = sapply(ii[1:n][bool.cl] - lag, min, ri[j]) - sapply(ii[1:n][bool.cl] - lag, min, ii[j])
      }
      b[u] = rgamma(1, shape = bshape[u] + num.type, rate = brate[u] + sum(tau) + (Ns[u] - num.type) * sum.ri)
    }
    storage[1:num.infection.classes,k] = b * N


    # iteration update
    if(!(k %% U)){
      print(k)
    }
  }

  return(storage)

}


