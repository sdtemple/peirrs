#' Simulate general stochastic epidemic model with different infection rates
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param betas numeric vector: infection rates
#' @param gammas numeric: removal rates
#' @param beta.sizes integers: subpopulation sizes for infection rates
#' @param gamma.sizes integers: subpopulation sizes for removal rates
#' @param m integer: positive shape
#' @param lag numeric: fixed exposure period
#'
#' @return numeric list: matrix of (infection times, removal times, classes and rates), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulate_sem_multitype <- function(betas,
                                   gammas,
                                   beta.sizes,
                                   gamma.sizes,
                                   m = 1,
                                   lag = 0) {

  # initialize vectors
  t <- 0
  Ns <- beta.sizes
  betaNs <- betas / sum(Ns)
  N <- sum(Ns)
  i <- rep(Inf, N)
  r <- rep(Inf, N)
  classes <- rep(NA, N)
  M <- rep(0, N)
  weights <- Ns / N
  zeroclass <- sample(1:length(Ns), size=1, prob=weights)
  i[1] <- t
  classes[1] <- zeroclass
  Ns[zeroclass] <- Ns[zeroclass] - 1
  itr <- 1
  e <- lag

  Ms <- gamma.sizes
  ratesB <- rep(NA, N)
  ratesG <- rep(NA, N)
  classesG <- rep(NA, N)
  if(sum(gamma.sizes) != sum(beta.sizes)){
    stop("Infection and removal rate classes do not add up")
  }
  currentsG <- rep(0, length(Ms))
  weightsG <- Ms / N
  zeroclassG <- sample(1:length(Ms), size=1, prob=weightsG)
  classesG[1] <- zeroclassG
  Ms[zeroclassG] <- Ms[zeroclassG] - 1
  ratesB[1] <- betas[zeroclass]
  ratesG[1] <- gammas[zeroclassG]
  currentsG[zeroclassG] <- currentsG[zeroclassG] + 1

  # simulate epidemic
  St <- sum(is.infinite(i))
  It <- sum(is.finite(i)) - sum(is.finite(r))
  Et <- 0
  Rt <- 0

  # recording the evolution
  Srecording <- c(St)
  Irecording <- c(It)
  Erecording <- c(Et)
  Rrecording <- c(Rt)
  Trecording <- c(0)
  ctr <- 1

  while ( (It > 0) || (Et > 0) ) {

    # closest infectious time after exposure
    min.time <- min(
      i[is.infinite(r) & is.finite(i) & (i > t)],
      Inf
    )
    # for updating the number recovering
    valid_indices <- which(is.infinite(r) & is.finite(i) & (i > t))
    index_of_min_in_subset <- which.min(i[valid_indices])
    arg.min.time <- valid_indices[index_of_min_in_subset]

    if(It == 0){
      # no infecteds but there are exposeds
      # the closest exposure wait
      t <- min.time + .Machine$double.eps
      # and update the number recovering
      sampled.classG <- classesG[arg.min.time]
      currentsG[sampled.classG] <- currentsG[sampled.classG] + 1
    } else{
      # simulate time
      irate <- It * sum(Ns * betaNs)
      rrate <- sum(currentsG * gammas)
      t <- t + rexp(1, rate = irate + rrate)

      if(t > min.time) {
        # update time to make an exposed infectious
        t <- min.time + .Machine$double.eps
        # and update the number recovering
        sampled.classG <- classesG[arg.min.time]
        currentsG[sampled.classG] <- currentsG[sampled.classG] + 1
      } else {
        # there is infection or removal before
        # simulate transition
        x <- rbinom(1, size = 1, prob = rrate / (irate + rrate))
        x = (x + 1) %% 2
        if(x){
          # infect a susceptible
          weights <- Ns * betaNs / (sum(Ns * betaNs))
          sampled.class <- sample(1:length(Ns), size=1, prob=weights)
          Ns[sampled.class] <- Ns[sampled.class] - 1
          itr <- itr + 1
          classes[itr] <- sampled.class
          i[itr] <- t + e # fixed exposure period
          ratesB[itr] <- betas[sampled.class]

          # give the infected a removal rate
          weightsG <- Ms / sum(Ms)
          sampled.classG <- sample(1:length(Ms), size=1, prob=weightsG)
          Ms[sampled.classG] <- Ms[sampled.classG] - 1
          classesG[itr] <- sampled.classG
          ratesG[itr] <- gammas[sampled.classG]
          if(e==0){
            # don't automatically update this
            # if there is an exposure delay
            currentsG[sampled.classG] <- currentsG[sampled.classG] + 1
          }

        } else{
          # remove an infected
          if(It > 1){

            # sample based on removal rates
            removal.bool = which(is.infinite(r) & is.finite(i) & (i <= t), arr.ind=T)
            removal.weights = ratesG[removal.bool]
            removal.weights = removal.weights / sum(removal.weights)
            argx = sample(removal.bool, size=1, prob=removal.weights)

            # & (i <= t) means can't be removed before infectious when exposed
          } else{
            # when epidemic is winding down
            # no more infecteds
            argx = which(is.infinite(r) & is.finite(i) & (i <= t))
            # & (i <= t) means can't be removed before infectious when exposed
          }
          M[argx] = M[argx] + 1
          if(M[argx] == m){ # after m renewals
            r[argx] = t
            argx.class = classesG[argx]
            currentsG[argx.class] = currentsG[argx.class] - 1
          }
        }
      }
    }

    # update (S,I) counts
    St = sum(is.infinite(i))
    It = sum(is.finite(i) & (i <= t)) - sum(is.finite(r))
    Rt = sum(is.finite(i) & is.finite(r))
    Et = sum(is.finite(i) & (i > t))
    if(St + Rt + Et + It != N){
      stop("S(t) + I(t) + E(t) + R(t) do not equal N")
    }
    # & (i <= t) delays the infectious period after exposure

    Srecording = c(Srecording, St)
    Irecording = c(Irecording, It)
    Erecording = c(Erecording, Et)
    Rrecording = c(Rrecording, Rt)
    Trecording = c(Trecording, t)
    ctr = ctr + 1

  }

  # assign the remaining non-infecteds to classes
  # infection classes
  next.class = Rt
  for(l in 1:length(Ns)){
    if(Ns[l] > 0){
      classes[(next.class+1):(next.class+Ns[l])] = l
      ratesB[(next.class+1):(next.class+Ns[l])] = betas[l]
      next.class = next.class + Ns[l]
    }
  }
  # removal classes
  next.class = Rt
  for(l in 1:length(Ms)){
    if(Ms[l]>0){
      classesG[(next.class+1):(next.class+Ms[l])] = l
      ratesG[(next.class+1):(next.class+Ms[l])] = gammas[l]
      next.class = next.class + Ms[l]
    }
  }

  # there should be no negatives
  # and ignore the Inf values
  ri.check <- r - i
  ri.check <- ri.check[is.finite(ri.check)]
  if( any( (ri.check) < 0)){
    stop("At least one r - i value is negative")
  }

  # formatting
  output = matrix(c(i,r,classes,ratesB,classesG,ratesG),
                  nrow = N,
                  ncol = 6,
                  byrow = FALSE)
  colnames(output) = c('i',
                       'r',
                       'infection.group',
                       'infection.rate',
                       'removal.group',
                       'removal.rate')

  recording = matrix(c(Srecording, Erecording, Irecording, Rrecording, Trecording),
                     nrow = ctr,
                     ncol = 5,
                     byrow = FALSE
  )
  colnames(recording) = c('St','Et','It', 'Rt', 'Time')

  return(list(matrix.time = output,
              matrix.record = recording
  )
  )
}
