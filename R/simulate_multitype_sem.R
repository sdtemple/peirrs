#' Simulate general stochastic epidemic model with different infection rates
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param betas numeric vector: infection rates
#' @param gamma numeric: removal rate
#' @param Ns integer: subpopulation sizes
#' @param m integer: positive shape
#'
#' @return matrix: infection times, removal times, infection classes
#'
#' @export
simulate_multitype_sem <- function(betas, gamma, Ns, m = 1){
  
  
  # initialize vectors
  t = 0
  betaNs = betas / sum(Ns)
  N <- sum(Ns)
  i = rep(Inf, N)
  r = rep(Inf, N)
  classes = rep(NA, N)
  M = rep(0, N)
  weights = Ns / N
  zeroclass = sample(1:length(Ns), size=1, prob=weights)
  i[1] = t
  classes[1] = zeroclass
  Ns[zeroclass] = Ns[zeroclass] - 1
  
  # simulate epidemic
  It = sum(is.finite(i)) - sum(is.finite(r))
  itr = 1
  while(It > 0){
    
    # simulate time
    irate = It * sum(Ns * betaNs)
    rrate = gamma * It
    t = t + rexp(1, rate = irate + rrate)
    
    # simulate transition
    x = rbinom(1, size = 1, prob = irate / (irate + rrate))
    if(x){
      # infect a susceptible
      weights = Ns * betaNs / (sum(Ns * betaNs))
      sampled.class = sample(1:length(Ns), size=1, prob=weights)
      Ns[sampled.class] = Ns[sampled.class] - 1
      itr = itr + 1
      classes[itr] = sampled.class
      i[itr] = t
    } else{
      # remove an infected
      if(It > 1){
        argx = sample(which(is.infinite(r) & is.finite(i), arr.ind = T), 1)
      } else{
        argx = which(is.infinite(r) & is.finite(i))
      }
      M[argx] = M[argx] + 1
      if(M[argx] == m){ # after m renewals
        r[argx] = t
      }
    }
    
    # update (I) counts
    St = sum(is.infinite(i))
    It = sum(is.finite(i)) - sum(is.finite(r))
    
  }
  
  # formatting
  # need to add a class
  output = matrix(c(i,r,classes),
                  nrow = N,
                  ncol = 3,
                  byrow = F)
  colnames(output) = c('i','r','classes')
  return(output)
}
