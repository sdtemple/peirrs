#' Simulate general stochastic epidemic model
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param N integer: population size
#' @param m integer: positive shape
#'
#' @return matrix: infection times, removal times
#'
#' @export
simulate_sem <- function(beta, gamma, N, m = 1){
  
  # initialize vectors
  t = 0
  betaN = beta / N
  i = rep(Inf, N)
  r = rep(Inf, N)
  M = rep(0, N)
  alpha = sample(N, 1)
  i[alpha] = t
  
  # simulate epidemic
  St = sum(is.infinite(i))
  It = sum(is.finite(i)) - sum(is.finite(r))
  while(It > 0){
    
    # simulate time
    irate = betaN * It * St
    rrate = gamma * It
    t = t + rexp(1, rate = irate + rrate)
    
    # simulate transition
    x = rbinom(1, size = 1, prob = irate / (irate + rrate))
    if(x){
      # infect a susceptible
      argx = sample(which(is.infinite(i), arr.ind = T), 1)
      i[argx] = t
    } else{
      # remove an infected
      if(It > 1){
        argx = sample(which(is.infinite(r) & is.finite(i), arr.ind = T), 1)
      } else{
        # when epidemic is winding down
        # no more infecteds
        argx = which(is.infinite(r) & is.finite(i))
      }
      M[argx] = M[argx] + 1
      if(M[argx] == m){ # after m renewals
        r[argx] = t
      }
    }
    
    # update (S,I) counts
    St = sum(is.infinite(i))
    It = sum(is.finite(i)) - sum(is.finite(r))
    
  }
  
  # formatting
  output = matrix(c(i,r),
                  nrow = N,
                  ncol = 2,
                  byrow = F)
  colnames(output) = c('i','r')
  return(output)
}
