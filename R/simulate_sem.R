#' Simulate general stochastic epidemic model
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param N integer: population size
#' @param m integer: positive shape
#' @param e numeric: fixed exposure period
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulate_sem <- function(beta, gamma, N, m = 1, e = 0){
  
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
  Et = 0
  Rt = 0
  
  # recording the evolution
  Srecording = c(St)
  Irecording = c(It)
  Erecording = c(Et)
  Rrecording = c(Rt)
  Trecording = c(0)
  ctr = 1
  
  while((It > 0) | (Et > 0)){
    
    # closest infectious time after exposure
    min.time = min(
      i[is.infinite(r) & is.finite(i) & (i > t)],
      Inf
    )
    
    if(It == 0){
      # no infecteds but there are exposeds
      # the closest exposure wait
      t = min.time + .Machine$double.eps
    } else{
      # simulate time
      irate = betaN * It * St
      rrate = gamma * It
      t = t + rexp(1, rate = irate + rrate)

      if(t > min.time){
        # update time to make an exposed infectious
        t = min.time + .Machine$double.eps
      } else{
        # there is infection or removal before
        # simulate transition
        x = rbinom(1, size = 1, prob = rrate / (irate + rrate))
        x = (x + 1) %% 2
        if(x){
          # infect a susceptible
          argx = sample(which(is.infinite(i), arr.ind = T), 1)
          i[argx] = t + e # fixed exposure period
        } else{
          # remove an infected
          if(It > 1){
            argx = sample(which(is.infinite(r) & is.finite(i) & (i <= t), arr.ind = T), 1)
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
  
  # formatting
  output = matrix(c(i,r),
                  nrow = N,
                  ncol = 2,
                  byrow = F)
  colnames(output) = c('i','r')
  
  recording = matrix(c(Srecording, Irecording, Erecording, Rrecording, Trecording),
                     nrow = ctr,
                     ncol = 5,
                     byrow = F
                     )
  colnames(recording) = c('St','It','Et', 'Rt', 'Time')
  
  return(list(matrix.time = output,
              matrix.record = recording
              )
         )
}
