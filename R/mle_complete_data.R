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
mle_complete_data = function(r, i, N){
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
