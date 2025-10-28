#' Simulate general stochastic epidemic model with formatting
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param beta numeric: infection rate
#' @param gamma numeric: removal rate
#' @param N integer: population size
#' @param m integer: positive shape
#' @param e numeric: fixed exposure period
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' @param min.sample.size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulator <- function(beta, 
                      gamma, 
                      N, 
                      m = 1, 
                      e = 0,
                     p = 0,
                     q = 1,
                     min.sample.size = 10
                     ){
  sample.size <- 0
  while(sample.size < min.sample.size){
    epi = simulate_sem(beta, gamma, N, m, e)
    epi$matrix.time = filter_sem(epi$matrix.time)
    epi$matrix.time = decomplete_sem(epi$matrix.time, p=p, q=q)
    epi$matrix.time = sort_sem(epi$matrix.time)
    sample.size = dim(epi$matrix.time)[1]
  }
  return(epi)
}
