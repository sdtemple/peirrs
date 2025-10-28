#' Simulate multi-type stochastic epidemic model with formatting
#'
#' Draw infectious periods for stochastic epidemic with different classes.
#'
#' @param betas numeric vector: infection rates
#' @param gamma numeric: removal rates
#' @param beta.sizes integers: subpopulation sizes for infection rates
#' @param gamma.sizes integers: subpopulation sizes for removal rates
#' @param m integer: positive shape
#' @param e numeric: fixed exposure period
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' @param min.sample.size integer
#'
#' @return numeric list: matrix of (infection times, removal times), matrix of (St, It, Et, Rt, Time)
#'
#' @export
simulator_multitype <- function(betas, 
                      gammas, 
                      beta.sizes,
                      gamma.sizes,
                      m = 1, 
                      e = 0,
                      p = 0,
                      q = 1,
                      min.sample.size = 10
                      ){
  sample.size <- 0
  while(sample.size < min.sample.size){
    epi = simulate_sem_multitype(betas, gammas, beta.sizes, gamma.sizes, m, e)
    epi$matrix.time = filter_sem(epi$matrix.time)
    epi$matrix.time = decomplete_sem(epi$matrix.time, p=p, q=q)
    epi$matrix.time = sort_sem(epi$matrix.time)
    sample.size = dim(epi$matrix.time)[1]
  }
  return(epi)
}
