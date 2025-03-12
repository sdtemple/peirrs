#' Posterior samples for infection and removal rates given incomplete data
#' 
#' Gibbs sample the infection and removal rates 
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' 
#' @return numeric list of posterior (infection.rate, removal.rate)
peirr_bayes <- function(r,i,N){
  return(list(infection.rate=0,
              removal.rate=0,
              R0=0
  ))
}