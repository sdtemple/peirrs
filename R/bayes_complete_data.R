#' Posterior parameters for infection and removal rates given complete data
#' 
#' Parameters for independent Gibbs sampling of the infection and removal rates 
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' @param infection.rate.rate.prior numeric
#' @param infection.rate.shape.prior numeric
#' @param removal.rate.rate.prior numeric
#' @param removal.rate.shape.prior numeric
#' 
#' @return numeric list of exact posterior and prior parameters
bayes_complete_data <- function(r,i,N,
                                infection.rate.rate.prior,
                                infection.rate.shape.prior,
                                removal.rate.rate.prior,
                                removal.rate.shape.prior
                                ){
  n = length(r)
  tau = matrix(0, nrow = n, ncol = N)
  for(j in 1:n){
    tau[j,1:n] = sapply(i, min, r[j]) - sapply(i, min, i[j])
  }
  tau[,(n+1):N] = r - i
  Ai = sum(tau)
  Ci = sum(r-i)
  # not exporting documentation
  # this function has limited interesting use
  return(list(infection.rate.rate=infection.rate.rate.prior+Ai,
              infection.rate.shape=infection.rate.shape.prior+n-1,
              infection.rate.rate.prior=infection.rate.rate.prior,
              infection.rate.shape.prior=infection.rate.shape.prior,
              removal.rate.rate=removal.rate.rate.prior+Ci,
              removal.rate.shape=removal.rate.shape.prior+n,
              removal.rate.rate.prior=removal.rate.rate.prior,
              removal.rate.shape.prior=removal.rate.shape.prior
  ))
}
