#' MLE for Removal Rate
#'
#' Compute the maximum likelihood estimate for the (global) removal rate.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#'
#' @return float: Removal rate estimate
#'
#' @export
mle_removal_rate <- function(r, i){
  ind <- (!is.na(r)) * (!is.na(i))
  ind <- which(ind == 1, arr.ind=TRUE)
  r <- r[ind]
  i <- i[ind]
  ri <- r - i
  return(length(ri) / sum(ri))
}
