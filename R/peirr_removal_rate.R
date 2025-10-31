#' Removal rate estimator
#'
#' Estimate the removal rate as the inverse of the mean or median infection period.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param med bool: TRUE for median, FALSE for mean
#'
#' @return float: Removal rate estimate
#'
#' @export
peirr_removal_rate <- function(r, i, med=TRUE){
  ind <- (!is.na(r)) * (!is.na(i))
  ind <- which(ind == 1, arr.ind=TRUE)
  r <- r[ind]
  i <- i[ind]
  ri <- r - i
  if(med){
    return(1 / median(ri))
  } else{
    return(length(ri) / sum(ri))
  }
}
