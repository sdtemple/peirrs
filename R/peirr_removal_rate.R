#' Removal rate estimator
#'
#' Estimate the removal rate as the inverse of the mean or median infection period.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#' @param median_gamma bool: TRUE for median, FALSE for mean
#'
#' @return float: Removal rate estimate
#'
#' @export
peirr_removal_rate <- function(removals, 
                                infections, 
                                median_gamma = TRU
                                ) {

  ind <- (!is.na(removals)) * (!is.na(infections))
  ind <- which(ind == 1, arr.ind=TRUE)
  removals <- removals[ind]
  infections <- infections[ind]
  complete_period <- removals - infections
  
  if (median_gamma) {
    return(1 / median(complete_period) * log(2))

  } else {
    return(length(complete_period) / sum(complete_period))
  }
}
