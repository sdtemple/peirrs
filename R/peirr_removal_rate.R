#' Removal rate estimator
#'
#' Estimate the removal rate as the inverse of the mean infectious period length.
#'
#' @param removals numeric vector: removal times
#' @param infections numeric vector: infection times
#'
#' @return float: Removal rate estimate
#'
#' @export
peirr_removal_rate <- function(removals, 
                                infections
                                ) {

  ind <- (!is.na(removals)) * (!is.na(infections))
  ind <- which(ind == 1, arr.ind=TRUE)
  removals <- removals[ind]
  infections <- infections[ind]
  complete_period <- removals - infections
  return(length(complete_period) / sum(complete_period))
}
