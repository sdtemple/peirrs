#' Filter epidemic data
#'
#' Keep cases only.
#'
#' @param epi matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
filter_sem <- function(epi){
  r <- epi[,2][is.finite(epi[,2]) | is.finite(epi[,1])]
  i <- epi[,1][is.finite(epi[,2]) | is.finite(epi[,1])]
  if(dim(epi)[2]==3){
    # multitype model
    classes <- epi[,3][is.finite(epi[,2]) | is.finite(epi[,1])]
    return(cbind(i,r,classes))
  } else{
    return(cbind(i,r))
  }
}
