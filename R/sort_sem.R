#' Sort epidemic data
#' 
#' Sort matrix by increasing removal times.
#' 
#' @param epi matrix: infection times, removal times, (optional: infection classes)
#' 
#' @return matrix: infection times, removal times, (optional: infection classes)
#' 
#' @export
sort_sem <- function(epi){
  r <- epi[,2]
  i <- epi[,1]
  ind <- order(r)
  r <- r[ind]
  i <- i[ind]
  if(dim(epi)[2]==3){
    # multitype model
    classes <- epi[,3][ind]
    return(cbind(i,r,classes))
  } else{
    return(cbind(i,r))
  }
}
