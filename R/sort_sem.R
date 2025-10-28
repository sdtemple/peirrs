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
  if(dim(epi)[2]==6){
    # multitype model
    classes <- epi[,3][ind]
    ratesB <- epi[,4][ind]
    classesG <- epi[,5][ind]
    ratesG <- epi[,6][ind]
    return(cbind(i,r,classes,ratesB,classesG,ratesG))
  } else{
    return(cbind(i,r))
  }
}
