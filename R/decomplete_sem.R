#' Missing epidemic data
#' 
#' Impute NAs for infection and removal times
#' 
#' @param epi infection and removal times
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' 
#' @return matrix: infection times, removal times, (optional: infection claseses)
#' 
#' @export 
decomplete_sem <- function(epi, p, q = 1){
  r <- epi[,2]
  i <- epi[,1]
  n <- length(r)
  for(j in 1:n){
    if(rbinom(1, 1, 1 - p)){
      if(rbinom(1, 1, q)){
        i[j] <- NA
      } else{
        r[j] <- NA
      }
    }
  }
  if(dim(epi)[2]==3){
    # multitype model
    return(cbind(i,r,epi[,3]))
  } else{
    return(cbind(i,r))
  }
}