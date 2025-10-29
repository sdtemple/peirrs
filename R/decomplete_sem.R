#' Missing epidemic data
#'
#' Impute NAs for infection and removal times
#'
#' @param epi infection and removal times
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
decomplete_sem <- function(epi, p = 0, q = 1){
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
  if(dim(epi)[2]==6){
    # multitype model
    classes <- epi[,3][is.finite(epi[,2]) | is.finite(epi[,1])]
    ratesB <- epi[,4][is.finite(epi[,2]) | is.finite(epi[,1])]
    classesG <- epi[,5][is.finite(epi[,2]) | is.finite(epi[,1])]
    ratesG <- epi[,6][is.finite(epi[,2]) | is.finite(epi[,1])]
    N = length(r)
    # formatting
    output = matrix(c(i,r,classes,ratesB,classesG,ratesG),
                    nrow = N,
                    ncol = 6,
                    byrow = F)
    colnames(output) = c('i',
                         'r',
                         'infection.group',
                         'infection.rate',
                         'removal.group',
                         'removal.rate')
    return(output)
  } else{
    return(cbind(i,r))
  }
}
