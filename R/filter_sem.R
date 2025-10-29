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
                    byrow = FALSE)
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
