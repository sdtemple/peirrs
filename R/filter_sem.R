#' Filter epidemic data
#'
#' Keep cases only.
#'
#' @param epi matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
filter_sem <- function(epi) {
  filtline <- is.finite(epi[, 2]) | is.finite(epi[, 1])
  r <- epi[, 2][filtline]
  i <- epi[, 1][filtline]
  if (dim(epi)[2] == 6) {
    # multitype model
    classes <- epi[, 3][filtline]
    ratesB <- epi[, 4][filtline]
    classesG <- epi[, 5][filtline]
    ratesG <- epi[, 6][filtline]
    N <- length(r)
    # formatting
    output <- matrix(c(i, r, classes, ratesB, classesG, ratesG),
                     nrow = N,
                     ncol = 6,
                     byrow = FALSE)
    colnames(output) <- c("i",
                          "r",
                          "infection.group",
                          "infection.rate",
                          "removal.group",
                          "removal.rate")
    output
  } else {
    cbind(i, r)
  }
}
