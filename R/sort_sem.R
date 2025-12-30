#' Sort epidemic data
#'
#' Sort matrix by increasing removal times.
#'
#' @param epidemic_time matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
sort_sem <- function(epidemic_time) {
  removals <- epidemic_time[, 2]
  infections <- epidemic_time[, 1]
  ind <- order(removals)
  removals <- removals[ind]
  infections <- infections[ind]
  if (dim(epidemic_time)[2] == 6) {
    # multitype model
    infection_classes <- epidemic_time[, 3][ind]
    infection_rates <- epidemic_time[, 4][ind]
    removal_classes <- epidemic_time[, 5][ind]
    removal_rates <- epidemic_time[, 6][ind]
    population_size <- length(removals)
    # formatting
    output <- matrix(c(infections, removals, infection_classes, infection_rates, removal_classes, removal_rates),
                     nrow = population_size,
                     ncol = 6,
                     byrow = FALSE)
    colnames(output) <- c("infection",
                          "removal",
                          "infection_class",
                          "infection_rate",
                          "removal_class",
                          "removal_rate")
    output
  } else {
    output <- cbind(infections, removals)
    colnames(output) <- c("infection", "removal")
    output
  }
}
