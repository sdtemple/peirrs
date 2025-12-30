#' Filter epidemic data
#'
#' Keep cases only.
#'
#' @param epidemic_time matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
filter_sem <- function(epidemic_time) {
  filter_by <- is.finite(epidemic_time[, 2]) | is.finite(epidemic_time[, 1])
  removals <- epidemic_time[, 2][filter_by]
  infections <- epidemic_time[, 1][filter_by]
  if (dim(epidemic_time)[2] == 6) {
    # multitype model
    infection_classes <- epidemic_time[, 3][filter_by]
    infection_rates <- epidemic_time[, 4][filter_by]
    removal_classes <- epidemic_time[, 5][filter_by]
    removal_rates <- epidemic_time[, 6][filter_by]
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
