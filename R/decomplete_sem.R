#' Missing epidemic data
#'
#' Impute NAs for infection and removal times
#'
#' @param epidemic_time matrix: infection and removal times
#' @param prop_complete numeric: expected proportion of complete pairs observed
#' @param prop_infection_missing numeric: probability infection time missing
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
decomplete_sem <- function(epidemic_time, prop_complete = 0, prop_infection_missing = 1) {
  filter_by <- is.finite(epidemic_time[, 2]) | is.finite(epidemic_time[, 1])
  removals <- epidemic_time[, 2]
  infections <- epidemic_time[, 1]
  population_size <- length(removals)
  for (j in 1:population_size) {
    if (rbinom(1, 1, 1 - prop_complete)) {
      if (rbinom(1, 1, prop_infection_missing)) {
        infections[j] <- NA
      } else {
        removals[j] <- NA
      }
    }
  }
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
