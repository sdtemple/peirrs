#' Insert missingness into epidemic data
#'
#' Randomly insert NAs for infection and removal times to simulate incomplete observation.
#'
#' @param matrix_time matrix: infection and removal times
#' @param prop_complete numeric: expected proportion of complete pairs observed
#' @param prop_infection_missing numeric: probability infection time missing
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @details
#' This function takes complete epidemic data and randomly inserts missing values (NAs)
#' to simulate incomplete observation scenarios.
#' The matrix_time is either 2 columns if from the \code{simulator()} function
#' or 6 columns if from the \code{simulator_multitype()} function.
#' For each individual, with probability \code{1 - prop_complete}, one time is set to NA.
#' When a time is set to NA, it is the infection time with probability \code{prop_infection_missing},
#' otherwise it is the removal time.
#'
#' @examples
#' # Basic SIR model
#' epidemic <- simulator(beta = 2, gamma = 1, population_size = 100)
#' incomplete_data <- decomplete_sem(epidemic$matrix_time, prop_complete = 0.5, prop_infection_missing = 0.8)
#'
#' # Multi-type model
#' beta <- c(3, 2)
#' gamma <- c(1, 1)
#' infection_class_sizes <- c(50, 50)
#' removal_class_sizes <- c(50, 50)
#' epidemic_mt <- simulator_multitype(beta, gamma, infection_class_sizes, removal_class_sizes)
#' incomplete_data_mt <- decomplete_sem(epidemic_mt$matrix_time, prop_complete = 0.6)
#'
#' @keywords internal
decomplete_sem <- function(matrix_time,
                            prop_complete = 0,
                            prop_infection_missing = 1) {
  filter_by <- is.finite(matrix_time[, 2]) | is.finite(matrix_time[, 1])
  removals <- matrix_time[, 2]
  infections <- matrix_time[, 1]
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
  if (dim(matrix_time)[2] == 6) {
    # multitype model
    infection_classes <- matrix_time[, 3][filter_by]
    infection_rates <- matrix_time[, 4][filter_by]
    removal_classes <- matrix_time[, 5][filter_by]
    removal_rates <- matrix_time[, 6][filter_by]
    population_size <- length(removals)
    # formatting
    output <- matrix(c(infections,
                        removals,
                        infection_classes,
                        infection_rates,
                        removal_classes,
                        removal_rates),
                        nrow = population_size,
                        ncol = 6,
                        byrow = FALSE)
    colnames(output) <- c("infection",
                          "removal",
                          "infection_class",
                          "infection_rate",
                          "removal_class",
                          "removal_rate"
                          )
    output
  } else {
    output <- cbind(infections, removals)
    colnames(output) <- c("infection", "removal")
    output
  }
}
