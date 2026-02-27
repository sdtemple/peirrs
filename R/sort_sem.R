#' Sort epidemic data
#'
#' Sort matrix by increasing removal times.
#'
#' @param matrix_time matrix: infection times, removal times, (optional: infection and removal classes)
#'
#' @return matrix: infection times, removal times, (optional: infection and removal classes)
#'
#' @details
#' The matrix_time is either 2 columns if from the \code{simulator()} function
#' or 6 columns if from the \code{simulator_multitype()} function.
#' For 2-column matrices: columns are infection and removal times.
#' For 6-column matrices: columns are infection time, removal time, infection class,
#' infection rate, removal class, and removal rate.
#'
#' @examples
#' # Basic SIR model
#' epidemic <- simulator(beta = 2, gamma = 1, population_size = 100)
#' sorted_data <- sort_sem(epidemic$matrix_time)
#'
#' # Multi-type model
#' beta <- c(3, 2)
#' gamma <- c(1, 1)
#' infection_class_sizes <- c(50, 50)
#' removal_class_sizes <- c(50, 50)
#' epidemic_mt <- simulator_multitype(beta, gamma, infection_class_sizes, removal_class_sizes)
#' sorted_data_mt <- sort_sem(epidemic_mt$matrix_time)
#'
#' @keywords internal
sort_sem <- function(matrix_time) {
  removals <- matrix_time[, 2]
  infections <- matrix_time[, 1]
  ind <- order(removals)
  removals <- removals[ind]
  infections <- infections[ind]
  if (dim(matrix_time)[2] == 6) {
    # multitype model
    infection_classes <- matrix_time[, 3][ind]
    infection_rates <- matrix_time[, 4][ind]
    removal_classes <- matrix_time[, 5][ind]
    removal_rates <- matrix_time[, 6][ind]
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
                        byrow = FALSE
                        )
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
