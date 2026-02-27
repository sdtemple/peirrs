#' Filter epidemic data
#'
#' Keep infected individuals only.
#'
#' @param matrix_time matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @details
#' Non-infected individuals are denoted as infinite values in the matrix_time.
#' The matrix_time is either 2 columns if from the \code{simulator()} function
#' or 6 columns if from the \code{simulator_multitype()} function.
#'
#' @examples
#' # Basic SIR model
#' epidemic <- simulator(beta = 2, gamma = 1, population_size = 100)
#' infected_only <- filter_sem(epidemic$matrix_time)
#'
#' # Multi-type model
#' beta <- c(3, 2)
#' gamma <- c(1, 1)
#' infection_class_sizes <- c(50, 50)
#' removal_class_sizes <- c(50, 50)
#' epidemic_mt <- simulator_multitype(beta, gamma, infection_class_sizes, removal_class_sizes)
#' infected_only_mt <- filter_sem(epidemic_mt$matrix_time)
#'
#' @keywords internal
filter_sem <- function(matrix_time) {

  filter_by <- is.finite(matrix_time[, 2]) | is.finite(matrix_time[, 1])
  removals <- matrix_time[, 2][filter_by]
  infections <- matrix_time[, 1][filter_by]

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
