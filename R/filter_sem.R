#' Filter epidemic data
#'
#' Keep cases only.
#'
#' @param matrix_time matrix: infection times, removal times, (optional: infection classes)
#'
#' @return matrix: infection times, removal times, (optional: infection classes)
#'
#' @export
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
