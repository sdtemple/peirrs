#' Simulate a Distance Matrix with Specified Statistical Properties
#'
#' Generates a random distance matrix by creating coordinates in 2D space,
#' computing pairwise distances, applying a transformation function, and
#' scaling the resulting values to match target mean and standard deviation.
#'
#' @param population_size Positive integer. Number of points (observations) for which to
#'   generate distance matrix.
#' @param kernel Function. Transformation function to apply to scaled distances.
#'   Must have a corresponding inverse function.
#' @param inverse_kernel Function. Inverse of the transformation function \code{kernel}.
#'   Must satisfy \code{inverse_kernel(kernel(x)) ≈ x}.
#' @param mu Numeric. Target mean for the transformed distance values.
#'   Default is 0.9.
#' @param sigma Numeric. Target standard deviation for the transformed distance
#'   values. Must be positive. Default is 0.01.
#' @param method Character. Distance metric to use. Passed to \code{dist()}.
#'   Default is 'euclidean'.
#' @param runif_max Numeric. Upper bound for random coordinate generation.
#'   Default is 100.
#' @param scalar Numeric. Scaling factor applied to distances before
#'   transformation. Default is -0.05.
#' @param num_tries Positive integer. Maximum number of iterations to attempt
#'   finding a valid distance matrix. Default is 1000.
#'
#' @return A numeric matrix of pairwise distances with dimensions N × N,
#'   where the transformed values have mean \code{mu} and standard deviation
#'   \code{sigma}.
#'
#' @details
#' The function iteratively:
#' \enumerate{
#'   \item Generates N random 2D coordinates
#'   \item Computes pairwise distances using the specified metric
#'   \item Applies transformation h and scales to target mean/sd
#'   \item Stops when all scaled values are non-negative
#'   \item Applies inverse transformation to obtain final distance matrix
#' }
#'
#' @examples
#' \dontrun{
#' D <- simulate_distance_matrix(population_size=10, kernel=exp, inverse_kernel=log, mu=0.9, sigma=0.01)
#' }
#'
#' @export
simulate_distance_matrix <- function(population_size,
                         kernel=exp,
                         inverse_kernel=log,
                         mu=0.9,
                         sigma=0.01,
                         method='euclidean',
                         runif_max=100,
                         scalar=-0.05,
                         num_tries=1000
                         ) {

  # some double checks, especially for user-defined inverse
  if (!is.function(kernel) || !is.function(inverse_kernel)) {
    stop("kernel and inverse_kernel must be functions")
  }
  val <- try(inverse_kernel(kernel(1)), silent = TRUE)
  if (abs(val - 1) > 1e-16) {
    stop("kernel and inverse_kernel must be inverses of each other")
  }
  val <- try(inverse_kernel(kernel(2)), silent = TRUE)
  if (abs(val - 2) > 1e-16) {
    stop("kernel and inverse_kernel must be inverses of each other")
  }
  val <- try(inverse_kernel(kernel(3)), silent = TRUE)
  if (abs(val - 3) > 1e-16) {
    stop("kernel and inverse_kernel must be inverses of each other")
  }
  if (population_size <= 0 || !is.numeric(population_size)) {
    stop("population_size must be a positive number")
  }
  if (sigma <= 0) {
    stop("sigma must be positive")
  }

  condition <- TRUE
  ctr <- 0
  while (condition) {
    ctr <- ctr + 1
    coords <- data.frame(
      x = runif(population_size, min = 0, max = runif_max),
      y = runif(population_size, min = 0, max = runif_max)
    )
    D <- as.matrix(dist(coords, method = method))
    W_h <- kernel(scalar * D)
    # Target means and standard deviations
    target_mean <- mu
    target_sd <- sigma

    # Original means and standard deviations
    original_values <- W_h[lower.tri(W_h)]
    original_mean <- mean(original_values)
    original_sd <- sd(original_values)

    # The transformation preserves the structure
    # while adjusting the scale and center
    W_scaled <- target_mean + (W_h - original_mean) * (target_sd / original_sd)

    # We don't want to scale by negatives
    if (all(W_scaled >= 0)) {
      condition <- FALSE
    }
    if (ctr >= num_tries) {
      stop("Maximum number of tries reached without finding a valid distance matrix.")
    }

  }

  # undo the transformation
  D_new <- inverse_kernel(W_scaled) / scalar
  diag(D_new) <- 0
  D_object <- as.dist(D_new) # to ensure it's a valid distance matrix
  return(D_new)
}
