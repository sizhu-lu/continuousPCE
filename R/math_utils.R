#' Matrix Multiply with Expansion
#'
#' Perform element-wise multiplication of two matrices `U` and `V`, with the
#' option to expand dimensions of either matrix in a specified way. The function
#' supports three expansion types: expanding `U`, expanding `V`, or expanding
#' both matrices to match their respective dimensions.
#'
#' @param U A numeric matrix of size (\eqn{n_{u}}, n) or of size (\eqn{n_{u}},
#'   \eqn{n_{v}}, n). When `expansion_type = "U" or "both"`, this matrix will be
#'   expanded.
#' @param V A numeric matrix of size (\eqn{n_{v}}, n) or of size (\eqn{n_{u}},
#'   \eqn{n_{v}}, n). When `expansion_type = "V" or "both"`, this matrix will be
#'   expanded.
#' @param expansion_type A character string specifying the type of expansion.
#'   One of:
#'   \itemize{
#'     \item "U": Expand the dimensions of `U` to match the number of rows in `V`.
#'     \item "V": Expand the dimensions of `V` to match the number of rows in `U`.
#'     \item "both" (default): Expand both `U` and `V` to match the respective dimensions.
#'   }
#'
#' @return A numeric array of size (\eqn{n_{u}}, \eqn{n_{v}}, n), where
#'   \eqn{n_{u}} and \eqn{n_{v}} are the dimensions of `U` and `V` respectively.
#'
#' @examples
#' # Example 1: Expanding "both"
#' U <- matrix(1:6, nrow = 2, ncol = 3)
#' V <- matrix(7:18, nrow = 4, ncol = 3)
#' result <- matrix_multiply_with_expansion(U, V, expansion_type = "both")
#'
#' @export
#'
matrix_multiply_with_expansion <- function(U, V, expansion_type = "both") {
  switch (expansion_type,
          "U" = {
            # U (n_u, n)
            # V (n_u, n_v, n)
            # return (n_u, n_v, n)
            dim(U) = c(dim(U)[1], 1, dim(U)[2])
            U = U[, rep(1, dim(V)[2]), ]
            return(U * V)
          },
          "V" = {
            # U (n_u, n_v, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
            dim(V) = c(1, dim(V))
            V = V[rep(1, dim(U)[1]), , ]
            return(U * V)
          },
          "both" = {
            # U (n_u, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
            dim(U) = c(dim(U)[1], 1, dim(U)[2])
            U = U[, rep(1, dim(V)[1]), ]
            dim(V) = c(1, dim(V))
            V = V[rep(1, dim(U)[1]), , ]
            return(U * V)
          })
}

#' Trapezoidal Integration for Vectors (Numerical Approximation)
#'
#' Perform trapezoidal numerical integration of `y` with respect to `x`, always
#' along the first axis (rows). This function computes the area under the curve
#' defined by `y` as a function of `x` using the trapezoidal rule.
#'
#' @param x A numeric vector representing the x-axis values, assumed to be
#'   ordered. The length of `x` should match the length of `y` (or be compatible
#'   if `y` is a matrix or array).
#' @param y A numeric vector, matrix (can be high dimensional), or array
#'   representing the y-axis values. If `y` is a matrix or array, the
#'   integration will be performed along the first axis (rows). So the first
#'   dimension of `y` should have the same length as `x`.
#'
#' @return A numeric value or vector (depending on the input) representing the
#'   result of the trapezoidal integration. If `y` is a matrix or array, the
#'   result will be a vector of integrated values along the first axis (rows).
#'
#' @examples
#' # Example 1: Simple trapezoidal integration for vectors
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 3, 4, 5, 6)
#' result <- trapz_vector(x, y)
#'
#' # Example 2: Integration with a matrix `y`
#' x <- c(1, 2, 3, 4, 5)
#' y <- matrix(1:15, nrow = 5, ncol = 3)
#' result <- trapz_vector(x, y)
#'
#' @export
#'
trapz_vector <- function(x, y) {
  idx = 2:length(x)
  if (is.null(dim(y))) {
    return(sum((x[idx] - x[idx - 1]) * (
      asub(y, idx, dims = 1) + asub(y, idx - 1, dims = 1)
    ) / 2))
  }
  return(colSums((x[idx] - x[idx - 1]) * (
    asub(y, idx, dims = 1) + asub(y, idx - 1, dims = 1)
  ) / 2))
}

#' Trapezoidal CDF (Cumulative Distribution Function) Calculation
#'
#' This function calculates the cumulative distribution function (CDF) of a set
#' of data points using the trapezoidal rule. The trapezoidal rule approximates
#' the integral of a function by summing up the areas of trapezoids formed
#' between consecutive points.
#'
#' @param x A numeric vector representing the independent variable values (e.g.,
#'   time or space). It should have a length of at least 2.
#' @param y A numeric vector or matrix representing the dependent variable
#'   values. If a vector is provided, it will be treated as a single series. If
#'   a matrix is provided, each column will be treated as a separate series of
#'   dependent variables.
#'
#' @return If `y` is a vector, returns a numeric vector containing the CDF
#'   values. If `y` is a matrix, returns a matrix where each column represents
#'   the CDF values for the corresponding series in `y`.
#'
#' @details The CDF is computed by numerically integrating `y` with respect to
#'   `x` using the trapezoidal rule. The integration is performed by summing the
#'   areas of trapezoids formed between consecutive `x` values. The function is
#'   vectorized to support both univariate and multivariate input.
#'
#' @examples
#' # Example with a simple vector input for y
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x)
#' cdf <- trapz_cdf(x, y)
#' plot(x, cdf, type = "l", main = "CDF using Trapezoidal Rule", xlab = "x", ylab = "CDF")
#'
#' # Example with a matrix input for y (multiple series)
#' y_matrix <- cbind(sin(x), cos(x))
#' cdf_matrix <- trapz_cdf(x, y_matrix)
#' matplot(x, cdf_matrix, type = "l", col = c("red", "blue"),
#' main = "CDF of Multiple Series", xlab = "x", ylab = "CDF")
#'
#' @seealso \code{\link{trapz_vector}} for basic trapezoidal rule integration.
#'
#' @export
#'
trapz_cdf <- function(x, y) {
  x = c(2 * x[1] - x[2], x)
  idx = 2:length(x)
  if (is.null(dim(y))) {
    y = c(0, y)
    # x (n_x)
    # y (n_x)
    # return (n_x)
    return(cumsum((x[idx] - x[idx - 1]) * (
      asub(y, idx, dims = 1) + asub(y, idx - 1, dims = 1)
    ) / 2))
  }
  y = rbind(rep(0, dim(y)[2]), y)
  # x (n_x)
  # y (n_x, n)
  # return (n_x, n)
  return(apply(((x[idx] - x[idx - 1]) * (
    asub(y, idx, dims = 1) + asub(y, idx - 1, dims = 1)
  ) / 2), 2, cumsum))
}
