#' Gaussian Copula Function
#'
#' This function calculates the value of a Gaussian copula given two uniform
#' random variables and a correlation coefficient. The Gaussian copula is used
#' to model the dependence structure between random variables, and it is defined
#' using the inverse cumulative distribution function (quantile function) of the
#' standard normal distribution.
#'
#' @param u A numeric vector or matrix of uniform random variables (u_i ∈ [0,
#'   1]). Each value in `u` represents a univariate uniform random variable.
#' @param v A numeric vector or matrix of uniform random variables (v_i ∈ [0,
#'   1]). Each value in `v` represents a univariate uniform random variable.
#' @param rho A numeric scalar between -1 and 1 representing the correlation
#'   parameter of the Gaussian copula.
#'
#' @return A numeric vector or matrix of copula values corresponding to the
#'   input values of `u` and `v`. The return value has the same dimensions as
#'   the inputs `u` and `v`.
#'
#' @details The Gaussian copula captures the dependency between random
#' variables. It is parameterized by a correlation coefficient `rho`, which
#' should lie between -1 and 1. The function computes the value of the Gaussian
#' copula based on the formula:
#'
#' \deqn{ C(u, v) = \frac{\exp\left( \frac{2 \rho x_u y_v - \rho^2 (x_u^2 +
#' y_v^2)}{2 (1 - \rho^2)} \right)}{\sqrt{1 - \rho^2}} }
#'
#' where \eqn{x_{u}} and \eqn{y_{v}} are the quantile functions (inverse CDF) of
#' the standard normal distribution.
#'
#' @examples
#' # Example usage with rho = 0.5
#' u <- runif(10)
#' v <- runif(10)
#' rho <- 0.5
#' copula_values <- gaussian_copula_function(u, v, rho)
#' plot(copula_values, main = "Gaussian Copula Values", xlab = "Index", ylab = "Copula Value")
#'
#' @seealso \code{\link{qnorm}} for the quantile function of the normal
#'   distribution.
#'
#' @export
#'
gaussian_copula_function <- function(u, v, rho) {
  # rho \in (-1, 1) is correlation
  x_u <- stats::qnorm(u)
  y_v <- stats::qnorm(v)
  c <- exp((2 * rho * x_u * y_v - rho^2 * (x_u^2 + y_v^2)) / (2 * (1 - rho^2))) / sqrt(1 - rho^2)
  return(c)
}

#' Farlie-Gumbel-Morgenstern (FGM) Copula Function
#'
#' This function computes the value of the Farlie-Gumbel-Morgenstern (FGM)
#' copula given two uniform random variables and a correlation parameter. The
#' FGM copula is a simple bivariate copula model used to describe dependence
#' between two random variables.
#'
#' @param u A numeric vector or matrix of uniform random variables (u_i ∈ [0,
#'   1]). Each value in `u` represents a univariate uniform random variable.
#' @param v A numeric vector or matrix of uniform random variables (v_i ∈ [0,
#'   1]). Each value in `v` represents a univariate uniform random variable.
#' @param rho A numeric scalar between -1 and 1 representing the correlation
#'   parameter of the FGM copula.
#'
#' @return A numeric vector or matrix of copula values corresponding to the
#'   input values of `u` and `v`. The return value has the same dimensions as
#'   the inputs `u` and `v`.
#'
#' @details The Farlie-Gumbel-Morgenstern (FGM) copula is given by the formula:
#'
#' \deqn{ C(u, v) = 1 + \rho \cdot (1 - 2u) \cdot (1 - 2v) }
#'
#' where \eqn{\rho} is the correlation parameter and u and v are the uniform
#' random variables. The FGM copula is relatively simple and is used for
#' modeling linear or near-linear dependencies between random variables.
#'
#' @examples
#' # Example usage with rho = 0.3
#' u <- runif(10)
#' v <- runif(10)
#' rho <- 0.3
#' copula_values <- fgm_copula_function(u, v, rho)
#' plot(copula_values, main = "FGM Copula Values", xlab = "Index", ylab = "Copula Value")
#'
#' @export
#'
fgm_copula_function <- function(u, v, rho) {
  # rho \in [-1, 1]
  c <- 1 + rho * (1 - 2 * u) * (1 - 2 * v)
  return(c)
}

#' Copula Function with Vectorization Support
#'
#' This function computes copula values for pairs of uniform random variables
#' using a specified copula type and expansion method. The function supports
#' vectorized operations for efficient computation of copula values for multiple
#' inputs at once. It can handle the 'independent', 'Gaussian', and
#' 'Farlie-Gumbel-Morgenstern (FGM)' copulas. The expansion type determines how
#' the copula is applied to the input arrays.
#'
#' @param U A numeric vector or matrix representing the first set of uniform
#'   random variables. If expansion type is 'U' or 'both', `U` is replicated
#'   along the second dimension.
#' @param V A numeric vector or matrix representing the second set of uniform
#'   random variables. If expansion type is 'V' or 'both', `V` is replicated
#'   along the first dimension.
#' @param expansion_type A character string indicating how to expand the copula
#'   inputs. Can be one of:
#'        - "U": Only the `U` vector is expanded (i.e., `U` is broadcasted to match `V`).
#'        - "V": Only the `V` vector is expanded (i.e., `V` is broadcasted to match `U`).
#'        - "both": Both `U` and `V` are expanded to create a grid of copula evaluations.
#'   Default is "both".
#' @param copula_type A character string specifying the type of copula to
#'   compute. Can be one of:
#'        - "independent": Assumes no dependence between `U` and `V`, returning 1 for all pairs.
#'        - "gaussian": Computes the Gaussian copula, which models correlation between `U` and `V`.
#'        - "fgm": Computes the Farlie-Gumbel-Morgenstern copula, another model for correlation between `U` and `V`.
#'   Default is "independent".
#' @param rho A numeric scalar between -1 and 1, representing the correlation
#'   parameter used in the Gaussian and FGM copulas. Default is 0, indicating no
#'   correlation.
#'
#' @return A numeric matrix or array of copula values. The dimensions of the
#'   return value depend on the `expansion_type`.
#' \itemize{
#'   \item If `expansion_type` is "U", the output will have dimensions (\eqn{n_v}, n).
#'   \item If "V", the output will have dimensions (\eqn{n_u}, n).
#'   \item If "both", the output will have dimensions (\eqn{n_u}, \eqn{n_v}, n).
#' }
#'
#' @details The function computes the copula values for pairs of uniform random
#'   variables `U` and `V`. It supports three types of copulas:
#'
#' \itemize{
#'   \item Independent copula: The copula value is always 1, representing no dependence between `U` and `V`.
#'   \item Gaussian copula: The copula value models correlation between `U` and `V`, with correlation parameter `rho`.
#'   \item FGM copula: The copula value is computed using the Farlie-Gumbel-Morgenstern copula formula with correlation `rho`.
#' }
#'
#'   The function also allows for different methods of expanding `U` and `V` to
#'   create grids of copula evaluations, making it flexible for different data
#'   structures.
#'
#' @examples
#' # Example: Compute Gaussian copula with correlation 0.5
#' U <- matrix(1:6, nrow = 2, ncol = 3)
#' V <- matrix(7:18, nrow = 4, ncol = 3)
#' copula_vals <- copula_function_vectorized(U, V, expansion_type = "both",
#' copula_type = "gaussian", rho = 0.5)
#'
#' # Example: Compute independent copula
#' copula_vals_independent <- copula_function_vectorized(U, V,
#' expansion_type = "both", copula_type = "independent")
#'
#' @export
#'
copula_function_vectorized <- function(U,
                                       V,
                                       expansion_type = "both",
                                       copula_type = "independent",
                                       rho = 0) {
  if (copula_type == "independent") {
    switch(
      expansion_type,
      "U" = return(array(1, dim(V))),
      # U (n)
      # V (n_s, n)
      # return (n_s, n)
      "V" = return(array(1, dim(U))),
      # U (n_s, n)
      # V (n)
      # return (n_s, n)
      "both" = return(array(1, c(dim(
        U
      )[1], dim(
        V
      ))))
      # U (n_u, n)
      # V (n_v, n)
      # return (n_u, n_v, n)
    )
  }
  U[U <= .Machine$double.eps] <- .Machine$double.eps
  U[U >= 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  V[V <= .Machine$double.eps] <- .Machine$double.eps
  V[V >= 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  switch(expansion_type,
         "U" = {
           # U (n)
           # V (n_v, n)
           # return (n_v, n)
           U <- t(matrix(rep(U, dim(V)[1]), ncol = dim(V)[1]))
         },
         "V" = {
           # U (n_u, n)
           # V (n)
           # return (n_u, n)
           V <- t(matrix(rep(V, dim(U)[1]), ncol = dim(U)[1]))
         },
         "both" = {
           # U (n_u, n)
           # V (n_v, n)
           # return (n_u, n_v, n)
           dim(U) <- c(dim(U)[1], 1, dim(U)[2])
           U <- U[, rep(1, dim(V)[1]), ]
           dim(V) <- c(1, dim(V)[1], dim(V)[2])
           V <- V[rep(1, dim(U)[1]), , ]
         })
  switch(copula_type, "gaussian" = return(gaussian_copula_function(U, V, rho)), "fgm" = return(fgm_copula_function(U, V, rho)))
}


#' Gradient of the Gaussian Copula
#'
#' This function computes the gradients of the Gaussian copula with respect to
#' both components `u` and `v`. The gradients represent how sensitive the copula
#' values are to changes in `u` and `v`. These gradients are useful in various
#' optimization and modeling tasks, such as in maximum likelihood estimation for
#' copula-based models.
#'
#' @param u A numeric vector of uniform random variables representing the first
#'   set of inputs. The values must be in the range (0, 1).
#' @param v A numeric vector of uniform random variables representing the second
#'   set of inputs. The values must be in the range (0, 1).
#' @param rho A numeric scalar between -1 and 1, representing the correlation
#'   parameter of the Gaussian copula.
#'
#' @return A list containing two components:
#' \itemize{
#'   \item \code{cu}: A numeric vector of the gradient of the Gaussian copula with respect to `u`.
#'   \item \code{cv}: A numeric vector of the gradient of the Gaussian copula with respect to `v`.
#' }
#'
#' @details The gradient of the Gaussian copula is calculated with respect to
#' both `u` and `v`. For each element of `u` and `v`, the function computes the
#' partial derivatives of the Gaussian copula with respect to the first and
#' second inputs, respectively. This is useful when working with the Gaussian
#' copula in optimization algorithms or when performing sensitivity analysis.
#'
#' The gradients are computed as follows: \deqn{c_u = \frac{\partial
#' C(u,v)}{\partial u} = \frac{C(u,v) \cdot \left( \rho y_v - \rho^2 x_u
#' \right)}{(1-\rho^2) dnorm(x_u)}} \deqn{c_v = \frac{\partial C(u,v)}{\partial
#' v} = \frac{C(u,v) \cdot \left( \rho x_u - \rho^2 y_v \right)}{(1-\rho^2)
#' dnorm(y_v)}} where \eqn{C(u, v)} is the Gaussian copula, and \eqn{x_{u}} and
#' \eqn{y_{v}} are the quantile transformations of `u` and `v`, respectively.
#'
#' @examples
#' # Example: Compute the gradient of the Gaussian copula with correlation 0.5
#' u <- runif(10)
#' v <- runif(10)
#' rho <- 0.5
#' gradients <- gaussian_copula_gradient(u, v, rho)
#'
#' # Display the gradient components
#' gradients$cu
#' gradients$cv
#'
#' @seealso \code{\link{gaussian_copula_function}} for computing the Gaussian
#'   copula.
#'
#' @export
#'
gaussian_copula_gradient <- function(u, v, rho) {
  # rho \in (-1, 1) is correlation
  x_u <- stats::qnorm(u)
  y_v <- stats::qnorm(v)
  c_u <- gaussian_copula_function(u, v, rho) * (rho * y_v - rho^2 * x_u) / (1 - rho^2) / stats::dnorm(x_u)
  c_v <- gaussian_copula_function(u, v, rho) * (rho * x_u - rho^2 * y_v) / (1 - rho^2) / stats::dnorm(y_v)
  return(list(cu = c_u, cv = c_v))
}

#' Gradient of the FGM Copula
#'
#' This function computes the gradients of the FGM (Farlie-Gumbel-Morgenstern)
#' copula with respect to both components `u` and `v`. The gradients represent
#' how sensitive the copula values are to changes in `u` and `v`. These
#' gradients are useful in various optimization tasks, such as in maximum
#' likelihood estimation for copula-based models or other applications involving
#' copula gradients.
#'
#' @param u A numeric vector of uniform random variables representing the first
#'   set of inputs. The values must be in the range (0, 1).
#' @param v A numeric vector of uniform random variables representing the second
#'   set of inputs. The values must be in the range (0, 1).
#' @param rho A numeric scalar between -1 and 1, representing the correlation
#'   parameter of the FGM copula.
#'
#' @return A list containing two components:
#' \itemize{
#'   \item \code{cu}: A numeric vector of the gradient of the FGM copula with respect to `u`.
#'   \item \code{cv}: A numeric vector of the gradient of the FGM copula with respect to `v`.
#' }
#'
#' @details The gradient of the FGM copula is calculated with respect to both
#' `u` and `v`. For each element of `u` and `v`, the function computes the
#' partial derivatives of the FGM copula with respect to the first and second
#' inputs, respectively. This is useful when working with the FGM copula in
#' optimization algorithms or sensitivity analysis.
#'
#' The gradients are computed as follows: \deqn{c_u = \frac{\partial
#' C(u,v)}{\partial u} = -2 \rho (1 - 2v)} \deqn{c_v = \frac{\partial
#' C(u,v)}{\partial v} = -2 \rho (1 - 2u)} where \eqn{C(u, v)} is the FGM
#' copula, and \eqn{\rho} is the correlation parameter of the copula.
#'
#' @examples
#' # Example: Compute the gradient of the FGM copula with correlation 0.5
#' u <- runif(10)
#' v <- runif(10)
#' rho <- 0.5
#' gradients <- fgm_copula_gradient(u, v, rho)
#'
#' # Display the gradient components
#' gradients$cu
#' gradients$cv
#'
#' @seealso \code{\link{fgm_copula_function}} for computing the FGM copula.
#'
#' @export
#'
fgm_copula_gradient <- function(u, v, rho) {
  # rho \in [-1, 1]
  c_u <- -2 * rho * (1 - 2 * v)
  c_v <- -2 * rho * (1 - 2 * u)
  return(list(cu = c_u, cv = c_v))
}

#' Vectorized Gradient of a Copula
#'
#' This function computes the gradients of a specified copula (either Gaussian
#' or FGM) with respect to both components `U` and `V` in a vectorized manner.
#' It efficiently computes the gradients for multiple pairs of variables `U` and
#' `V` simultaneously. The gradients represent the partial derivatives of the
#' copula with respect to its inputs and are useful in optimization tasks, such
#' as maximum likelihood estimation or sensitivity analysis.
#'
#' @param U A numeric matrix of uniform random variables representing the first
#'   set of inputs. The values must be in the range (0, 1).
#' @param V A numeric matrix of uniform random variables representing the second
#'   set of inputs. The values must be in the range (0, 1).
#' @param copula_type A character string specifying the copula type to use.
#'   Possible values are:
#'   - "independent": For the independent copula (no dependency between `U` and `V`).
#'   - "gaussian": For the Gaussian copula.
#'   - "fgm": For the FGM (Farlie-Gumbel-Morgenstern) copula.
#'   Default is `'independent'`.
#' @param rho A numeric scalar between -1 and 1, representing the correlation
#'   parameter for copulas that require it (i.e., Gaussian and FGM). Default is
#'   0.
#'
#' @return A list containing two components:
#' \itemize{
#'   \item \code{cu}: A numeric array representing the gradient of the copula with respect to `U`.
#'   \item \code{cv}: A numeric array representing the gradient of the copula with respect to `V`.
#' }
#'   The arrays will have dimensions corresponding to the number of elements in
#'   `U` and `V`.
#'
#' @details This function computes the gradients of different copulas
#' (independent, Gaussian, and FGM) in a vectorized fashion for an array of
#' inputs `U` and `V`. It efficiently handles the computation for multiple
#' variable pairs simultaneously, which can be crucial for high-dimensional
#' tasks.
#'
#' For the independent copula, the gradients are zero as there is no dependence
#' between `U` and `V`. For the Gaussian and FGM copulas, the gradients are
#' computed using their respective formulas (via the `gaussian_copula_gradient`
#' and `fgm_copula_gradient` functions). The function applies bounds to ensure
#' that the values in `U` and `V` lie within the valid range (0, 1) before
#' performing the gradient calculation.
#'
#' @examples
#' # Example: Compute the gradient of the Gaussian copula with correlation 0.5
#' U <- matrix(runif(10), nrow = 5, ncol = 2)
#' V <- matrix(runif(10), nrow = 5, ncol = 2)
#' rho <- 0.5
#' gradients <- copula_gradient_vectorized(U, V, copula_type = "gaussian", rho = rho)
#'
#' # Display the gradient components for the first pair of inputs
#' gradients$cu[1, , 1]
#' gradients$cv[1, , 1]
#'
#' @seealso \code{\link{gaussian_copula_gradient}},
#'   \code{\link{fgm_copula_gradient}} for gradient calculations of individual
#'   copulas.
#'
#' @export
#'
copula_gradient_vectorized <- function(U,
                                       V,
                                       copula_type = "independent",
                                       rho = 0) {
  # U (n_u, n)
  # V (n_v, n)
  if (copula_type == "independent") {
    return(list(cu = array(0, c(
      dim(U)[1], dim(V)[1], dim(U)[2]
    )), cv = array(0, c(
      dim(U)[1], dim(V)[1], dim(U)[2]
    ))))
  }
  U[U <= .Machine$double.eps] <- .Machine$double.eps
  U[U >= 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  V[V <= .Machine$double.eps] <- .Machine$double.eps
  V[V >= 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  # Expoand to match dimensions
  dim(U) <- c(dim(U)[1], 1, dim(U)[2])
  U <- U[, rep(1, dim(V)[1]), ]
  dim(V) <- c(1, dim(V)[1], dim(V)[2])
  V <- V[rep(1, dim(U)[1]), , ]
  switch(copula_type, "gaussian" = return(gaussian_copula_gradient(U, V, rho)), "fgm" = return(gaussian_copula_gradient(U, V, rho)))
}
