#' Vectorized Weighting Function
#'
#' This function generates a matrix of weights where each entry is set to 1,
#' representing a simple weighting scheme. It computes the matrix in a
#' vectorized fashion, efficiently handling multiple input sizes. The function
#' is primarily useful when a constant weight (1) is required between two sets
#' of values, such as in certain copula-related calculations or in exploratory
#' data analysis.
#'
#' @param s1 A numeric vector representing the first set of input values (e.g.,
#'   one dimension of the data).
#' @param s0 A numeric vector representing the second set of input values (e.g.,
#'   another dimension of the data).
#'
#' @return A numeric array with dimensions `length(s1) x length(s0)`. All values
#'   in the array are set to 1.
#'
#' @examples
#' # Example: Generate a weighting matrix for two input vectors
#' s1 <- 1:3
#' s0 <- c(0.1, 0.5, 0.9)
#' W <- identity_weighting_function_vectorized(s1, s0)
#'
#' # Display the resulting weight matrix
#' print(W)
#'
#' @export
#'
identity_weighting_function_vectorized <- function(s1, s0) {
  # return (n_s1, n_s0)
  return(array(1, c(length(s1), length(s0))))
}

#' Vectorized G Function
#'
#' This function constructs a three-dimensional array by combining three
#' matrices. The resulting array can be used in various multi-dimensional data
#' processing or for computing interactions between two sets of input vectors.
#' The function utilizes vectorized operations for efficiency, making it
#' suitable for large datasets.
#'
#' @param s1 A numeric vector of length `n_s1` representing the first set of
#'   input values.
#' @param s0 A numeric vector of length `n_s0` representing the second set of
#'   input values.
#'
#' @return A three-dimensional numeric array of dimensions `n_s1 x n_s0 x 3`,
#'   containing:
#'   \itemize{
#'     \item The first slice of the array, corresponding to `s1` repeated for each element of `s0`.
#'     \item The second slice of the array, corresponding to `s0` repeated for each element of `s1`.
#'     \item A third slice filled with ones, representing a constant dimension.
#'   }
#'
#' @details The function generates a three-dimensional array by combining three
#' components: 1. A matrix where each element of `s1` is repeated for each value
#' in `s0`. 2. A matrix where each element of `s0` is repeated for each value in
#' `s1`. 3. A matrix filled with ones, serving as a constant dimension.
#'
#' This operation can be useful in scenarios where multi-dimensional
#' interactions or transformations of two sets of data are required, such as in
#' copula or statistical modeling.
#'
#' @examples
#' # Example: Create a 3D array combining two vectors
#' s1 <- 1:3
#' s0 <- c(0.1, 0.5, 0.9)
#' result <- identity_g_function_vectorized(s1, s0)
#'
#' # Display the resulting 3D array
#' print(result)
#'
#' @export
#'
identity_g_function_vectorized <- function(s1, s0) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, dim_g)
  return(abind(replicate(length(s0), s1), t(replicate(length(
    s1
  ), s0)), array(1, c(
    length(s1), length(s0)
  )), along = 3))
}

#' Principal Score Prediction (Vectorized)
#'
#' This function performs vectorized prediction of principal scores based on a
#' given bandwidth object and predictor variables. It uses the kernel density
#' estimation function `npcdens` from the `np` package to predict the target
#' variable given a set of predictors.
#'
#' @param bw A bandwidth object generated from the `npcdensbw` function. This
#'   contains the bandwidth parameter for kernel density estimation.
#' @param S A numeric vector of length \eqn{n_s}, containing the target variable
#'   values for which predictions will be made.
#' @param X A numeric matrix with dimensions (n, dims(X)), or data frame
#'   containing the predictor variables for which principal score predictions
#'   are to be made.
#'
#' @return A numeric matrix of predicted values with dimensions (\eqn{n_s},  n),
#'   where each row corresponds to predictions for the respective values in S.
#'
#' @details This function leverages the `npcdens` function from the `np` package
#' to perform kernel density estimation using the bandwidth object provided in
#' `bw`. The target variable `S` is predicted based on the predictors `X` using
#' the nonparametric kernel density estimate, where each row of the resulting
#' output corresponds to a predicted value for a corresponding target value in
#' `S`. The computation is vectorized for efficiency.
#'
#' @examples
#' # Example: Compute principal score predictions
#' set.seed(42)
#' X <- matrix(rnorm(200), ncol=2)  # Predictor variables
#' S <- rnorm(100)  # Target variable
#'
#' # Compute bandwidth using the npcdensbw function
#' bw <- np::npcdensbw(xdat=X, ydat=S)
#'
#' # Compute predicted values for principal scores
#' predictions <- principal_score_predict_vectorized(bw, S, X)
#' print(predictions)
#'
#' @seealso \code{\link[np]{npcdens}} for kernel density estimation and
#'   \code{\link[np]{npcdensbw}} for bandwidth selection.
#'
#' @importFrom np npcdens
#'
#' @export
#'
principal_score_predict_vectorized <- function(bw, S, X) {
  # S (n_s, )
  # X (n, dim_x)
  # return (n_s, n)
  return(t(array(
    stats::fitted(npcdens(
      bws = bw,
      exdat = do.call(rbind, replicate(length(S), X, simplify = FALSE)),
      eydat = cbind(as.vector(t(
        replicate(dim(X)[1], S)
      )))
    )), c(dim(X)[1], length(S))
  )))
}

#' Vectorized Function for Linear Model Prediction
#'
#' This function performs vectorized predictions of the target variable based on
#' a given linear model (`lm_model`) and predictor variables. It applies the
#' model to a set of predictors and target values, returning the predictions in
#' a matrix format.
#'
#' @param lm_model A linear model object created by \code{\link{lm}}. This model
#'   is used to generate predictions for the target variable.
#' @param S A numeric vector of length \eqn{n_s}, containing the S predictor
#'   variables (independent variables) for which predictions will be made.
#' @param X A numeric matrix with dimensions (n, dims(X)) or data frame
#'   containing the X predictor variables (independent variables) for which the
#'   target variable predictions will be made.
#'
#' @return A numeric matrix of predicted values with dimensions (\eqn{n_s}, n),
#'   where each row corresponds to a predicted value for a respective target
#'   value in S and X.
#'
#' @details This function takes a linear model (`lm_model`), a predictor
#' variables S and X. It then generates predictions using the provided linear
#' model. The prediction is vectorized, meaning that it can handle multiple
#' predictor variables S and X efficiently.
#'
#' @examples
#' # Example:
#' set.seed(42)
#' X <- matrix(rnorm(100), ncol=2)  # Predictor variables
#' S <- rnorm(50)  # Predictor variables
#' Y <- rnorm(50) # Target variables
#' S_X <- data.frame(cbind(S, data.frame(X)))
#' colnames(S_X) <- c('S',  paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
#' lm_model <- lm(Y ~., data=S_X)
#' S_new = rnorm(100)
#' # Compute predicted values for the target variable
#' predictions <- lm_prediction_vectorized(lm_model, S_new, X)
#' print(dim(predictions))
#'
#' @seealso \code{\link{lm}} for linear model fitting.
#'
#' @export
#'
lm_prediction_vectorized <- function(lm_model, S, X) {
  # return (n_s, n)
  newdata = data.frame(cbind(as.vector(t(
    replicate(dim(X)[1], S)
  )), do.call(
    rbind, replicate(length(S), X, simplify = FALSE)
  )))
  colnames(newdata) = c('S', paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
  return(t(array(
    stats::predict.lm(lm_model, newdata = newdata),
    c(dim(X)[1], length(S))
  )))
}

#' Vectorized \eqn{W * G * G^T} Computation
#'
#' This function computes a vectorized version of a weighting and generalized
#' transformation between two sets of values `s1` and `s0`, using the specified
#' weighting function and `g` transformation function.
#'
#' @param s1 A numeric vector of size \eqn{n_{s1}} representing the first set of
#'   values.
#' @param s0 A numeric vector of size \eqn{n_{s0}} representing the second set
#'   of values.
#' @param weighting_function_vectorized A function that computes the weights
#'   based on `s1` and `s0`. The default is
#'   \code{identity_weighting_function_vectorized}.
#' @param g_function_vectorized A function that computes the generalized group
#'   transformation based on `s1` and `s0`. The default is
#'   \code{identity_g_function_vectorized}.
#'
#' @return A numeric array of size (\eqn{n_{s1}, n_{s0}, dim_g, dim_g}), where
#'   \eqn{dim_g} corresponds to the third dimension of the output from the
#'   `g_function_vectorized`.
#'
#' @details The function performs a series of vectorized operations between the
#' two input vectors `s1` and `s0`, applying the specified weighting function
#' (`weighting_function_vectorized`) and the generalized transformation function
#' (`g_function_vectorized`) to these values. The results are then combined in a
#' multi-dimensional array.
#'
#' @examples
#' # Example: Compute vectorized weighted GGGT
#' s1 <- rnorm(10)  # First set of values
#' s0 <- rnorm(5)   # Second set of values
#' result <- wggt_vectorized(s1, s0)
#' print(dim(result))
#'
#' @export
#'
wggt_vectorized <- function(s1,
                            s0,
                            weighting_function_vectorized = identity_weighting_function_vectorized,
                            g_function_vectorized = identity_g_function_vectorized) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, dim_g, dim_g)
  # to confirm
  w = weighting_function_vectorized(s1, s0)
  g = g_function_vectorized(s1, s0)
  dim(w) = c(dim(w), 1, 1)
  w = w[, , rep(1, dim(g)[3]), rep(1, dim(g)[3])]
  g1 = g
  g2 = g
  dim(g1) = c(dim(g1), 1)
  dim(g2) = c(dim(g1)[1], dim(g1)[2], 1, dim(g1)[3])
  g1 = g1[, , , rep(1, dim(g)[3])]
  g2 = g2[, , rep(1, dim(g)[3]), ]
  return(w * g1 * g2)
}

#' Vectorized \eqn{W * G} Computation
#'
#' This function computes a vectorized version of the weighting and
#' transformation between two sets of values `s1` and `s0`, using the specified
#' weighting function and transformation function.
#'
#' @param s1 A numeric vector of size \eqn{n_{s1}} representing the first set of
#'   values.
#' @param s0 A numeric vector of size \eqn{n_{s0}} representing the second set
#'   of values.
#' @param weighting_function_vectorized A function that computes the weights
#'   based on `s1` and `s0`. The default is
#'   \code{identity_weighting_function_vectorized}.
#' @param g_function_vectorized A function that computes the generalized
#'   transformation based on `s1` and `s0`. The default is
#'   \code{identity_g_function_vectorized}.
#'
#' @return A numeric array of size (\eqn{n_{s1}, n_{s0}, dim_g}), where
#'   \eqn{dim_g} corresponds to the third dimension of the output from the
#'   `g_function_vectorized`.
#'
#' @details The function performs a series of vectorized operations between the
#' two input vectors `s1` and `s0`, applying the specified weighting function
#' (`weighting_function_vectorized`) and the transformation function
#' (`g_function_vectorized`) to these values. The results are then combined in a
#' multi-dimensional array.
#'
#' @examples
#' # Example: Compute vectorized weighted transformation
#' s1 <- rnorm(10)  # First set of values (e.g., target variable)
#' s0 <- rnorm(5)   # Second set of values (e.g., reference variable)
#' result <- wg_vectorized(s1, s0)
#' print(dim(result))
#'
#' @export
#'
wg_vectorized <- function(s1,
                          s0,
                          weighting_function_vectorized = identity_weighting_function_vectorized,
                          g_function_vectorized = identity_g_function_vectorized) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, dim_g)
  w = weighting_function_vectorized(s1, s0)
  g = g_function_vectorized(s1, s0)
  dim(w) = c(dim(w), 1)
  w = w[, , rep(1, dim(g)[3])]
  return(w * g)
}
