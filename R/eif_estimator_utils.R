#' Conditional Integral Calculation for Copula Models
#'
#' This function computes a conditional integral based on a copula model,
#' typically used in financial or econometric models involving copula-based
#' dependencies.
#'
#' @param s1 A numeric vector of size \eqn{n_{s1}} representing the first set of
#'   values.
#' @param s0 A numeric vector of size \eqn{n_{s0}} representing the second set
#'   of values.
#' @param psp_s1 A matrix of size (\eqn{n_{s1}}, n), representing the
#'   probabilities associated with the first set of values s1.
#' @param psp_s0 A matrix of size (\eqn{n_{s0}}, n), representing the
#'   probabilities associated with the second set of values s0.
#' @param psp_s1_s0 A numeric matrix of size (\eqn{n_{s1}}, \eqn{n_{s0}}, n)
#'   representing joint probabilities for s1 and s0 .
#' @param psp_s1_cdf A matrix of size (\eqn{n_{s1}}, n) representing the
#'   cumulative distribution function (CDF) of s1.
#' @param psp_s0_cdf A matrix of size (\eqn{n_{s0}}, n) representing the
#'   cumulative distribution function (CDF) of s0.
#' @param Z A numeric vector of size n  representing the first condition (e.g.,
#'   treatment assignment indicator).
#' @param Tp A numeric value representing the total probability for the positive
#'   outcome (e.g., treatment or event probability).
#' @param S A numeric vector of size n , typically representing time or a
#'   sequential index (used in conjunction with `Z`).
#' @param copula_type A character string indicating the type of copula to be
#'   used. Defaults to "independent". Other options are "gaussian" and "fgm".
#' @param rho A numeric value between -1 and 1 representing the correlation
#'   parameter for copula models. Default is 0 (for independence).
#'
#' @return A numeric array of size (\eqn{n_{s1}}, \eqn{n_{s0}}, n), representing
#'   the weighted and transformed probabilities of the copula model after
#'   applying the integral calculations.
#'
#' @details The function computes the conditional integral for a copula-based
#' model with specified dependencies between `s1` and `s0`. It incorporates the
#' given probability values for both `s1` and `s0`, their CDFs, and joint
#' probabilities. It also applies the specified copula (independent, gaussian,
#' or fgm) to model dependencies between the variables. The integral is
#' calculated through a series of matrix operations, utilizing the copula
#' function, its gradient, and expansions for handling interactions between the
#' variables.
#'
#' Specifically:
#' - The function calculates the copula weights (`l1_B`) and the copula gradient.
#' - It computes the contributions from `s1` and `s0` (`part2_s1` and `part2_s0`) based on these copula weights and gradients.
#' - The final result is a product of the copula probability matrix (`psp_s1_s0`) and the computed values (`lp`), accounting for both conditions `Z` and `Tp`.
#'
#' @examples
#' # Example: Compute the conditional integral with the independent copula
#' s1 <- rnorm(20)  # First set of values (e.g., target variable)
#' s0 <- rnorm(5)   # Second set of values (e.g., reference variable)
#' psp_s1 <- matrix(runif(200), 20, 10)  # Probability values for s1
#' psp_s0 <- matrix(runif(50), 5, 10)   # Probability values for s0
#' psp_s1_s0 <- matrix_multiply_with_expansion(psp_s1, psp_s0, 'both')  # Joint probability matrix
#' psp_s1_cdf <- trapz_cdf(s1, psp_s1)  # CDF of s1
#' psp_s0_cdf <- trapz_cdf(s0, psp_s0)  # CDF of s0
#' Z <- sample(0:1, 10, replace = TRUE)  # Treatment assignment indicator
#' Tp <- 0.7  # Event probability
#' S <- 1:10  # Sequential index
#'
#' # Call the int2 function
#' result <- int2(s1, s0, psp_s1, psp_s0, psp_s1_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S)
#' print(dim(result))
#'
#' @export
#'
int2 <- function(s1,
                 s0,
                 psp_s1,
                 psp_s0,
                 psp_s1_s0,
                 psp_s1_cdf,
                 psp_s0_cdf,
                 Z,
                 Tp,
                 S,
                 copula_type = 'independent',
                 rho = 0) {
  # return(n_s1, n_s0, n)
  l1_B = copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'both', copula_type, rho)
  copula_gradient = copula_gradient_vectorized(psp_s1_cdf, psp_s0_cdf, copula_type, rho)
  part2_s1 = l1_B + matrix_multiply_with_expansion(psp_s1_cdf - 1 * outer(s1, S, FUN =
                                                                            ">="), copula_gradient$cu, 'U')
  # part2_s1 (n_s1, n_s0, n)
  part2_s0 = l1_B + matrix_multiply_with_expansion(copula_gradient$cv, psp_s0_cdf - 1 * outer(s0, S, FUN =
                                                                                                ">="), 'V')
  # part2_s0 (n_s1, n_s0, n)
  lp = l1_B - sweep(part2_s1, MARGIN = 3, Z / Tp, `*`) - sweep(part2_s0, MARGIN =
                                                                 3, (1 - Z) / (1 - Tp), `*`)
  return(lp * psp_s1_s0)
}

#' Compute the Integral over s0
#'
#' This function computes a component of an integral related to the copula-based
#' modeling framework. It involves calculating a weighted sum over s0, using the
#' given copula and correlation parameters. The result is used in the context of
#' a model that involves score-based functions or copula-based likelihood
#' estimation.
#'
#' @param s1 A numeric vector of size \eqn{n_{s1}} representing the first set of
#'   values.
#' @param psp_s0 A matrix of size (\eqn{n_{s0}}, n), representing the
#'   probabilities associated with the second set of values s0.
#' @param psp_s1_cdf A matrix of size (\eqn{n_{s1}}, n) representing the
#'   cumulative distribution function (CDF) of s1.
#' @param psp_s0_cdf A matrix of size (\eqn{n_{s0}}, n) representing the
#'   cumulative distribution function (CDF) of s0.
#' @param Z A numeric vector of size n  representing the first condition (e.g.,
#'   treatment assignment indicator).
#' @param Tp A numeric value representing the total probability for the positive
#'   outcome (e.g., treatment or event probability).
#' @param S A numeric vector of size n , typically representing time or a
#'   sequential index (used in conjunction with `Z`).
#' @param copula_type A character string indicating the type of copula to be
#'   used. Defaults to "independent". Other options are "gaussian" and "fgm".
#' @param rho A numeric value between -1 and 1 representing the correlation
#'   parameter for copula models. Default is 0 (for independence).
#'
#' @return A matrix of size (\eqn{n_{s0}}, n).
#'
#' @examples
#' # Example usage of the function
#' s1 <- rnorm(20)  # First set of values (e.g., target variable)
#' s0 <- rnorm(5)   # Second set of values (e.g., reference variable)
#' psp_s1 <- matrix(runif(100), 20, 10)   # Probability values for s1
#' psp_s0 <- matrix(runif(50), 5, 10)   # Probability values for s0
#' psp_s1_cdf <- trapz_cdf(s1, psp_s1)  # CDF of s1
#' psp_s0_cdf <- trapz_cdf(s0, psp_s0)  # CDF of s0
#' Z <- sample(0:1, 10, replace = TRUE)  # Treatment assignment indicator
#' Tp <- 0.7  # Event probability
#' S <- 1:10  # Sequential index
#'
#' # Call the function
#' result <- int_s0(s1, psp_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S)
#' print(dim(result))
#'
#' @export
#'
int_s0 <- function(s1,
                   psp_s0,
                   psp_s1_cdf,
                   psp_s0_cdf,
                   Z,
                   Tp,
                   S,
                   copula_type = 'independent',
                   rho = 0) {
  # return (n_s0, n)
  psp_s1_cdf_S = psp_s1_cdf[cbind(sapply(S, function(s) {
    max(which(s1 <= s))
  }), 1:length(S))]
  # psp_s1_cdf_S (n)
  part2 = sweep(
    copula_function_vectorized(psp_s1_cdf_S, psp_s0_cdf, 'U', copula_type, rho),
    MARGIN = 2,
    Z / Tp,
    `*`
  )
  return(part2 * psp_s0)
}

#' Compute the Integral over s0
#'
#' This function computes a component of an integral related to the copula-based
#' modeling framework. It involves calculating a weighted sum over s0, using the
#' given copula and correlation parameters. The result is used in the context of
#' a model that involves score-based functions or copula-based likelihood
#' estimation.
#'
#' @param s0 A numeric vector of size \eqn{n_{s0}} representing the first set of
#'   values.
#' @param psp_s1 A matrix of size (\eqn{n_{s1}}, n), representing the
#'   probabilities associated with the second set of values s1.
#' @param psp_s1_cdf A matrix of size (\eqn{n_{s1}}, n) representing the
#'   cumulative distribution function (CDF) of s1.
#' @param psp_s0_cdf A matrix of size (\eqn{n_{s0}}, n) representing the
#'   cumulative distribution function (CDF) of s0.
#' @param Z A numeric vector of size n  representing the first condition (e.g.,
#'   treatment assignment indicator).
#' @param Tp A numeric value representing the total probability for the positive
#'   outcome (e.g., treatment or event probability).
#' @param S A numeric vector of size n , typically representing time or a
#'   sequential index (used in conjunction with `Z`).
#' @param copula_type A character string indicating the type of copula to be
#'   used. Defaults to "independent". Other options are "gaussian" and "fgm".
#' @param rho A numeric value between -1 and 1 representing the correlation
#'   parameter for copula models. Default is 0 (for independence).
#'
#' @return A matrix of size (\eqn{n_{s1}}, n).
#'
#' @examples
#' # Example usage of the function
#' s1 <- rnorm(20)  # First set of values (e.g., target variable)
#' s0 <- rnorm(5)   # Second set of values (e.g., reference variable)
#' psp_s1 <- matrix(runif(100), 20, 10)   # Probability values for s1
#' psp_s0 <- matrix(runif(50), 5, 10)   # Probability values for s0
#' psp_s1_cdf <- trapz_cdf(s1, psp_s1)  # CDF of s1
#' psp_s0_cdf <- trapz_cdf(s0, psp_s0)  # CDF of s0
#' Z <- sample(0:1, 10, replace = TRUE)  # Treatment assignment indicator
#' Tp <- 0.7  # Event probability
#' S <- 1:10  # Sequential index
#'
#' # Call the function
#' result <- int_s1(s0, psp_s1, psp_s1_cdf, psp_s0_cdf, Z, Tp, S)
#' print(dim(result))
#'
#' @export
#'
int_s1 <- function(s0,
                   psp_s1,
                   psp_s1_cdf,
                   psp_s0_cdf,
                   Z,
                   Tp,
                   S,
                   copula_type = 'independent',
                   rho = 0) {
  # return (n_s1, n)
  psp_s0_cdf_S = psp_s0_cdf[cbind(sapply(S, function(s) {
    max(which(s0 <= s))
  }), 1:length(S))]
  # psp_s0_cdf_S (n)
  part2 = sweep(
    copula_function_vectorized(psp_s1_cdf, psp_s0_cdf_S, 'V', copula_type, rho),
    MARGIN = 2,
    (1 - Z) / (1 - Tp),
    `*`
  )
  return(part2 * psp_s1)
}

B_int2_function <- function(s1,
                            s0,
                            psp_s1,
                            psp_s0,
                            psp_s1_s0,
                            psp_s1_cdf,
                            psp_s0_cdf,
                            Z,
                            Tp,
                            S,
                            copula_type = 'independent',
                            rho = 0,
                            weighting_function_vectorized = identity_weighting_function_vectorized,
                            g_function_vectorized = identity_g_function_vectorized) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # X (n, dim_x)
  # return (n_s1, n_s0, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1,
                         s0,
                         weighting_function_vectorized,
                         g_function_vectorized)
  # wggt (n_s1, n_s0, dim_g, dim_g)
  lp = int2(s1,
            s0,
            psp_s1,
            psp_s0,
            psp_s1_s0,
            psp_s1_cdf,
            psp_s0_cdf,
            Z,
            Tp,
            S,
            copula_type,
            rho)
  dim(wggt) = c(dim(wggt)[1], dim(wggt)[2], 1, dim(wggt)[3], dim(wggt)[4])
  wggt = wggt[, , rep(1, length(Z)), , ]
  dim(lp) = c(dim(lp), 1, 1)
  lp = lp[, , , rep(1, dim(wggt)[4]), rep(1, dim(wggt)[5])]
  return(wggt * lp)
}

B_ints0_function <- function(s1,
                             s0,
                             psp_s0,
                             psp_s1_cdf,
                             psp_s0_cdf,
                             Z,
                             Tp,
                             S,
                             copula_type = 'independent',
                             rho = 0,
                             weighting_function_vectorized = identity_weighting_function_vectorized,
                             g_function_vectorized = identity_g_function_vectorized) {
  # return (n_s0, n, dim_g, dim_g)
  wggt = aperm(
    wggt_vectorized(S, s0, weighting_function_vectorized, g_function_vectorized),
    c(2, 1, 3, 4)
  )
  # wggt (n_s0, n, dim_g, dim_g)
  part2 = int_s0(s1, psp_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  dim(part2) = c(dim(part2), 1, 1)
  part2 = part2[, , rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * part2)
}

B_ints1_function <- function(s1,
                             s0,
                             psp_s1,
                             psp_s1_cdf,
                             psp_s0_cdf,
                             Z,
                             Tp,
                             S,
                             copula_type = 'independent',
                             rho = 0,
                             weighting_function_vectorized = identity_weighting_function_vectorized,
                             g_function_vectorized = identity_g_function_vectorized) {
  # return (n_s1, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1, S, weighting_function_vectorized, g_function_vectorized)
  # wggt (n_s1, n, dim_g, dim_g)
  part2 = int_s1(s0, psp_s1, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  dim(part2) = c(dim(part2), 1, 1)
  part2 = part2[, , rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * part2)
}


vec_int2_function <- function(mu1_lm,
                              mu0_lm,
                              s1,
                              s0,
                              psp_s1,
                              psp_s0,
                              psp_s1_s0,
                              psp_s1_cdf,
                              psp_s0_cdf,
                              Z,
                              Tp,
                              X,
                              S,
                              copula_type = 'independent',
                              rho = 0,
                              weighting_function_vectorized = identity_weighting_function_vectorized,
                              g_function_vectorized = identity_g_function_vectorized) {
  mu1 = lm_prediction_vectorized(mu1_lm, s1, X)
  # mu1 (n_s1, n)
  mu0 = lm_prediction_vectorized(mu0_lm, s0, X)
  # mu1 (n_s0, n)
  wg = wg_vectorized(s1,
                     s0,
                     weighting_function_vectorized,
                     g_function_vectorized)
  # wg (n_s1, n_s0, dim_g)
  lp = int2(s1,
            s0,
            psp_s1,
            psp_s0,
            psp_s1_s0,
            psp_s1_cdf,
            psp_s0_cdf,
            Z,
            Tp,
            S,
            copula_type,
            rho)
  dim(wg) = c(dim(wg)[1], dim(wg)[2], 1, dim(wg)[3])
  wg = wg[, , rep(1, length(Z)), ]
  dim(lp) = c(dim(lp), 1)
  lp = lp[, , , rep(1, dim(wg)[4])]
  # lp (n_s1, n_s0, n)
  wg_lp  = wg * lp
  # wglp (n_s1, n_s0, n, dim_g)
  dim(mu1) = c(dim(mu1)[1], 1, dim(mu1)[2], 1)
  dim(mu0) = c(1, dim(mu0)[1], dim(mu0)[2], 1)
  mu1 = mu1[, rep(1, length(s0)), , rep(1, dim(wg_lp)[4])]
  mu0 = mu0[rep(1, length(s1)), , , rep(1, dim(wg_lp)[4])]
  return(list(eta1 = wg_lp * mu1, eta0 = wg_lp * mu0))
}


vec_ints0_function <- function(mu0_lm,
                               s1,
                               s0,
                               psp_s0,
                               psp_s1_cdf,
                               psp_s0_cdf,
                               Z,
                               Tp,
                               X,
                               S,
                               Y,
                               copula_type = 'independent',
                               rho = 0,
                               weighting_function_vectorized = identity_weighting_function_vectorized,
                               g_function_vectorized = identity_g_function_vectorized) {
  wg = aperm(
    wg_vectorized(S, s0, weighting_function_vectorized, g_function_vectorized),
    c(2, 1, 3)
  )
  # wg (n_s0, n_s1 = n, dim_g)
  l = int_s0(s1, psp_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # l(n_s0, n)
  # for eta1
  l_eta1 = sweep(l, MARGIN = 2, Y, `*`)
  dim(l_eta1) = c(dim(l_eta1), 1)
  l_eta1 = l_eta1[, , rep(1, dim(wg)[3])]
  # for eta0
  mu0 = lm_prediction_vectorized(mu0_lm, s0, X)
  # mu0 (n_s0, n)
  l_eta0 = l * mu0
  dim(l_eta0) = c(dim(l_eta0), 1)
  l_eta0 = l_eta0[, , rep(1, dim(wg)[3])]
  return(list(eta1 = wg * l_eta1, eta0 = wg * l_eta0))
}

vec_ints1_function <- function(mu1_lm,
                               s1,
                               s0,
                               psp_s1,
                               psp_s1_cdf,
                               psp_s0_cdf,
                               Z,
                               Tp,
                               X,
                               S,
                               Y,
                               copula_type = 'independent',
                               rho = 0,
                               weighting_function_vectorized = identity_weighting_function_vectorized,
                               g_function_vectorized = identity_g_function_vectorized) {
  wg = wg_vectorized(s1, S, weighting_function_vectorized, g_function_vectorized)
  # wg (n_s1, n, dim_g)
  l = int_s1(s0, psp_s1, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # l (n_s1, n)
  # for eta1
  mu1 = lm_prediction_vectorized(mu1_lm, s1, X)
  # mu1 (n_s1, n)
  l_eta1 = l * mu1
  dim(l_eta1) = c(dim(l_eta1), 1)
  l_eta1 = l_eta1[, , rep(1, dim(wg)[3])]
  # for eta0
  l_eta0 = sweep(l, MARGIN = 2, Y, `*`)
  dim(l_eta0) = c(dim(l_eta0), 1)
  l_eta0 = l_eta0[, , rep(1, dim(wg)[3])]
  return(list(eta1 = wg * l_eta1, eta0 = wg * l_eta0))
}
