#' Point estimator for principal causal effects
#'
#' This function estimates principal causal effects using principal densities, treatment
#' probabilities, and outcome models. The method combines an efficient influence
#' function estimator (eif), a treatment probability plus principal density estimator (tp_pd), and a
#' principal density plus outcome modeling estimator (pd_om). The function integrates over possible
#' values of the post-treatment variable to compute point estimates for the coefficients eta1 and eta0.
#'
#' @param Z A numeric vector indicating the treatment assignment (1 for treatment,
#' 0 for control).
#' @param X A matrix of covariates.
#' @param S A numeric vector representing the post-treatment variable.
#' @param Y A numeric vector of outcomes.
#' @param n_divisions Integer specifying the number of divisions for integration.
#' Default is 100.
#' @param copula_type A string indicating the type of copula to be used.
#' Default is 'gaussian'.
#' @param rho A numeric parameter for the copula, the correlation coefficient, default is 0.
#' @param weighting_function_vectorized A vectorized weighting function used in
#' the calculation of estimators. Default is `identity_weighting_function_vectorized`.
#' @param g_function_vectorized A vectorized g-function used in the calculation of
#' estimators. Default is `identity_g_function_vectorized`.
#'
#' @return A numeric vector of the point estimates for the projection coefficients.
#' The vector contains the differences between the estimators for eta1 and eta0, for eif, tp_pd, pd_om estimators
#'
#' @examples
#' # Example usage
#' set.seed(42)
#' Z <- sample(0:1, 100, replace = TRUE)
#' X <- matrix(rnorm(100*3), ncol = 3)
#' S <- rnorm(100)
#' Y <- rnorm(100)
#' result <- point_estimator(Z, X, S, Y)
#' print(result)
#'
#' @import abind
#' @importFrom np npcdensbw
#' @importFrom pracma linspace inv
#'
#' @export
#'
point_estimator <- function(Z,
                            X,
                            S,
                            Y,
                            n_divisions = 100,
                            copula_type = 'gaussian',
                            rho = 0,
                            weighting_function_vectorized = identity_weighting_function_vectorized,
                            g_function_vectorized = identity_g_function_vectorized) {
  # 1. principal score model
  bw_treated <- npcdensbw(xdat = X[which(Z == 1), ], ydat = S[which(Z == 1)])
  bw_control <- npcdensbw(xdat = X[which(Z == 0), ], ydat = S[which(Z == 0)])

  # 2. treatment probability model
  treatment_probability_logit <- stats::glm(Z ~ as.matrix(X), family = 'binomial')
  Tp <- stats::predict(treatment_probability_logit, type = 'response')

  # 3. outcome model
  S_X <- data.frame(cbind(S, data.frame(X)))
  colnames(S_X) <- c('S', paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
  mu1_lm <- stats::lm(Y ~ ., data = S_X, weights = Z)
  Mu1 <- stats::predict(mu1_lm)
  mu0_lm <- stats::lm(Y ~ ., data = S_X, weights = 1 - Z)
  Mu0 <- stats::predict(mu0_lm)

  # integrating points
  lower <- min(S)
  upper <- max(S)
  s1 <- linspace(lower, upper, n_divisions)
  s0 <- linspace(lower, upper, n_divisions)

  # conditional density and cdfs
  psp_s1 = principal_score_predict_vectorized(bw_treated, s1, X)
  # psp_s1 (n_s1, n)
  psp_s0 = principal_score_predict_vectorized(bw_control, s0, X)
  # psp_s0 (n_s0, n)
  psp_s1_s0 = matrix_multiply_with_expansion(psp_s1, psp_s0, 'both')
  # psp_s1_s0 (n_s1, n_s0, n)
  psp_s1_cdf = trapz_cdf(s1, psp_s1)
  # psp_s1_cdf (n_s1, n)
  psp_s0_cdf = trapz_cdf(s0, psp_s0)
  # psp_s0_cdf (n_s0, n)


  # eif_estimator
  B_int2_integrand <- colMeans(trapz_vector(s0, trapz_vector(
    s1,
    B_int2_function(
      s1,
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
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  )))
  B_ints0_integrand <- colMeans(trapz_vector(
    s0,
    B_ints0_function(
      s1,
      s0,
      psp_s0,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))
  B_ints1_integrand <- colMeans(trapz_vector(
    s1,
    B_ints1_function(
      s1,
      s0,
      psp_s1,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))

  vec_int2 <- vec_int2_function(
    mu1_lm,
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
    copula_type,
    rho,
    weighting_function_vectorized,
    g_function_vectorized
  )
  vec_int_to_s0 <- vec_ints0_function(
    mu0_lm,
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
    copula_type,
    rho,
    weighting_function_vectorized,
    g_function_vectorized
  )
  vec_int_to_s1 <- vec_ints1_function(
    mu1_lm,
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
    copula_type,
    rho,
    weighting_function_vectorized,
    g_function_vectorized
  )
  # for eta_1
  vec_int2_integrand_eta1 <- colMeans(trapz_vector(s0, trapz_vector(s1, vec_int2$eta1)))
  vec_int_observed_integrand_eta1 <- colMeans(trapz_vector(s0, vec_int_to_s0$eta1))
  vec_int_counterfactual_integrand_eta1 <- colMeans(trapz_vector(s1, vec_int_to_s1$eta1))
  # for eta_0
  vec_int2_integrand_eta0 <- colMeans(trapz_vector(s1, trapz_vector(s0, vec_int2$eta0)))
  vec_int_observed_integrand_eta0 <- colMeans(trapz_vector(s1, vec_int_to_s1$eta0))
  vec_int_counterfactual_integrand_eta0 <- colMeans(trapz_vector(s0, vec_int_to_s0$eta0))

  B_eif <- B_int2_integrand + B_ints0_integrand + B_ints1_integrand
  vec_eif_eta1 <- vec_int2_integrand_eta1 + vec_int_observed_integrand_eta1 + vec_int_counterfactual_integrand_eta1
  vec_eif_eta0 <- vec_int2_integrand_eta0 + vec_int_observed_integrand_eta0 + vec_int_counterfactual_integrand_eta0
  eif_est_eta1 <- inv(B_eif) %*% vec_eif_eta1
  eif_est_eta0 <- inv(B_eif) %*% vec_eif_eta0

  # tp_pd_estimator
  B_tp_pd_eta1 <- colMeans(trapz_vector(
    s0,
    B_ints0_tp_pd_function(
      s1,
      s0,
      psp_s0,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))
  B_tp_pd_eta0 <- colMeans(trapz_vector(
    s1,
    B_ints1_tp_pd_function(
      s1,
      s0,
      psp_s1,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))
  vec_tp_pd_eta1 <- colMeans(trapz_vector(
    s0,
    vec_ints0_tp_pd_function(
      s1,
      s0,
      psp_s0,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      Y,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))
  vec_tp_pd_eta0 <- colMeans(trapz_vector(
    s1,
    vec_ints1_tp_pd_function(
      s1,
      s0,
      psp_s1,
      psp_s1_cdf,
      psp_s0_cdf,
      Z,
      Tp,
      S,
      Y,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  ))
  tp_pd_est_eta1 <- inv(B_tp_pd_eta1) %*% vec_tp_pd_eta1
  tp_pd_est_eta0 <- inv(B_tp_pd_eta0) %*% vec_tp_pd_eta0

  # pd_om_estimator
  B_pd_om <- colMeans(trapz_vector(s0, trapz_vector(
    s1,
    B_int2_pd_om_function(
      s1,
      s0,
      psp_s1_s0,
      psp_s1_cdf,
      psp_s0_cdf,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  )))
  vec_int2_pd_om <- vec_int2_pd_om_function(
    mu1_lm,
    mu0_lm,
    s1,
    s0,
    psp_s1_s0,
    psp_s1_cdf,
    psp_s0_cdf,
    X,
    copula_type,
    rho,
    weighting_function_vectorized,
    g_function_vectorized
  )
  # for eta1
  vec_pd_om_eta1 <- colMeans(trapz_vector(s0, trapz_vector(s1, vec_int2_pd_om$eta1)))
  # for eta0
  vec_pd_om_eta0 <- colMeans(trapz_vector(s1, trapz_vector(s0, vec_int2_pd_om$eta0)))
  pd_om_est_eta1 <- inv(B_pd_om) %*% vec_pd_om_eta1
  pd_om_est_eta0 <- inv(B_pd_om) %*% vec_pd_om_eta0

  result <- as.numeric(cbind(
    t(eif_est_eta1 - eif_est_eta0),
    t(tp_pd_est_eta1 - tp_pd_est_eta0),
    t(pd_om_est_eta1 - pd_om_est_eta0)
  ))

  return(result)
}


#' Nonparametric Bootstrap Estimation for Treatment Effects
#'
#' This function performs nonparametric bootstrap to estimate the standard errors
#' of projection coefficients of principal causal effect estimators. The bootstrap samples are drawn with replacement
#' from the data, and for each sample, the `point_estimator` function is used to
#' estimate the projection coefficients. The function returns both the point estimates and
#' the bootstrap standard errors.
#'
#' @param Z A numeric vector indicating the treatment assignment (1 for treatment,
#' 0 for control).
#' @param X A matrix of covariates.
#' @param S A numeric vector representing the post-treatment variable.
#' @param Y A numeric vector of outcomes.
#' @param n_boot Integer specifying the number of bootstrap iterations.
#' Default is 500.
#' @param n_divisions Integer specifying the number of divisions for integration.
#' Default is 100.
#' @param copula_type A string indicating the type of copula to be used.
#' Default is 'gaussian'.
#' @param rho A numeric parameter for the copula, the correlation coefficient, default is 0.
#' @param weighting_function_vectorized A vectorized weighting function used in
#' the calculation of estimators. Default is `identity_weighting_function_vectorized`.
#' @param g_function_vectorized A vectorized g-function used in the calculation of
#' estimators. Default is `identity_g_function_vectorized`.
#'
#' @return A list with the following elements:
#' \item{point_est}{The point estimates for projection coefficients,
#' as returned by the `point_estimator` function.}
#' \item{boot_est}{A matrix of bootstrap estimates across the `n_boot` iterations.}
#' \item{res}{A matrix with two rows: the first row contains the point estimates,
#' and the second row contains the bootstrap standard errors.}
#'
#' @examples
#' # Example usage
#' set.seed(42)
#' Z <- sample(0:1, 100, replace = TRUE)
#' X <- matrix(rnorm(100*3), ncol = 3)
#' S <- rnorm(100)
#' Y <- rnorm(100)
#' boot_result <- boot(Z, X, S, Y, n_boot=5)
#' print(boot_result$res)
#'
#' @export
#'
boot <- function(Z,
                 X,
                 S,
                 Y,
                 n_boot = 500,
                 n_divisions = 100,
                 copula_type = 'gaussian',
                 rho = 0,
                 weighting_function_vectorized = identity_weighting_function_vectorized,
                 g_function_vectorized = identity_g_function_vectorized) {
  point_est <- point_estimator(
    Z,
    X,
    S,
    Y,
    n_divisions,
    copula_type,
    rho,
    weighting_function_vectorized,
    g_function_vectorized
  )

  # nonparametric bootstrap
  n <- length(Z)
  X <- as.matrix(X)
  boot_est <- replicate(n_boot, {
    id_boot = sample(1:n, n, replace = TRUE)
    point_estimator(
      Z[id_boot],
      X[id_boot, ],
      S[id_boot],
      Y[id_boot],
      n_divisions,
      copula_type,
      rho,
      weighting_function_vectorized,
      g_function_vectorized
    )
  })

  boot_se <- apply(data.frame(boot_est), 1, stats::sd)

  res <- rbind(point_est, boot_se)
  rownames(res) <- c("est", "boot_se")
  return(list(
    point_est = point_est,
    boot_est = boot_est,
    res = res
  ))
}
