# for eta1
B_ints0_tp_pd_function <- function(s1,
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
  wggt = aperm(
    wggt_vectorized(S, s0, weighting_function_vectorized, g_function_vectorized),
    c(2, 1, 3, 4)
  )
  # wggt (n_s0, n, dim_g, dim_g)
  joint_score = int_s0(s1, psp_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # joint_score (n_s0, n)
  dim(joint_score) = c(dim(joint_score), 1, 1)
  joint_score = joint_score[, , rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * joint_score)
}

vec_ints0_tp_pd_function <- function(s1,
                                     s0,
                                     psp_s0,
                                     psp_s1_cdf,
                                     psp_s0_cdf,
                                     Z,
                                     Tp,
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
  # wg (n_s0, n, dim_g)
  joint_score = int_s0(s1, psp_s0, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # joint_score (n_s0, n)
  z_joint_score = sweep(joint_score, MARGIN = 2, Y, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1)
  z_joint_score = z_joint_score[, , rep(1, dim(wg)[3])]
  return(wg * z_joint_score)
}

# for eta0
B_ints1_tp_pd_function <- function(s1,
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
  wggt = wggt_vectorized(s1, S, weighting_function_vectorized, g_function_vectorized)
  # wggt (n_s1, n, dim_g, dim_g)
  joint_score = int_s1(s0, psp_s1, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # joint_score (n_s1, n)
  dim(joint_score) = c(dim(joint_score), 1, 1)
  joint_score = joint_score[, , rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * joint_score)
}

vec_ints1_tp_pd_function <- function(s1,
                                     s0,
                                     psp_s1,
                                     psp_s1_cdf,
                                     psp_s0_cdf,
                                     Z,
                                     Tp,
                                     S,
                                     Y,
                                     copula_type = 'independent',
                                     rho = 0,
                                     weighting_function_vectorized = identity_weighting_function_vectorized,
                                     g_function_vectorized = identity_g_function_vectorized) {
  wg = wg_vectorized(s1, S, weighting_function_vectorized, g_function_vectorized)
  # wg (n_s1, n, dim_g)
  joint_score = int_s1(s0, psp_s1, psp_s1_cdf, psp_s0_cdf, Z, Tp, S, copula_type, rho)
  # joint_score (n_s1, n)
  z_joint_score = sweep(joint_score, MARGIN = 2, Y, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1)
  z_joint_score = z_joint_score[, , rep(1, dim(wg)[3])]
  return(wg * z_joint_score)
}
