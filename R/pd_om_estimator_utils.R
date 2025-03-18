int2_pd_om <- function(psp_s1_s0,
                       psp_s1_cdf,
                       psp_s0_cdf,
                       copula_type = 'independent',
                       rho = 0) {
  # return(n_s1, n_s0, n)
  l1_B = copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'both', copula_type, rho)
  return(l1_B * psp_s1_s0)
}

B_int2_pd_om_function <- function(s1,
                                  s0,
                                  psp_s1_s0,
                                  psp_s1_cdf,
                                  psp_s0_cdf,
                                  copula_type = 'independent',
                                  rho = 0,
                                  weighting_function_vectorized = identity_weighting_function_vectorized,
                                  g_function_vectorized = identity_g_function_vectorized) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1,
                         s0,
                         weighting_function_vectorized,
                         g_function_vectorized)
  # wggt (n_s1, n_s0, dim_g, dim_g)
  l1 = int2_pd_om(psp_s1_s0, psp_s1_cdf, psp_s0_cdf, copula_type, rho)
  # l1 (n_s1, n_s0, n)
  dim(wggt) = c(dim(wggt)[1], dim(wggt)[2], 1, dim(wggt)[3], dim(wggt)[4])
  wggt = wggt[, , rep(1, dim(l1)[3]), , ]
  dim(l1) = c(dim(l1), 1, 1)
  l1 = l1[, , , rep(1, dim(wggt)[4]), rep(1, dim(wggt)[5])]
  return(wggt * l1)
}

vec_int2_pd_om_function <- function(mu1_lm,
                                    mu0_lm,
                                    s1,
                                    s0,
                                    psp_s1_s0,
                                    psp_s1_cdf,
                                    psp_s0_cdf,
                                    X,
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
  # l1
  l1 = int2_pd_om(psp_s1_s0, psp_s1_cdf, psp_s0_cdf, copula_type, rho)
  # l1 (n_s1, n_s0, n)
  dim(wg) = c(dim(wg)[1], dim(wg)[2], 1, dim(wg)[3])
  wg = wg[, , rep(1, dim(X)[1]), ]
  dim(l1) = c(dim(l1), 1)
  l1 = l1[, , , rep(1, dim(wg)[4])]
  dim(mu1) = c(dim(mu1)[1], 1, dim(mu1)[2], 1)
  dim(mu0) = c(1, dim(mu0)[1], dim(mu0)[2], 1)
  mu1 = mu1[, rep(1, length(s0)), , rep(1, dim(wg)[4])]
  mu0 = mu0[rep(1, length(s1)), , , rep(1, dim(wg)[4])]
  return(list(eta1 = wg * l1 * mu1, eta0 = wg * l1 * mu0))
}
