#' @import np
#' @import ks
#' @import abind
#' @import plyr
#' @import dplyr
NULL
#> NULL

#' Two matrix multiplication with different first dimension
#' @param U A matrix of size (n_u, n)
#' @param V A matrix of size (n_v, n)
#' @return A matrix of size (n_u, n_v, n)
matrix_multiply_with_expansion <- function(U, V, expansion_type = 'both') {
  switch (expansion_type,
          'U' = {
            # U (n_u, n)
            # V (n_u, n_v, n)
            # return (n_u, n_v, n)
            dim(U) = c(dim(U)[1], 1, dim(U)[2])
            U = U[, rep(1, dim(V)[2]), ]
            return(U * V)
          },
          'V' = {
            # U (n_u, n_v, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
            dim(V) = c(1, dim(V))
            V = V[rep(1, dim(U)[1]), , ]
            return(U * V)
          },
          'both' = {
            # U (n_u, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
            dim(U) = c(dim(U)[1], 1, dim(U)[2])
            U = U[, rep(1, dim(V)[1]), ]
            dim(V) = c(1, dim(V))
            V = V[rep(1, dim(U)[1]), , ]
            return(U * V)
          }
  )     
  # U (n_u, n)
  # V (n_v, n)
  # return (n_u, n_v, n)
  
  dim(V) = c(1, dim(V))
  V = V[rep(1, dim(U)[1]), , ]
  return(U * V)
}

#' Compute integral in vector form
#' @param x A vector 
#' @param y A vector or a matrix
#' @return if y is a vector, return a number else return a vector
trapz_vector <- function(x, y, axis=1) {
  idx = 2:length(x)
  if(is.null(dim(y))) {
    return(sum((x[idx] - x[idx-1]) * (asub(y, idx, dims=axis) + asub(y, idx-1, dims=axis)) / 2))
  }
  return(colSums((x[idx] - x[idx-1]) * (asub(y, idx, dims=axis) + asub(y, idx-1, dims=axis)) / 2))
}

#' Compute cdf in vector form
#' @param x A vector 
#' @param y A vector or a matrix
#' @return if y is a vector, return a number else return a vector
trapz_cdf <- function(x, y, axis=1) {
  idx = 2:length(x)
  if(is.null(dim(y))) {
    # x (n_x)
    # y (n_x)
    # return (n_x)
    return(c(0, cumsum((x[idx] - x[idx-1]) * (asub(y, idx, dims=axis) + asub(y, idx-1, dims=axis)) / 2)))
  }
  # x (n_x)
  # y (n_x, n)
  # return (n_x, n)
  return(rbind(rep(0, dim(y)[2]), apply(((x[idx] - x[idx-1]) * (asub(y, idx, dims=axis) + asub(y, idx-1, dims=axis)) / 2), 2, cumsum)))
}

# ------------------------------- #
# --- User provided functions --- #
# ------------------------------- #

# User provided
# Gaussian copula
gaussian_copula_function <- function(u, v, rho) { 
  # rho \in (-1, 1) is correlation
  x_u <- qnorm(u)
  y_v <- qnorm(v)
  c <- (sqrt(1 - rho ^ 2) * exp(-(rho^2 * (x_u^2 + y_v^2) - 2*rho*x_u*y_v) / (2*(1 - rho ^ 2))))
  return(c)
}

# FGM copula (vectorized version to be added)
fgm_copula_function <- function(u, v, rho) {
  # rho \in [-1, 1]
  c <- 1 + rho * (1 - 2*u) * (1 - 2*v)
  return(c)
}

# copula functions vectorized
# independent copula
copula_function_vectorized <- function(U, V, expansion_type='both', copula_type='independent', rho=0) {
  if (copula_type == 'independent') {
    switch (expansion_type,
            'U' = return(array(1, dim(V))),
            # U (n)
            # V (n_s, n)
            # return (n_s, n)
            'V' = return(array(1, dim(U))),
            # U (n_s, n)
            # V (n)
            # return (n_s, n)
            'both' = return(array(1, c(dim(U)[1], dim(V))))
            # U (n_u, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
    )
  }
  U[U <= 0] = .Machine$double.eps
  U[U >= 1] = 1 - .Machine$double.eps
  V[V <= 0] = .Machine$double.eps
  V[V >= 1] = 1 - .Machine$double.eps
  switch (expansion_type,
          'U' = {
            # U (n)
            # V (n_s, n)
            # return (n_s, n)
            U = t(matrix(rep(U, dim(V)[1]), ncol=dim(V)[1]))
          },
          'V' = {
            # U (n_s, n)
            # V (n)
            # return (n_s, n)
            V = t(matrix(rep(V, dim(U)[1]), ncol=dim(U)[1]))
          },
          'both' = {
            # U (n_u, n)
            # V (n_v, n)
            # return (n_u, n_v, n)
            dim(U) = c(dim(U)[1], 1, dim(U)[2])
            U = U[, rep(1, dim(V)[1]), ]
            dim(V) = c(1, dim(V)[1], dim(V)[2])
            V = V[rep(1, dim(U)[1]),,]
          }
  )
  switch (copula_type,
          'gaussian' = return(gaussian_copula_function(U, V, rho)),
          'fgm'= return(fgm_copula_function(U, V, rho))
  )
}


# User provided
gaussian_copula_gradient <- function(u, v, rho) { 
  # rho \in (-1, 1) is correlation
  x_u <- qnorm(u)
  y_v <- qnorm(v)
  c_u <- gaussian_copula_function(u, v, rho) * (rho * y_v - rho ^ 2 * x_u) / (1 - rho ^ 2) / dnorm(x_u)
  c_v <- gaussian_copula_function(u, v, rho) * (rho * x_u - rho ^ 2 * y_v) / (1 - rho ^ 2) / dnorm(y_v)
  return(list(cu=c_u, cv=c_v))
}

fgm_copula_gradient <- function(u, v, rho) {
  # rho \in [-1, 1]
  c_u <- -2 * rho * (1 - 2 * v)
  c_v <- -2 * rho * (1 - 2 * u)
  return(list(cu=c_u, cv=c_v))
}

copula_gradient_vectorized <- function(U, V, copula_type='independent', rho=0) {
  # U (n_u, n)
  # V (n_v, n)
  if (copula_type == 'independent') {
    return(list(cu=array(0, c(dim(U)[1], dim(V)[1], dim(U)[2])), 
                cv=array(0, c(dim(U)[1], dim(V)[1], dim(U)[2]))))
  }
  U[U <= 0] = .Machine$double.eps
  U[U >= 1] = 1 - .Machine$double.eps
  V[V <= 0] = .Machine$double.eps
  V[V >= 1] = 1 - .Machine$double.eps
  dim(U) = c(dim(U)[1], 1, dim(U)[2])
  U = U[, rep(1, dim(V)[1]), ]
  dim(V) = c(1, dim(V)[1], dim(V)[2])
  V = V[rep(1, dim(U)[1]),,]
  switch (copula_type,
          'gaussian' = return(gaussian_copula_gradient(U, V, rho)),
          'fgm'= return(gaussian_copula_gradient(U, V, rho))
  )
  
}

# User provided
weighting_function_vectorized <- function(s1, s0) {
  # return (n_s1, n_s0)
  return(array(1, c(length(s1), length(s0))))
}

# User provided
g_function_vectorized <- function(s1, s0) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, dim_g)
  return(abind(replicate(length(s0), s1), t(replicate(length(s1), s0)), array(1, c(length(s1), length(s0))), along=3))
}


# --------------------------------------- #
# --- Useful functions for estimation --- #
# --------------------------------------- #

principal_score_models <- function(S, X) {
  return(npcdensbw(xdat=X, ydat=S))
}

principal_score_predict_array <- function(bw, S, X) {
  return(fitted(npcdens(bws=bw, exdat=X, eydat=S)))
}

principal_score_predict_vectorized <- function(bw, S, X) {
  # S (n_s, )
  # X (n, dim_x)
  # return (n_s, n)
  return(t(array(fitted(npcdens(
    bws=bw,
    exdat=do.call(rbind, replicate(length(S), X, simplify=FALSE)),
    eydat=cbind(as.vector(t(replicate(dim(X)[1], S)))))), c(dim(X)[1], length(S)))))
}



mu_function_array <- function(mu_model, S, X) {
  newdata = data.frame(cbind(S, data.frame(X)))
  colnames(newdata) = c('S',  paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
  return(predict.lm(mu_model, newdata=newdata))
}

mu_function_vectorized <- function(mu_model, S, X) {
  # return (n_s, n)
  newdata = data.frame(cbind(as.vector(t(replicate(dim(X)[1], S))),
                             do.call(rbind, replicate(length(S), X, simplify=FALSE))))
  colnames(newdata) = c('S',  paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
  return(t(array(predict.lm(mu_model, newdata=newdata), c(dim(X)[1], length(S)))))
}

wggt_vectorized <- function(s1, s0) {
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
  return(w * g1 *g2)
}

wg_vectorized <- function(s1, s0) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # return (n_s1, n_s0, dim_g)
  w = weighting_function_vectorized(s1, s0)
  g = g_function_vectorized(s1, s0)
  dim(w) = c(dim(w), 1)
  w = w[, , rep(1, dim(g)[3])]
  return(w * g)
}

# ---------------------------------------------- #
# --- functions to compute the eif estimator --- #
# ---------------------------------------------- #
int2 <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, S, copula_type='independent', rho=0) {
  # return(n_s1, n_s0, n)
  psp_s1 = principal_score_predict_vectorized(bw_treated, s1, X)
  # psp_s1 (n_s1, n)
  psp_s0 = principal_score_predict_vectorized(bw_control, s0, X)
  # psp_s0 (n_s0, n)
  psp_s1_cdf = trapz_cdf(s1, psp_s1, axis=1)
  # psp_s1_cdf (n_s1, n)
  psp_s0_cdf = trapz_cdf(s0, psp_s0, axis=1)
  # psp_s0_cdf (n_s0, n)
  l1_B = copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'both', copula_type, rho)
  copula_gradient = copula_gradient_vectorized(psp_s1_cdf, psp_s0_cdf, copula_type, rho)
  part2_s1 = l1_B + matrix_multiply_with_expansion(psp_s1_cdf - 1 * outer(s1, S, FUN=">="), copula_gradient$cu, 'U')
  # part2_s1 (n_s1, n_s0, n)
  part2_s0 = l1_B + matrix_multiply_with_expansion(copula_gradient$cv, psp_s0_cdf - 1 * outer(s0, S, FUN=">="), 'V')
  # part2_s0 (n_s1, n_s0, n)
  lp = l1_B - sweep(part2_s1, MARGIN=3, Z / Tp, `*`) - sweep(part2_s0, MARGIN=3, (1 - Z) / (1 - Tp), `*`)
  return(lp * matrix_multiply_with_expansion(psp_s1, psp_s0, 'both'))
}

int_s0 <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type='independent', rho=0) {
  # s1 = S
  # return (n_s0, n)
  psp_s1 = principal_score_predict_array(bw_treated, s1, X)
  # psp_s1 (n)
  psp_s1_cdf = trapz_cdf(s1, psp_s1)
  # psp_s1_cdf (n)
  psp_s0 = principal_score_predict_vectorized(bw_control, s0, X)
  # psp_s0 (n_s0, n)
  psp_s0_cdf = trapz_cdf(s0, psp_s0, axis=1)
  # psp_s0_cdf (n_s0, n)
  part2 = sweep(copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'U', copula_type, rho), MARGIN=2, Z / Tp, `*`)
  return(part2 * psp_s0)
}

int_s1 <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type='independent', rho=0) {
  # s0 = S
  # return (n_s1, n)
  psp_s1 = principal_score_predict_vectorized(bw_treated, s1, X)
  # psp_s1 (n_s1, n)
  psp_s1_cdf = trapz_cdf(s1, psp_s1, axis=1)
  # psp_s1_cdf (n_s1, n)
  psp_s0 = principal_score_predict_array(bw_control, s0, X)
  # psp_s0 (n)
  psp_s0_cdf = trapz_cdf(s0, psp_s0)
  # psp_s0_cdf (n)
  part2 = sweep(copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'V', copula_type, rho), MARGIN=2, (1 - Z) / (1 - Tp), `*`)
  # part2 (n_s1, n)
  return(part2 * psp_s1)
}

B_int2_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, S, copula_type='independent', rho=0) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # X (n, dim_x)
  # return (n_s1, n_s0, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1, s0)
  # wggt (n_s1, n_s0, dim_g, dim_g)
  lp = int2(bw_treated, bw_control, s1, s0, Z, Tp, X, S, copula_type, rho)
  dim(wggt) = c(dim(wggt)[1], dim(wggt)[2], 1, dim(wggt)[3], dim(wggt)[4])
  wggt = wggt[,,rep(1, length(Z)),,]
  dim(lp) = c(dim(lp), 1, 1)
  lp = lp[,,,rep(1, dim(wggt)[4]),rep(1, dim(wggt)[5])]
  return(wggt * lp)
}

B_ints0_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type='independent', rho=0) {
  # s1 = S
  # return (n_s0, n, dim_g, dim_g)
  wggt = aperm(wggt_vectorized(s1, s0), c(2, 1, 3, 4))
  part2 = int_s0(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  dim(part2) = c(dim(part2), 1, 1)
  part2 = part2[,,rep(1, dim(wggt)[3]),rep(1, dim(wggt)[4])]
  return(wggt * part2)
}

B_ints1_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho) {
  # return (n_s1, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1, s0)
  # wggt (n_s1, n_s0=n, dim_g, dim_g)
  part2 = int_s1(bw_treated, bw_control, s1, s0, Z, Tp, X)
  dim(part2) = c(dim(part2), 1, 1)
  part2 = part2[,,rep(1, dim(wggt)[3]),rep(1, dim(wggt)[4])]
  return(wggt * part2)
}


vec_int2_function <- function(bw_treated, bw_control, mu1_lm, mu0_lm, s1, s0, Z, Tp, X, S, copula_type='independent', rho=0) {
  mu1 = mu_function_vectorized(mu1_lm, s1, X)
  # mu1 (n_s1, n)
  mu0 = mu_function_vectorized(mu0_lm, s0, X)
  # mu1 (n_s0, n)
  wg = wg_vectorized(s1, s0)
  # wg (n_s1, n_s0, dim_g)
  lp = int2(bw_treated, bw_control, s1, s0, Z, Tp, X, S, copula_type, rho)
  dim(wg) = c(dim(wg)[1], dim(wg)[2], 1, dim(wg)[3])
  wg = wg[,,rep(1, length(Z)),]
  dim(lp) = c(dim(lp), 1)
  lp = lp[,,,rep(1, dim(wg)[4])]
  # lp (n_s1, n_s0, n)
  wg_lp  = wg * lp
  # wglp (n_s1, n_s0, n, dim_g)
  dim(mu1) = c(dim(mu1)[1], 1, dim(mu1)[2], 1)
  dim(mu0) = c(1, dim(mu0)[1], dim(mu0)[2], 1)
  mu1 = mu1[, rep(1, length(s0)),,rep(1, dim(wg_lp)[4])]
  mu0 = mu0[rep(1, length(s1)),,,rep(1, dim(wg_lp)[4])]
  return(list(eta1=wg_lp * mu1, 
              eta0=wg_lp * mu0))
}


vec_int_to_s0_function <- function(bw_treated, bw_control, mu0_lm, s1, s0, Z, Tp, X, Y, copula_type='independent', rho=0) {
  # (n_s1 = n, )
  wg = aperm(wg_vectorized(s1, s0), c(2, 1, 3))
  # wg (n_s0, n_s1 = n, dim_g)
  l = int_s0(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # l(n_s0, n)
  # for eta1
  l_eta1 = sweep(l, MARGIN=2, Y, `*`)
  dim(l_eta1) = c(dim(l_eta1), 1)
  l_eta1 = l_eta1[,,rep(1, dim(wg)[3])]
  # for eta0
  mu0 = mu_function_vectorized(mu0_lm, s0, X)
  # mu0 (n_s0, n)
  l_eta0 = l * mu0
  dim(l_eta0) = c(dim(l_eta0), 1)
  l_eta0 = l_eta0[,,rep(1, dim(wg)[3])]
  return(list(eta1=wg * l_eta1, 
              eta0=wg * l_eta0))
}

vec_int_to_s1_function <- function(bw_treated, bw_control, mu1_lm, s1, s0, Z, Tp, X, Y, copula_type='independent', rho=0) {
  # (n_s0 = n, )
  wg = wg_vectorized(s1, s0)
  # wg (n_s1, n_s0 = n, dim_g)
  l = int_s1(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # l (n_s1, n)
  # for eta1
  mu1 = mu_function_vectorized(mu1_lm, s1, X)
  # mu1 (n_s1, n)
  l_eta1 = l * mu1
  dim(l_eta1) = c(dim(l_eta1), 1)
  l_eta1 = l_eta1[,,rep(1, dim(wg)[3])]
  # for eta0
  l_eta0 = sweep(l, MARGIN=2, Y, `*`)
  dim(l_eta0) = c(dim(l_eta0), 1)
  l_eta0 = l_eta0[,,rep(1, dim(wg)[3])]
  return(list(eta1=wg * l_eta1, 
              eta0=wg * l_eta0))
}

# ------------------------------------------------ #
# --- functions to compute the tp_ps estimator --- #
# ------------------------------------------------ #

# for eta1
B_ints0_tp_ps_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, P1, copula_type='independent', rho=0) {
  wggt = aperm(wggt_vectorized(s1, s0), c(2, 1, 3, 4))
  # wggt (n_s0, n_s1=n, dim_g, dim_g)
  joint_score = int_s0(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # joint_score (n_s0, n)
  psp_s1 = principal_score_predict_array(bw_treated, s1, X)
  # psp_s1 (n)
  z_joint_score = sweep(joint_score, MARGIN=2, psp_s1 / P1, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1, 1)
  z_joint_score=z_joint_score[,,rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * z_joint_score)
}

vec_ints0_tp_ps_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, Y, P1, copula_type='independent', rho=0) {
  # (n_s1 = n, )
  wg = aperm(wg_vectorized(s1, s0), c(2, 1, 3))
  # wg (n_s0, n_s1 = n, dim_g)
  joint_score = int_s0(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # joint_score (n_s0, n)
  psp_s1 = principal_score_predict_array(bw_treated, s1, X)
  # psp_s1 (n)
  z_joint_score = sweep(joint_score, MARGIN=2, Y * psp_s1 / P1, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1)
  z_joint_score=z_joint_score[,,rep(1, dim(wg)[3])]
  return(wg * z_joint_score)
}

# for eta0
B_ints1_tp_ps_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, P0, copula_type='independent', rho=0) {
  wggt = wggt_vectorized(s1, s0)
  # wggt (n_s1, n_s0=n, dim_g, dim_g)
  joint_score = int_s1(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # joint_score (n_s1, n)
  psp_s0 = principal_score_predict_array(bw_control, s0, X)
  # psp_s0 (n)
  z_joint_score = sweep(joint_score, MARGIN=2, psp_s0 / P0, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1, 1)
  z_joint_score=z_joint_score[,,rep(1, dim(wggt)[3]), rep(1, dim(wggt)[4])]
  return(wggt * z_joint_score)
}

vec_ints1_tp_ps_function <- function(bw_treated, bw_control, s1, s0, Z, Tp, X, Y, P0, copula_type='independent', rho=0) {
  # (n_s0 = n, )
  wg = wg_vectorized(s1, s0)
  # wg (n_s1, n_s0 = n, dim_g)
  joint_score = int_s1(bw_treated, bw_control, s1, s0, Z, Tp, X, copula_type, rho)
  # joint_score (n_s1, n)
  psp_s0 = principal_score_predict_array(bw_control, s0, X)
  # psp_s0 (n)
  z_joint_score = sweep(joint_score, MARGIN=2, Y * psp_s0 / P0, `*`)
  dim(z_joint_score) = c(dim(z_joint_score), 1)
  z_joint_score=z_joint_score[,,rep(1, dim(wg)[3])]
  return(wg * z_joint_score)
}

# ------------------------------------------------ #
# --- functions to compute the ps_om estimator --- #
# ------------------------------------------------ #
int2_ps_om <- function(bw_treated, bw_control, s1, s0, X, copula_type='independent', rho=0) {
  # return(n_s1, n_s0, n)
  psp_s1 = principal_score_predict_vectorized(bw_treated, s1, X)
  # psp_s1 (n_s1, n)
  psp_s1_cdf = trapz_cdf(s1, psp_s1, axis=1)
  # psp_s1_cdf (n_s1, n)
  psp_s0 = principal_score_predict_vectorized(bw_control, s0, X)
  # psp_s0 (n_s0, n)
  psp_s0_cdf = trapz_cdf(s0, psp_s0, axis=1)
  # psp_s0_cdf (n_s0, n)
  l1_B = copula_function_vectorized(psp_s1_cdf, psp_s0_cdf, 'both', copula_type, rho)
  return(l1_B * matrix_multiply_with_expansion(psp_s1, psp_s0, 'both'))
}

B_int2_ps_om_function <- function(bw_treated, bw_control, s1, s0, X, copula_type='independent', rho=0) {
  # s1 (n_s1, )
  # s0 (n_s0, )
  # X (n, dim_x)
  # return (n_s1, n_s0, n, dim_g, dim_g)
  wggt = wggt_vectorized(s1, s0)
  # wggt (n_s1, n_s0, dim_g, dim_g)
  l1 = int2_ps_om(bw_treated, bw_control, s1, s0, X, copula_type, rho)
  # l1 (n_s1, n_s0, n)
  dim(wggt) = c(dim(wggt)[1], dim(wggt)[2], 1, dim(wggt)[3], dim(wggt)[4])
  wggt = wggt[,,rep(1, dim(X)[1]),,]
  dim(l1) = c(dim(l1), 1, 1)
  l1 = l1[,,,rep(1, dim(wggt)[4]),rep(1, dim(wggt)[5])]
  return(wggt * l1)
}

vec_int2_ps_om_function <- function(bw_treated, bw_control, mu1_lm, mu0_lm, s1, s0, X, copula_type='independent', rho=0) {
  mu1 = mu_function_vectorized(mu1_lm, s1, X)
  # mu1 (n_s1, n)
  mu0 = mu_function_vectorized(mu0_lm, s0, X)
  # mu1 (n_s0, n)
  wg = wg_vectorized(s1, s0)
  # wg (n_s1, n_s0, dim_g)
  # l1
  l1 = int2_ps_om(bw_treated, bw_control, s1, s0, X, copula_type, rho)
  # l1 (n_s1, n_s0, n)
  dim(wg) = c(dim(wg)[1], dim(wg)[2], 1, dim(wg)[3])
  wg = wg[,,rep(1, dim(X)[1]),]
  dim(l1) = c(dim(l1), 1)
  l1 = l1[,,,rep(1, dim(wg)[4])]
  dim(mu1) = c(dim(mu1)[1], 1, dim(mu1)[2], 1)
  dim(mu0) = c(1, dim(mu0)[1], dim(mu0)[2], 1)
  mu1 = mu1[, rep(1, length(s0)),,rep(1, dim(wg)[4])]
  mu0 = mu0[rep(1, length(s1)),,,rep(1, dim(wg)[4])]
  return(list(eta1=wg * l1 * mu1,
              eta0=wg * l1 * mu0))
}

# ----------------------- #
# --- point estimator --- #
# ----------------------- #

#' Point estimator
#' @param Z to be added
#' @return estimated_result to be added
#' @examples point_estimator()
#' @export
point_estimator <- function(Z, X, S, Y, n_divisions=100,
                            copula_type='gaussian',
                            rho=0,
                            weighting_function_vectorized=weighting_function_vectorized,
                            g_function_vectorized=g_function_vectorized) {
  
  # 1. principal score model
  bw_treated <- principal_score_models(S[which(Z == 1)], X[which(Z == 1),])
  bw_control <- principal_score_models(S[which(Z == 0)], X[which(Z == 0),])
  
  P1 <- principal_score_predict_array(bw_treated, S, X) + .Machine$double.eps
  P0 <- principal_score_predict_array(bw_control, S, X) + .Machine$double.eps
  
  # 2. treatment probability model
  treatment_probability_logit <- glm(Z ~ as.matrix(X), family='binomial')
  Tp <- predict(treatment_probability_logit, type='response')
  
  # 3. outcome model
  S_X <- data.frame(cbind(S, data.frame(X)))
  colnames(S_X) <- c('S',  paste(rep('X', dim(X)[2]), 1:dim(X)[2], sep = ""))
  mu1_lm <- lm(Y ~., data=S_X, weights = Z)
  Mu1 <- predict(mu1_lm)
  mu0_lm <- lm(Y ~., data=S_X, weights = 1-Z)
  Mu0 <- predict(mu0_lm)
  
  # integrals
  lower <- min(S) - 3 * std(S)
  upper <- max(S) + 3 * std(S)
  s1 <- linspace(lower, upper, n_divisions)
  s0 <- linspace(lower, upper, n_divisions)
  
  # eif_estimator
  B_int2_integrand <- colMeans(trapz_vector(s0, trapz_vector(s1, B_int2_function(bw_treated, bw_control, s1, s0, Z, Tp, X, S, copula_type, rho))))
  B_ints0_integrand <- colMeans(trapz_vector(s0, B_ints0_function(bw_treated, bw_control, S, s0, Z, Tp, X, copula_type, rho)))
  B_ints1_integrand <- colMeans(trapz_vector(s1, B_ints1_function(bw_treated, bw_control, s1, S, Z, Tp, X, copula_type, rho)))
  vec_int2 <- vec_int2_function(bw_treated, bw_control, mu1_lm, mu0_lm, s1, s0, Z, Tp, X, S, copula_type, rho)
  vec_int_to_s0 <- vec_int_to_s0_function(bw_treated, bw_control, mu0_lm, S, s0, Z, Tp, X, Y, copula_type, rho)
  vec_int_to_s1 <- vec_int_to_s1_function(bw_treated, bw_control, mu1_lm, s1, S, Z, Tp, X, Y, copula_type, rho)
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
  
  # tp_ps_estimator
  B_tp_ps_eta1 <- colMeans(trapz_vector(s0, B_ints0_tp_ps_function(bw_treated, bw_control, S, s0, Z, Tp, X, P1, copula_type, rho)))
  B_tp_ps_eta0 <- colMeans(trapz_vector(s1, B_ints1_tp_ps_function(bw_treated, bw_control, s1, S, Z, Tp, X, P0, copula_type, rho)))
  vec_tp_ps_eta1 <- colMeans(trapz_vector(s0, vec_ints0_tp_ps_function(bw_treated, bw_control, S, s0, Z, Tp, X, Y, P1, copula_type, rho)))
  vec_tp_ps_eta0 <- colMeans(trapz_vector(s1, vec_ints1_tp_ps_function(bw_treated, bw_control, s1, S, Z, Tp, X, Y, P0, copula_type, rho)))
  tp_ps_est_eta1 <- inv(B_tp_ps_eta1) %*% vec_tp_ps_eta1
  tp_ps_est_eta0 <- inv(B_tp_ps_eta0) %*% vec_tp_ps_eta0
  
  # ps_om_estimator
  B_ps_om <- colMeans(trapz_vector(s0, trapz_vector(s1, B_int2_ps_om_function(bw_treated, bw_control, s1, s0, X, copula_type, rho))))
  vec_int2_ps_om <- vec_int2_ps_om_function(bw_treated, bw_control, mu1_lm, mu0_lm, s1, s0, X, copula_type, rho)
  # for eta1
  vec_ps_om_eta1 <- colMeans(trapz_vector(s0, trapz_vector(s1, vec_int2_ps_om$eta1)))
  # for eta0
  vec_ps_om_eta0 <- colMeans(trapz_vector(s1, trapz_vector(s0, vec_int2_ps_om$eta0)))
  ps_om_est_eta1 <- inv(B_ps_om) %*% vec_ps_om_eta1
  ps_om_est_eta0 <- inv(B_ps_om) %*% vec_ps_om_eta0
  
  # eta1_result <- data.frame(cbind(eif_est_eta1, tp_ps_est_eta1, ps_om_est_eta1))
  # colnames(eta1_result) <- c("eif", "tp_ps", "ps_om")
  # rownames(eta1_result) <- c("coef_s1", "coef_s0", "intercept")
  # 
  # eta0_result <- data.frame(cbind(eif_est_eta0, tp_ps_est_eta0, ps_om_est_eta0))
  # colnames(eta0_result) <- c("eif", "tp_ps", "ps_om")
  # rownames(eta0_result) <- c("coef_s1", "coef_s0", "intercept")
  
  # return a data_frame
  # tau_result <- data.frame(cbind(eif_est_eta1 - eif_est_eta0, 
  #                                tp_ps_est_eta1 - tp_ps_est_eta0, 
  #                                ps_om_est_eta1 - ps_om_est_eta0))
  # colnames(tau_result) <- c("eif", "tp_ps", "ps_om")
  # rownames(tau_result) <- c("coef_s1", "coef_s0", "intercept")
  
  # return a list
  tau_result <- as.numeric(cbind(t(eif_est_eta1 - eif_est_eta0),
                                 t(tp_ps_est_eta1 - tp_ps_est_eta0),
                                 t(ps_om_est_eta1 - ps_om_est_eta0)))
  
  return(tau_result)
}

#' Bootstrap variance estimator
#' @param Z to be added
#' @return estimated_result to be added
#' @examples boot()
#' @export
boot <- function(Z, X, S, Y, n_boot=500, n_divisions=100,
                 copula_type='gaussian',
                 rho=0,
                 weighting_function_vectorized=weighting_function_vectorized,
                 g_function_vectorized=g_function_vectorized) {
  
  point_est <- point_estimator(Z, X, S, Y, n_divisions, copula_type, rho, 
                               weighting_function_vectorized, g_function_vectorized=g_function_vectorized)
  
  # nonparametric bootstrap
  n <- length(Z)
  X <- as.matrix(X)
  boot_est <- replicate(n_boot, 
                        {id_boot = sample(1:n, n, replace = TRUE)
                        point_estimator(Z[id_boot], X[id_boot, ], S[id_boot], Y[id_boot], n_divisions, copula_type, rho,
                                        weighting_function_vectorized, g_function_vectorized=g_function_vectorized)})
  
  boot_se <- apply(data.frame(boot_est), 1, sd)
  
  res <- rbind(point_est, boot_se)
  rownames(res) <- c("est", "boot_se")
  return(list(point_est=point_est,
              boot_est=boot_est,
              res=res))
}



# estimation function for one sample
get_results_one_mc_sample <- function(n=500, n_divisions=100) {
  
  # generate data
  n = n
  p = 2
  X = data.frame(matrix(rnorm(n*p, 0, 1), n, p))
  # tilde_X = ((X + 0.25) ^ 2 - 1) / sqrt(2)
  
  # tp model
  Z = rbinom(n, 1, 1/2)
  # ps model
  S1 = 0.6 * X[,1] + 0.8 * rnorm(n)
  S0 = 0.6 * X[,2] + 0.8 * rnorm(n)
  # outcome model
  Y1 = S1 + X[,1] + rnorm(n)
  Y0 = S0 + X[,2] + rnorm(n)
  # observed data
  S = Z * S1 + (1-Z) * S0
  Y = Z * Y1 + (1-Z) * Y0
  data = list(X=X, Z=Z, S=S, Y=Y)
  df = data.frame(data)
  
  if (nrow(X) == 0) {stop("X data is empty!")}
  
  # estimate
  result <- point_estimator(Z, X, S, Y, n_divisions, copula_type='gaussian', rho=0)
  return(result)
}