fit_rvfl <- function(x, y, nb_hidden = 5, nodes_sim = c("sobol", "halton", "unif"),
                     activ = c("relu", "sigmoid", "tanh",
                               "leakyrelu", "elu", "linear"),
                     a = 0.01, lambda = 10^seq(from = -5, to = 2, length.out = 100),
                     seed = 1)
{
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  if (nb_hidden > 0)
  {
    list_xreg <- create_new_predictors(x = x, nb_hidden = nb_hidden,
                                       method = nodes_sim, activ = activ,
                                       a = a, seed = seed)
    X <- list_xreg$predictors
  } else {
    X <- x
  }

  ym <- mean(y)
  centered_y <- y - ym

  x_scaled <- my_scale(X)
  X <- x_scaled$res
  Xscale <- x_scaled$xsd
  Xm <- x_scaled$xm

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)

  d <- Xs$d
  lscoef <- crossprod(Xs$vt, rhs/d)

  lsfit <- X %*% lscoef
  resids <- centered_y - lsfit

  k <- length(lambda)
  dx <- length(d)
  div <- d^2 + rep(lambda, rep(dx, k))
  a <- drop(d * rhs)/div
  dim(a) <- c(dx, k)
  coef <- crossprod(Xs$vt, a)

  dimnames(coef) <- list(names(Xscale), format(lambda))
  GCV <- colSums((centered_y - X %*% coef)^2)/(nrow(X) - colSums(matrix(d^2/div,
                                                         dx)))^2

  if (nb_hidden > 0)
  {
    return(list(coef = drop(coef), scales = Xscale,
                lambda = lambda, ym = ym, xm = Xm, GCV = GCV,
                list_xreg = list_xreg))
  } else {
    return(list(coef = drop(coef), scales = Xscale,
                lambda = lambda, ym = ym, xm = Xm, GCV = GCV))
  }

}
