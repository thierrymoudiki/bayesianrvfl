fit_rvfl <- function(x, y, nb_hidden = 5,
                     lambda = 10^seq(from = -5, to = 2, length.out = 100))
{
  stopifnot(nb_hidden > 0)
  x <- as.matrix(x)
  y <- as.vector(y)

  ym <- mean(y)
  centered_y <- y - ym

  list_xreg <- create_new_predictors(x = x, nb_hidden = nb_hidden)

  x_scaled <- my_scale(list_xreg$predictors)
  X <- x_scaled$res

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  k <- length(lambda)
  dx <- length(d)
  div <- d^2 + rep(lambda, rep(dx, k))
  a <- drop(d * rhs)/div
  dim(a) <- c(dx, k)
  coef <- crossprod(Xs$vt, a)

  #dimnames(coef) <- list(names(Xscale), format(lambda))
  centered_y_hat <- X %*% coef
  #GCV <- colSums((centered_y -  centered_y_hat)^2)/(nrow(X) - colSums(matrix(d^2/div,
  #                                                       dx)))^2

    return(list(coef = drop(coef), scales = x_scaled$xsd,
                lambda = lambda, ym = ym, xm = x_scaled$xm, nb_hidden = nb_hidden,
                nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                #GCV = GCV,
                fitted_values = drop(ym +  centered_y_hat)))

}
