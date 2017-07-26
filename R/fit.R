fit_rvfl <- function(x, y, nb_hidden = 5,
                     lambda = 10^seq(from = -5, to = 2, length.out = 100),
                     compute_Sigma = FALSE)
{
  stopifnot(nb_hidden > 0)
  x <- as.matrix(x)
  y <- as.vector(y)
  nlambda <- length(lambda)

  ym <- mean(y)
  centered_y <- y - ym

  list_xreg <- create_new_predictors(x = x, nb_hidden = nb_hidden)

  x_scaled <- my_scale(list_xreg$predictors)
  X <- x_scaled$res

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  nb_di <- length(d)
  div <- d^2 + rep(lambda, rep(nb_di, nlambda))
  a <- drop(d * rhs)/div
  dim(a) <- c(nb_di, nlambda)

    if (nlambda == 1)
    {
      vt <- Xs$vt
      coef <- crossprod(vt, a)
      centered_y_hat <- X %*% coef

        if (compute_Sigma == TRUE)
        {
          rhsX <- crossprod(Xs$u, X)
          aX <- drop(d * rhsX)/div
          Sigma <- diag(ncol(X)) - crossprod(vt, aX)
          rownames(Sigma) <- colnames(X)
          return(list(coef = drop(coef), scales = x_scaled$xsd,
                      Sigma = Sigma, lambda = lambda, ym = ym,
                      xm = x_scaled$xm, nb_hidden = nb_hidden,
                      nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                      fitted_values = drop(ym +  centered_y_hat),
                      compute_Sigma = compute_Sigma))
        } else { #else: compute_Sigma == FALSE && nlambda == 1

          return(list(coef = drop(coef), scales = x_scaled$xsd,
                      lambda = lambda, ym = ym, xm = x_scaled$xm,
                      nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm,
                      nn_scales = list_xreg$nn_scales,
                      fitted_values = drop(ym +  centered_y_hat),
                      compute_Sigma = compute_Sigma))
        }

    } else { #else: nlambda > 1

        coef <- crossprod(Xs$vt, a)
        colnames(coef) <- lambda
        centered_y_hat <- X %*% coef
        fitted_values <- drop(ym +  centered_y_hat)
        colnames(fitted_values) <- lambda

        if (compute_Sigma == TRUE) #compute_Sigma == TRUE && nlambda > 1
        {
            rhsX <- crossprod(Xs$u, X)
            Sigma <- foreach(i = 1:nlambda)%do%{
              div_i <- d^2 + rep(lambda[i], rep(nb_di, 1))
              aX_i <- drop(d * rhsX)/div_i
              Sigma <- diag(ncol(X)) - crossprod(Xs$vt, aX_i)
              rownames(Sigma) <- colnames(X)
              Sigma
            }
            names(Sigma) <- lambda
            return(list(coef = drop(coef), scales = x_scaled$xsd, Sigma = Sigma,
                      lambda = lambda, ym = ym, xm = x_scaled$xm, nb_hidden = nb_hidden,
                      nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                      fitted_values = fitted_values,
                      compute_Sigma = compute_Sigma))
        } else { #else: compute_Sigma == FALSE && length(lambda) == 1

            return(list(coef = drop(coef), scales = x_scaled$xsd,
                      lambda = lambda, ym = ym, xm = x_scaled$xm,
                      nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                      fitted_values = drop(ym +  centered_y_hat),
                      compute_Sigma = compute_Sigma))
        }
    }
}

set.seed(129)

n <- 7 ; p <- 2
X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
fit_rvfl(x = X, y = y)

(fit_obj <- fit_rvfl(x = X, y = y, lambda = c(0.01, 0.05) , compute_Sigma = TRUE))
predict_rvfl(fit_obj, newx = X)

(fit_obj <- fit_rvfl(x = X, y = y, lambda = c(0.01, 0.05) , compute_Sigma = FALSE))
predict_rvfl(fit_obj, newx = X)

(fit_obj <- fit_rvfl(x = X, y = y, lambda = 0.01 , compute_Sigma = TRUE))
predict_rvfl(fit_obj, newx = X)

(fit_obj <- fit_rvfl(x = X, y = y, lambda = 0.01 , compute_Sigma = FALSE))
predict_rvfl(fit_obj, newx = X)
