predict_rvfl <- function(fit_obj, newx, ci = NULL, graph = FALSE)
{

  if (is.vector(newx)) newx <- t(newx)

  newx <- create_new_predictors(x = newx,
                                nb_hidden = fit_obj$nb_hidden,
                                nn_xm = fit_obj$nn_xm,
                                nn_scales = fit_obj$nn_scales,
                                activ = fit_obj$activ,
                                nodes_sim = fit_obj$nodes_sim)$predictors
  xm <- as.vector(fit_obj$xm)
  scales <- as.vector(fit_obj$scales)
  scaled_newx <- my_scale(x = as.matrix(newx), xm = xm,
                          xsd = scales)
  n <- nrow(scaled_newx)

  res <- drop(scaled_newx%*%as.matrix(fit_obj$coef) + fit_obj$ym)

  lambda <- fit_obj$lambda
  nlambda <- length(lambda)
  compute_Sigma <- fit_obj$compute_Sigma

  if(is.matrix(res))
    colnames(res) <- lambda

  if(nlambda == 1)
  {
    if (compute_Sigma == TRUE)
    {
      return (list(mean = res,
                   sd = sqrt(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma, scaled_newx) +
                                    lambda*diag(n)))))
    } else {
      return (res)
    }
  } else { # nlambda > 1
    if (compute_Sigma == TRUE)
    {
      i <- NULL
      `%op%` <-  foreach::`%do%`
      sd <- foreach::foreach(i = 1:nlambda, .combine = cbind)%op%{
        sqrt(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma[[i]], scaled_newx) +
                    lambda[i]*diag(n)))
      }
      colnames(sd) <- lambda
      return (list(mean = res,
                   sd = sd))
    } else {
      return (res)
    }
  }
}

predict_matern52 <- function(fit_obj, newx, ci = NULL)
{
   if (is.vector(newx)) newx <- t(newx)

    y <- fit_obj$y
    xm <- fit_obj$xm
    ym <- fit_obj$ym
    xsd <- fit_obj$xsd
    X <- fit_obj$scaled_x
    sigma <- fit_obj$sigma
    l <- fit_obj$l
    lambda_krls <- fit_obj$lambda_krls
    K <- fit_obj$K
    mat_coefs <- fit_obj$mat_coefs
    compute_Sigma <- fit_obj$compute_Sigma
    scaled_newx <- my_scale(newx, xm = xm, xsd = xsd)
    inv_method <- fit_obj$inv_method
    n_newx <- nrow(newx)

    if (is.vector(scaled_newx)){
      K_star <- matern52_kxy_cpp(x = X, y = scaled_newx,
                                 sigma = sigma, l = l)
    } else {
      K_star <- sapply(1:n_newx, function (i) matern52_kxy_cpp(x = X, y = scaled_newx[i, ],
                                 sigma = sigma, l = l))
    }

    preds <- drop(crossprod(K_star, mat_coefs)) + ym

    if (compute_Sigma == TRUE)
    {
      K_star2 <- matern52_kxx_cpp(x = scaled_newx,
                                  sigma = sigma, l = l)
      shrinked_K <- switch(
        inv_method,
        "chol" = chol2inv(chol(K + lambda_krls * diag(nrow(K)))),
        "ginv" = bayesianrvfl::my_ginv(K + lambda_krls * diag(nrow(K))))
      Sigma <- sqrt(diag(K_star2 - crossprod(K_star, shrinked_K)%*%K_star))
      return (list(mean = preds,
                   sd = Sigma))
    } else {
      return(preds)
    }
}
