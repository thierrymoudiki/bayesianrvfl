predict_rvfl <- function(fit_obj, newx, ci = NULL, graph = FALSE)
{
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
