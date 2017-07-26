predict_rvfl <- function(fit_obj, newx)
{
  newx <- create_new_predictors(x = newx,
                                nb_hidden = fit_obj$nb_hidden,
                                nn_xm = fit_obj$nn_xm,
                                nn_scales = fit_obj$nn_scales)$predictors
  xm <- as.vector(fit_obj$xm)
  scales <- as.vector(fit_obj$scales)

  res <- drop(my_scale(x = as.matrix(newx), xm = xm,
           xsd = scales)%*%as.matrix(fit_obj$coef) + fit_obj$ym)

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
                   sd = sqrt(diag(newx%*%tcrossprod(fit_obj$Sigma, newx) +
                                    lambda*diag(nrow(newx))))))
    } else {
      return (res)
    }
  } else { # nlambda > 1
    if (compute_Sigma == TRUE)
    {
      nSigma <- length(fit_obj$Sigma)
      sd <- foreach(i = 1:nSigma, .combine = cbind)%do%{
        sqrt(diag(newx%*%tcrossprod(fit_obj$Sigma[[i]], newx) +
                    lambda[i]*diag(nrow(newx))))
      }
      colnames(sd) <- lambda
      return (list(mean = res,
                   sd = sd))
    } else {
      return (res)
    }
  }
}
