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

  if(is.matrix(res))
    colnames(res) <- fit_obj$lambda

  return (res)
}
