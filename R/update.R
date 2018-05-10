update_params <- function(fit_obj, newx, newy)
{
  stopifnot(is.vector(newx))
  newx <- t(newx)
  Dn <- fit_obj$Dn
  beta_n <- fit_obj$coef
  ym <- fit_obj$ym

  newy <- newy - ym
  newx <- create_new_predictors(x = as.vector(newx),
                                nb_hidden = fit_obj$nb_hidden,
                                nn_xm = fit_obj$nn_xm,
                                nn_scales = fit_obj$nn_scales)$predictors
  xm <- as.vector(fit_obj$xm)
  scales <- as.vector(fit_obj$scales)
  scaled_newx <- my_scale(x = newx, xm = xm,
                          xsd = scales)

}
