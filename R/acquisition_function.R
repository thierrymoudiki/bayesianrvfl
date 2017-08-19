acq_function <- function(fit_obj, x_vec, f_best) {
  # fit_obj object prediction at point x_vec
  # --> Need to call fit_rvfl on training with compute_Sigma = FALSE
  fit_obj_pred <- predict_rvfl(fit_obj = fit_obj,
                newx = matrix(x_vec, nrow = 1),
                ci = NULL, graph = FALSE)

  # Predicted mean at point x_vec
  mu_x_vec <- fit_obj_pred$mean

  # Predicted std. dev at point x_vec
  # Avoid std dev. = 0
  Sigma_x_vec <- pmax(fit_obj_pred$sd, 1e-9)

  # Probability of improvement
  gamma_x_vec <- (f_best - mu_x_vec) / Sigma_x_vec

  # Expected improvement
  return(Sigma_x_vec*(gamma_x_vec*pnorm(gamma_x_vec, mean = 0, sd = 1) +
                        dnorm(gamma_x_vec, mean = 0, sd = 1)))
}
