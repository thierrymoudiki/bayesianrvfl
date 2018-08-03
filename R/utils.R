find_next_param_by_es <- function(x)
{
  n <- 10000
  x <- matrix(x, nrow = 1)
  pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                                                activ = activation_function,
                                                                nb_hidden = best_nb_hidden,
                                                                lambda = best_lam,
                                                                compute_Sigma = TRUE),
                                         newx = x)

  corr_mat <- diag(pred_obj$sd)

  norm_density_value <- mvtnorm::dmvnorm(x, mean = pred_obj$mean,
                              sigma = corr_mat,
                              log = FALSE)

  norm_pdf_value <- mvtnorm::pmvnorm(lower = rep(-Inf, ncol(x)),
                            upper = x,
                            mean = pred_obj$mean,
                            corr = corr_mat)

  p <- n*((1 - norm_pdf_value)^(n-1))*norm_density_value

  q <- dunif(x, min = lower, max = upper,
             log = FALSE)

  return (sum(p*log(p/q)))
}
find_next_param_by_es <- compiler::cmpfun(find_next_param_by_es)
