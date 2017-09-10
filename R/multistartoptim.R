multistartnlminb <- function(objective, lower, upper, nbstartx, ...)
{
  stopifnot(length(lower) == length(upper))
  stopifnot(prod(lower <= upper) == 1)
  nbpar <- length(lower)
  opt_pars <- matrix(0, nrow = nbstartx, ncol = nbpar)
  opt_vals <- rep(0, nbstartx)

  rep_1_nbstartx <- rep(1, nbstartx)
  lower_bound <- tcrossprod(rep_1_nbstartx, lower)
  upper_bound <- tcrossprod(rep_1_nbstartx, upper)
  matrix_pars <- lower_bound + randtoolbox::sobol(n = nbstartx, dim = nbpar)*(upper_bound - lower_bound)
  `%op%` <-  foreach::`%do%`

  i <- NULL
    opt_pars <- foreach::foreach(i = 1:nbstartx, .combine = "rbind")%op%
    {
      res <- try(stats::nlminb(start = matrix_pars[i, ],
                        objective = objective,
                        lower = lower, upper = upper, ...),
                 silent = TRUE)
      if (class(res)[1] == "try-error")
      {
        rep(1000, nbpar)
      } else {
        res$par
      }
    }

  opt_vals <- apply(X = opt_pars, MARGIN = 1, FUN = objective)

  i_opt <- which.min(opt_vals)

  res <- list(as.numeric(opt_pars[i_opt, ]), as.numeric(opt_vals[i_opt]))
  names(res) <- c("par", "value")

  return(res)
}
