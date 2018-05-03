
# 0 - utils ---------------------------------------------------------------

  # find regularization parameter and number of nodes with GCV
  find_lam_nbhidden <- function(x, y, vec_nb_hidden = 1:100,
                                lams = 10^seq(-20, 1, length.out = 100))
  {
    mat_GCV <- sapply(vec_nb_hidden,
                      function(i) fit_rvfl(x = x, y = y,
                                           nb_hidden = i, lambda = lams)$GCV)
    #colnames(mat_GCV) <- paste0('nb_hidden=', vec_nb_hidden)

    best_coords <- which(mat_GCV == min(mat_GCV), arr.ind = TRUE)

    return(list(best_lambda = lams[best_coords[1]],
                best_nb_hidden = vec_nb_hidden[best_coords[2]]))
  }

  # check if the set of parameter has already been found by the algo
  param_is_found <- function(mat, vec)
  {
    res <- sum(sapply(1:nrow(mat),
                      function(i) identical(mat[i, ], vec)))
    return(ifelse(res >= 1, TRUE, FALSE))
  }

  # scaled branin function for testing
  braninsc <- function(xx)
  {
    x1_bar <- 15*xx[1] - 5
    x2_bar <- 15*xx[2]

    term1 <- (x2_bar - (5.1/(4*pi^2)) * x1_bar^2 + (5/pi)*x1_bar - 6)^2
    term2 <- 10*(1-1/(8*pi))*cos(x1_bar)
    z <- (term1 + term2 - 44.81) / 51.95
    return(z)
  }

# 1 - optimization functions ---------------------------------------------------------------

# 1 - 1 bayesian optimization ---------------------------------------------------------------

  rvfl_bayes_opt <- function(objective, lower, upper,
                             type_acq = c("ei", "ucb"),
                             kappa = 1.96,
                             nb_init = 10, nb_iter = 25,
                             type_optim = c("nlminb", "DEoptim", "msnlminb"),
                             seed = 123,
                             verbose = TRUE,
                             record_points = FALSE, ...)
  {
    OF <- function(y) objective(y, ...)
    nb_is_found <- 0
    dim_xx <- length(lower)
    stopifnot(dim_xx == length(upper))
    type_acq <- match.arg(type_acq)
    type_optim <- match.arg(type_optim)

    rep_1_nb_init <- rep(1, nb_init)
    lower_mat_init <- tcrossprod(rep_1_nb_init, lower)
    upper_mat_init <- tcrossprod(rep_1_nb_init, upper)

    set.seed(seed)
    parameters <- lower_mat_init + (upper_mat_init - lower_mat_init)*matrix(runif(nb_init*dim_xx),
                                                                            nrow = nb_init, ncol = dim_xx)
    scores <- apply(parameters, 1, OF)

    parameters <- parameters[!is.na(scores),]
    scores <- scores[!is.na(scores)]

    best_params <- bayesianrvfl::find_lam_nbhidden(parameters, scores)
    best_lam <- best_params$best_lambda
    best_nb_hidden <- best_params$best_nb_hidden

    if (verbose == TRUE)
    {
      cat("----- GCV parameters", "\n")
      cat("\n")
      cat("selected regularization parameter", "\n")
      print(best_lam)
      cat("\n")
      cat("selected number of hidden nodes", "\n")
      print(best_nb_hidden)
      cat("\n")
    }

        find_next_param_by_ei <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                                                        nb_hidden = best_nb_hidden,
                                                                        lambda = best_lam,
                                                                        compute_Sigma = TRUE),
                                                 newx = x)
          mu_hat <- pred_obj$mean
          sigma_hat <- pred_obj$sd
          gamma_hat <- (min(scores) - mu_hat)/sigma_hat
          res <- -sigma_hat*(gamma_hat*pnorm(gamma_hat) + dnorm(gamma_hat))
          return (ifelse(is.na(res), 100, res))
        }

        find_next_param_by_ucb <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                                                        nb_hidden = best_nb_hidden,
                                                                        lambda = best_lam,
                                                                        compute_Sigma = TRUE),
                                                 newx = x)
          return (-(pred_obj$mean - kappa*pred_obj$sd))
        }

    find_next_param <- switch(type_acq,
                              "ei" = find_next_param_by_ei,
                              "ucb" = find_next_param_by_ucb)

    if (verbose == FALSE)
    {
      pb <- txtProgressBar(min = 1, max = nb_iter, style = 3)
    }

    for (iter in 1:nb_iter)
    {
      if (verbose == TRUE)
      {
        cat("\n")
        cat("----- iteration #", iter, "\n")
        cat("\n")
      }

      if (type_optim == "nlminb")
      {
        set.seed(iter + 1)
        next_param <- suppressWarnings(stats::nlminb(start = lower + (upper-lower)*runif(length(lower)),
                                                     objective = find_next_param,
                                                     lower = lower, upper = upper)$par)
      }

      if (type_optim == "msnlminb")
      {
       next_param <- suppressWarnings(bayesianrvfl::msnlminb(objective = find_next_param,
                                                    lower = lower,
                                                    upper = upper,
                                               nb_iter = 10)$par)
      }

      if (type_optim == "DEoptim")
      {
       next_param <- suppressWarnings(DEoptim::DEoptim(fn = find_next_param,
                                                    lower = lower,
                                                    upper = upper,
                                                    control = DEoptim::DEoptim.control(trace = FALSE,
                                                    parallelType = 0, itermax = 25))$optim$bestmem)
      }

      if (param_is_found(parameters, next_param) == TRUE)
      {
        nb_is_found <- nb_is_found + 1
        set.seed(iter + 1)
        next_param <- lower + (upper - lower)*runif(dim_xx)
      }

      if (verbose == TRUE)
      {
        cat("next_param", "\n")
        print(next_param)
        cat("\n")
      }

      current_score <- OF(next_param)
      if (verbose == TRUE)
      {
        cat("score", "\n")
        print(current_score)
        cat("\n")
      }

      parameters <- rbind(parameters, next_param)
      scores <- c(scores, current_score)
      if (verbose == TRUE)
      {
        index_min <- which.min(scores)
        best_param <- parameters[index_min,]
        cat("current best param", "\n")
        print(best_param)
        cat("\n")
        cat("current best score", "\n")
        print(scores[index_min])
        cat("\n")
      }

      if (verbose == FALSE) setTxtProgressBar(pb, iter)
    }
    if (verbose == FALSE) close(pb)

    # final best params
    index_min <- which.min(scores)
    best_param <- parameters[index_min,]

    if (record_points == FALSE){
      return(list(index_min = index_min,
                  nb_is_found = nb_is_found,
                  best_param = best_param,
                  best_value = scores[index_min]))
    } else {

      points_found <- cbind.data.frame(parameters, scores)
      n_params <- ncol(parameters)
      colnames(points_found) <- c(paste0("param", 1:n_params), "score")

      return(list(index_min = index_min,
                  nb_is_found = nb_is_found,
                  best_param = best_param,
                  best_value = scores[index_min],
                  points_found = cbind(parameters, scores)))
    }

  }

# 1 - 2 multistart nlminb ---------------------------------------------------------------

  msnlminb <- function(objective, nb_iter = 100, lower, upper, cl = NULL,...)
  {
    OF <- function(y) objective(y, ...)
    rep_1_nb_iter <- rep(1, nb_iter)
    lower_mat <- tcrossprod(rep_1_nb_iter, lower)
    upper_mat <- tcrossprod(rep_1_nb_iter, upper)

    starting_points <- lower_mat + (upper_mat -
                                      lower_mat)*randtoolbox::sobol(n = nb_iter,
                                                                    dim = length(lower))
    nb_iter <- nrow(starting_points)

    if (is.null(cl)){
      res <- lapply(1:nb_iter, function (i)  nlminb(start = starting_points[i, ],
                                                    objective = OF,
                                                    lower = lower, upper = upper, ...))

    } else {
      cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
      doSNOW::registerDoSNOW(cl_SOCK)

      pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      i <- j <- NULL
      res <- foreach::foreach(i = 1:nb_iter, .packages = "doSNOW",
                              .options.snow = opts,
                              .verbose = FALSE,
                              .errorhandling = "stop")%dopar%{
                                opt <- nlminb(start = starting_points[i, ], objective = OF,
                                              lower = lower, upper = upper, ...)
                                opt
                              }
      stopCluster(cl_SOCK)
    }

    index_opt <- which.min(sapply(1:nb_iter,
                                  function (i) res[[i]]$objective))

    return(res[[index_opt]])
  }

# 1 - 3 random search ---------------------------------------------------------------

  random_search_opt <- function(objective, nb_iter = 100, lower, upper,
                                sim = c("sobol", "unif"), seed = 123,
                                cl = NULL, ...)
  {
    OF <- function(y) objective(y, ...)
    rep_1_nb_iter <- rep(1, nb_iter)
    lower_mat <- tcrossprod(rep_1_nb_iter, lower)
    upper_mat <- tcrossprod(rep_1_nb_iter, upper)
    sim <- match.arg(sim)

      if (sim == "sobol")
      {
        searched_points <- lower_mat + (upper_mat -
                                          lower_mat)*randtoolbox::sobol(n = nb_iter,
                                                                        dim = length(lower))
      }

      if (sim == "unif")
      {
        set.seed(seed)
        searched_points <- lower_mat + (upper_mat -
                                          lower_mat)*matrix(runif(nb_iter*length(lower)),
                                                            nrow = nrow(lower_mat))
      }

    if (is.null(cl)){
      pb <- txtProgressBar(min = 1, max = nb_iter, style = 3)

      `%op%` <-  foreach::`%do%`

      res <- foreach::foreach(i = 1:nb_iter, .combine = c,
                              .verbose = FALSE,
                              .errorhandling = "stop")%op%{

                                setTxtProgressBar(pb, i)

                                OF(searched_points[i, ])
                              }

      close(pb)

    } else {
      cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
      doSNOW::registerDoSNOW(cl_SOCK)

      pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      i <- j <- NULL
      res <- suppressWarnings(foreach::foreach(i = 1:nb_iter, .combine = c,
                              .packages = "doSNOW",
                              .options.snow = opts,
                              .export = ...,
                              .verbose = FALSE,
                              .errorhandling = "stop")%dopar%{
                                 OF(searched_points[i, ])
                              })
      snow::stopCluster(cl_SOCK)
    }

    index_opt <- which.min(res)

    return(list(par = searched_points[index_opt, ],
                objective = res[index_opt]))
  }

