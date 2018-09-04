# 1 - functions ---------------------------------------------------------------

# G finds an approximation of the objective (cv error) for stratified
# subsets of dataset (x, y), and we have G(x, s = 1) = f(x) (we want to avoid s = 1)
# x: covariates (typically, parameters of the objective function)
# y: response (objective function value at x)
# tune_grid: named hyperparams vector for method `method` (tuneGrid)
# s: fraction of data in x to be trained (tune_grid x s)
# method: fitting model to be optimized (from the caret list)
# number: number of folds for the cross-validation
# repeats: number of repeats for the cross-validation procedure
# seed: random seed for selection the fraction of data to be trained
# ...: additional parameters to be passed to caret::train or caret::trainControl
G <- function(x, y,
              s = 0.5,
              tune_grid = data.frame(sigma = 1e03, C = 1e-03),
              method = "svmRadial",
              metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
              number = 5L, repeats = 1L, seed = 123,
              allowParallel = FALSE, ...)
{
  stopifnot(0 < s && s <= 1)
  stopifnot(nrow(x) == length(y))

    if (s < 1)
    {
      # fraction (s) of the dataset to be trained
      set.seed(seed)
      training_index <- drop(caret::createDataPartition(y = 1:nrow(x),
                                                        p = s,
                                                        list = FALSE))
      # subset the dataset
      x_train <- x[training_index, ]
      colnames(x_train) <- 1:ncol(x_train)
      y_train <- y[training_index]

    } else { # s == 1
      # subset the dataset
      x_train <- x
      colnames(x_train) <- 1:ncol(x_train)
      y_train <- y
    }

  # caret train parameters
  trControl <- caret::trainControl(method = "repeatedcv",
                                   number = number,
                                   repeats = repeats,
                                   allowParallel = allowParallel,
                                   ...)

  set.seed(123) # doesn't change
  res <- suppressWarnings(caret::train(x = x_train, y = y_train,
                                       method = method, trControl = trControl,
                                       tuneGrid = tune_grid, metric = metric,
                                       verbose = FALSE, ...))

    # return
    if (is.factor(y)){
      return(list(results_details = res,
                  results_metric = res$results$Accuracy))
    } else {
      return(list(results_details = res,
                  results_metric = res$results$RMSE))
    }
}

# initial design of hyperparameters (parameters of the objective)
initial_design <- function(G, # G(x, y, s) finds an approximation of the objective (cv error) for stratified
                           # subsets of dataset (x, y), and we have G(x, s = 1) = f(x) (we want to avoid s = 1)
                           x, y, # x = covariates, y = response
                           method = "svmRadial",
                           metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                           lower_df, upper_df, # bounds for hyperparams search, (x, y, s) cv error
                           nb_init = 10, s_max = 0.5, nb_s = 4L,
                           seed = 123, cl = 4)
{
  stopifnot(0 < s_max && s_max <= 1)
  dim_xx <- ncol(lower_df)
  stopifnot(dim_xx == ncol(upper_df))

  # initial design for hyperparams search ----------------------------------------------------

  s_vec <- cbind.data.frame(s = s_max*randtoolbox::sobol(n = nb_s),
                            join_var = 1)

  rep_1_nb_init <- rep(1, nb_init)

  lower_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(lower_df))
  upper_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(upper_df))

  set.seed(seed)
  df_hyper_params <- lower_mat_init + (upper_mat_init - lower_mat_init)*matrix(runif(nb_init*dim_xx),
                                                                  nrow = nb_init,
                                                                  ncol = dim_xx)
  colnames(df_hyper_params) <- names(lower_df)

  parameters <- cbind.data.frame(df_hyper_params,
                                 join_var = 1)

  init_design <- dplyr::select(dplyr::full_join(parameters,
                                                s_vec, by = "join_var"),
                               -join_var)
  # hyperparams without 's'
  df_hyper_params_ <- init_design[, 1:(ncol(init_design)-1)]

  # cv error on initial design ----------------------------------------------------

  allowParallel <- !is.null(cl) && cl > 0
  if (allowParallel)
  {
    cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
    doSNOW::registerDoSNOW(cl_SOCK)
    `%op%` <-  foreach::`%dopar%`
  }  else {
    `%op%` <-  foreach::`%do%`
  }

  nrows_init_design <- nrow(init_design)
  pb <- txtProgressBar(min = 0, max = nrows_init_design, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
      objective <- foreach::foreach(i = 1:nrows_init_design,
                                    .combine = c, .packages = "doSNOW",
                                    .options.snow = opts)%op%{
        G(x = x, y = y, s = init_design$s[i],
          tune_grid = df_hyper_params_[i, ],
          method = method, metric = metric,
           seed = i*100,
          allowParallel = TRUE)$results_metric
      }
  parallel::stopCluster(cl_SOCK)

  return(list(hyper_params = df_hyper_params,
              s_vec = dplyr::select(s_vec, s)$s,
              init_design = cbind.data.frame(init_design,
                                           cv_error = objective)))
}

# expand.grid for data frames
expand_grid_df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
expand_grid_df <- compiler::cmpfun(expand_grid_df)


# 2 - tests ---------------------------------------------------------------
#library(bayesianrvfl)

#  # on random dataset
#   n <- 500 ; p <- 5
#   X <- matrix(rnorm(n * p), n, p) # no intercept!
#   y <- rnorm(n)
#     (res <- G(x = X, y = y))
#
#  # on iris dataset
#   data(iris)
#   TrainData <- iris[,1:4]
#   TrainClasses <- iris[,5]
#   (res <- G(x = TrainData, y = TrainClasses))
#
#  # initial design
#   (df <- initial_design(G = G, x = TrainData, y = TrainClasses,
#                       lower_df = data.frame(sigma = 1e-05, C = 1e-05),
#                       upper_df = data.frame(sigma = 1e05, C = 1e05),
#                       nb_init = 5, nb_s = 4, seed = 12))
#    dplyr::filter(df$init_design, s == 0.375)
#    dplyr::filter(df$init_design, s == 0.25)
#
# # use other activ
#  (best_params <- bayesianrvfl::find_lam_nbhidden(x = df$init_design[, -ncol(df$init_design)],
#                    y = df$init_design[, "cv_error"], activ = "tanh"))
#
# # use other activ
# fit_obj <- bayesianrvfl::fit_rvfl(x = df$init_design[, -ncol(df$init_design)],
#                        y = df$init_design[, "cv_error"],
#                        nb_hidden = best_params$best_nb_hidden,
#                        lambda = best_params$best_lambda,
#                        activ = "tanh",
#                        compute_Sigma = TRUE)
#
# OF <- function(x) bayesianrvfl::predict_rvfl(fit_obj, newx = cbind(matrix(x, nrow = 1), 1))$mean
#
# # # projected cv error on whole dataset (high variance...)
#  index_best <- which.min(apply(df$hyper_params, 1, OF))
#  (x_best <- df$hyper_params[index_best, ])
#  (f_best <- OF(x_best))

bayes_opt_largedata <- function(x, y, G,
                                df_init_design = NULL,
                                lower_df = data.frame(sigma = 1e-05, C = 1e-05), # hyperparams for 'fit_method'
                                upper_df = data.frame(sigma = 1e06, C = 1e06), # hyperparams for 'fit_method'
                                nb_init = 10, s_max = 0.5, nb_s = 4L, # params for constructing the initial design
                                seed = 123, # for initial design?
                                fit_method = "svmRadial", # caret training methods on (x, y)
                                metric = ifelse(is.factor(y), "Accuracy", "RMSE"), # metric for cv error (on caret training methods)
                                type_acq = c("ei", "ucb"), kappa = 1.96, # acquistion function
                                method = c("standard", "direct_online", "polyak_online"), # '*_online' available for rvfl only
                                surrogate_model = c("rvfl", "matern52", "rf"),
                                optim_surr = c("GCV", "loglik"),
                                activation_function = c("relu", "tanh", "sigmoid"),
                                type_optim = c("nlminb", "DEoptim", "msnlminb"),
                                nb_iter = 25,
                                verbose = TRUE,
                                show_progress = TRUE,
                                record_points = FALSE, cl = 4, ...)
{
  nb_is_found <- 0
  dim_xx <- length(lower)
  stopifnot(dim_xx == length(upper))
  type_acq <- match.arg(type_acq)
  type_optim <- match.arg(type_optim)
  method <- match.arg(method)
  activation_function <- match.arg(activation_function)
  optim_surr <- match.arg(optim_surr)
  surrogate_model <- match.arg(surrogate_model)
  lower <- as.numeric(lower_df)
  upper <- as.numeric(upper_df)
  next_param_cv_error <- rep(0, nb_s)

  ## --- define initial design
  if(is.null(df_init_design)) {
    if (verbose == TRUE)
    {
      cat("\n", "----- define initial design...", "\n")
    }
    df_init_design <- initial_design(G = G, x = x, y = y,
                                     method = fit_method, metric = metric,
                                     lower_df = lower_df, upper_df = upper_df,
                                     nb_init = nb_init, s_max = s_max, nb_s = nb_s,
                                     seed = seed, cl = cl)
    }

  if (verbose == TRUE)
  {
    cat("\n", "----- initial design", "\n")
    print(df_init_design$init_design)
    cat("\n")
  }

    # initial design for hyperparams
    x_design <- df_init_design$init_design[, -ncol(df_init_design$init_design)]

    # cv error associated to hyperparams on the initial design
    y_design <- df_init_design$init_design[, "cv_error"]

  ## --- end define initial design

  ## --- define objective function
  if (optim_surr == "GCV")
  {
    # optimal params for the surrogate
    best_params <- bayesianrvfl::find_lam_nbhidden(x = x_design,
                                                   y = y_design,
                                                   activ = activation_function)
    best_lam <- best_params$best_lambda
    best_nb_hidden <- best_params$best_nb_hidden
  }

  if (optim_surr == "loglik")
  {
    if (surrogate_model == "rvfl")
    {
      # optimal params for the surrogate
      best_params <- bayesianrvfl::min_loglik(x = x_design,
                                              y = y_design,
                                              nodes_sim = "sobol",
                                              activ = activation_function)
      best_lam <- best_params$par[2]
      best_nb_hidden <- floor(best_params$par[1])

      if (verbose == TRUE)
      {
        cat("\n")
        if (optim_surr == "GCV") cat("----- GCV parameters", "\n")
        if (optim_surr == "loglik") cat("----- loglik parameters", "\n")
        cat("\n")
        cat("selected regularization parameter", "\n")
        print(best_lam)
        cat("\n")
        cat("selected number of hidden nodes", "\n")
        print(best_nb_hidden)
        cat("\n")
      }
    }

    if (surrogate_model == "matern52")
    {
      # optimal params for the surrogate
      best_params <- bayesianrvfl::min_loglik(x = parameters, y = scores,
                                              surrogate_model = "matern52")
      best_sigma <- best_params$par[1]
      best_l <- best_params$par[2]
      best_lam <- best_params$par[3]

      if (verbose == TRUE)
      {
        if (optim_surr == "loglik") cat("----- loglik parameters", "\n")
        cat("\n")
        cat("selected sigma", "\n")
        print(best_sigma)
        cat("\n")
        cat("selected lengthscale", "\n")
        print(best_l)
        cat("\n")
        cat("selected regularization parameter", "\n")
        print(best_lam)
        cat("\n")
      }
    }
  }

    # objective function = approx of f = G(x, s = 1)
    fit_obj_surr <- bayesianrvfl::fit_rvfl(x = x_design,
                                      y = y_design,
                                      nb_hidden = best_nb_hidden,
                                      lambda = best_lam,
                                      activ = activation_function,
                                      compute_Sigma = TRUE)
    OF <- function(x, ...) bayesianrvfl::predict_rvfl(fit_obj_surr,
                                                      newx = cbind(matrix(x, nrow = 1), 1))$mean

    ## --- end define objective function
    parameters <- df_init_design$hyper_params
    scores <- apply(parameters, 1, OF) # scores calculated with G_hat(x, s = 1)
    parameters <- as.matrix(parameters[!is.na(scores),])
    scores <- scores[!is.na(scores)]

    if (verbose == TRUE)
    {
      cat("scores", "\n")
      print(cbind.data.frame(parameters, scores))
      cat("\n")
    }

  # method == "standard" ----------------------------------------------------

  if (method == "standard")
  {
      if (surrogate_model == "rvfl")
      {
        find_next_param_by_ei <- function(x)
        {
          x_data <- as.matrix(expand_grid_df(parameters,
                                   data.frame(df_init_design$s_vec)))

          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = x_data,
                                                                        y = rep(scores,
                                                                                nrow(x_data)/length(scores)),
                                                                        activ = activation_function,
                                                                        nb_hidden = best_nb_hidden,
                                                                        lambda = best_lam,
                                                                        compute_Sigma = TRUE),
                                                 newx = cbind(matrix(x, nrow = 1), 1))

          mu_hat <- pred_obj$mean
          sigma_hat <- pred_obj$sd
          gamma_hat <- (min(scores) - mu_hat)/sigma_hat
          res <- -sigma_hat*(gamma_hat*pnorm(gamma_hat) + dnorm(gamma_hat))
          return (ifelse(is.na(res), 100, res))
        }
        find_next_param_by_ei <- compiler::cmpfun(find_next_param_by_ei)

        find_next_param_by_ucb <- function(x)
        {
          x_data <- as.matrix(expand_grid_df(parameters,
                                             data.frame(df$s_vec)))

          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = x_data,
                                                                        y = rep(scores,
                                                                                nrow(x_data)/length(scores)),
                                                                        activ = activation_function,
                                                                        nb_hidden = best_nb_hidden,
                                                                        lambda = best_lam,
                                                                        compute_Sigma = TRUE),
                                                 newx = cbind(matrix(x, nrow = 1), 1))

          return (-(pred_obj$mean - kappa*pred_obj$sd))
        }
        find_next_param_by_ucb <- compiler::cmpfun(find_next_param_by_ucb)

        find_next_param <- switch(type_acq,
                                  "ei" = find_next_param_by_ei,
                                  "ucb" = find_next_param_by_ucb)

        if (verbose == FALSE && show_progress == TRUE)
        {
          pb <- txtProgressBar(min = 1, max = nb_iter, style = 3)
        }

        # main loop
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
                                                                                               parallelType = 0,
                                                                                               itermax = 25))$optim$bestmem)
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

          #current_score <- OF(next_param) # consider also replacing this with evaluations of G for s < s_max
          df_next_param <- data.frame(matrix(next_param, nrow = 1))
          colnames(df_next_param) <- colnames(df_init_design$hyper_params)
          # n_cores <- parallel::detectCores()
          # cl <- snow::makeCluster(n_cores, type="SOCK")
          # doSNOW::registerDoSNOW(cl)

            for(j in 1:length(df_init_design$s_vec)){
              next_param_cv_error[j] <- G(x = x, y = y, tune_grid = df_next_param,
                                    s = df_init_design$s_vec[j], seed = iter*100,
                                    allowParallel = FALSE)$results_metric
            }

          # parallel::stopCluster(cl)
          next_param_design <- as.matrix(expand_grid_df(df_next_param,
                                         data.frame(df_init_design$s_vec)))
          colnames(next_param_design) <- colnames(x_design)

          # update x_design
          x_design <- rbind.data.frame(x_design, next_param_design)
          y_design <- c(y_design, next_param_cv_error)

          #objective function = approx of f = G(x, s = 1)
          OF <- function(x, ...) bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = x_design,
                                                                                   y = y_design,
                                                                                   nb_hidden = best_nb_hidden,
                                                                                   lambda = best_lam,
                                                                                   activ = activation_function,
                                                                                   compute_Sigma = TRUE),
                                                            newx = cbind(matrix(x, nrow = 1), 1))$mean
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

          if (verbose == FALSE && show_progress == TRUE) setTxtProgressBar(pb, iter)
        }
        if (verbose == FALSE && show_progress == TRUE) close(pb)

      } # end: if (surrogate_model == "rvfl")

      if (surrogate_model == "matern52")
      {
        find_next_param_by_ei <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_matern52(bayesianrvfl::fit_matern52(x = parameters,
                                                                                y = scores,
                                                                                sigma = best_sigma,
                                                                                l = best_l,
                                                                                lambda_krls = best_lam,
                                                                                compute_Sigma = TRUE),
                                                     newx = x)
          mu_hat <- pred_obj$mean
          sigma_hat <- pred_obj$sd
          gamma_hat <- (min(scores) - mu_hat)/sigma_hat
          res <- -sigma_hat*(gamma_hat*pnorm(gamma_hat) + dnorm(gamma_hat))
          return (ifelse(is.na(res), 100, res))
        }
        find_next_param_by_ei <- compiler::cmpfun(find_next_param_by_ei)

        find_next_param_by_ucb <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_matern52(bayesianrvfl::fit_matern52(x = parameters, y = scores,
                                                                                sigma = best_sigma,
                                                                                l = best_l,
                                                                                lambda_krls = best_lam,
                                                                                compute_Sigma = TRUE),
                                                     newx = x)
          return (-(pred_obj$mean - kappa*pred_obj$sd))
        }
        find_next_param_by_ucb <- compiler::cmpfun(find_next_param_by_ucb)

        find_next_param <- switch(type_acq,
                                  "ei" = find_next_param_by_ei,
                                  "ucb" = find_next_param_by_ucb)

        if (verbose == FALSE && show_progress == TRUE)
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

          # NOT CORRECT?
          current_score <- OF(next_param) # consider also replacing this with evaluation of G for s < s_max
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

          if (verbose == FALSE && show_progress == TRUE) setTxtProgressBar(pb, iter)
        }
        if (verbose == FALSE && show_progress == TRUE) close(pb)

      } # end: if (surrogate_model == "matern52")
  }

  # method == "direct_online" ----------------------------------------------------

  if (method == "direct_online")
  {
    if (surrogate_model == "rvfl")
    {

      fit_obj <- bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                        nb_hidden = best_nb_hidden,
                                        activ = activation_function,
                                        lambda = best_lam,
                                        method = "chol",
                                        compute_Sigma = TRUE)

      find_next_param_by_ei <- function(x)
      {
        x <- matrix(x, nrow = 1)
        pred_obj <- bayesianrvfl::predict_rvfl(fit_obj,
                                               newx = x)
        mu_hat <- pred_obj$mean
        sigma_hat <- pred_obj$sd
        gamma_hat <- (min(scores) - mu_hat)/sigma_hat
        res <- -sigma_hat*(gamma_hat*pnorm(gamma_hat) + dnorm(gamma_hat))
        return (ifelse(is.na(res), 100, res))
      }
      find_next_param_by_ei <- compiler::cmpfun(find_next_param_by_ei)

      find_next_param_by_ucb <- function(x)
      {
        x <- matrix(x, nrow = 1)
        pred_obj <- bayesianrvfl::predict_rvfl(fit_obj,
                                               newx = x)
        return (-(pred_obj$mean - kappa*pred_obj$sd))
      }
      find_next_param_by_ucb <- compiler::cmpfun(find_next_param_by_ucb)

      find_next_param <- switch(type_acq,
                                "ei" = find_next_param_by_ei,
                                "ucb" = find_next_param_by_ucb)

      if (verbose == FALSE && show_progress == TRUE)
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

        # NOT CORRECT?
        current_score <- OF(next_param) # consider also replacing this with evaluation of G for s < s_max
        if (verbose == TRUE)
        {
          cat("score", "\n")
          print(current_score)
          cat("\n")
        }

        parameters <- rbind(parameters, next_param)
        scores <- c(scores, current_score)
        fit_obj <- bayesianrvfl::update_params(fit_obj = fit_obj,
                                               newx = next_param, newy = current_score,
                                               method = "direct")

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

        if (verbose == FALSE && show_progress == TRUE) setTxtProgressBar(pb, iter)
      }
      if (verbose == FALSE && show_progress == TRUE) close(pb)
    }

    if (surrogate_model == "matern52")
    {stop("not implemented")}
  }

  # method == "polyak_online" ----------------------------------------------------

  if (method == "polyak_online")
  {
    if (surrogate_model == "rvfl")
    {

      fit_obj <- bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                        nb_hidden = best_nb_hidden,
                                        activ = activation_function,
                                        lambda = best_lam,
                                        method = "chol",
                                        compute_Sigma = TRUE)
      # with averaged coeffs
      fit_obj2 <- fit_obj
      mat_coefs <- matrix(0, ncol = nb_iter + 1,
                          nrow = length(fit_obj$coef))
      colnames(mat_coefs) <- 1:(nb_iter + 1)
      mat_coefs[ , 1] <- fit_obj$coef
      #current_coefs <- fit_obj$coef

      find_next_param_by_ei <- function(x)
      {
        x <- matrix(x, nrow = 1)
        pred_obj <- bayesianrvfl::predict_rvfl(fit_obj2,
                                               newx = x)
        mu_hat <- pred_obj$mean
        sigma_hat <- pred_obj$sd
        gamma_hat <- (min(scores) - mu_hat)/sigma_hat
        res <- -sigma_hat*(gamma_hat*pnorm(gamma_hat) + dnorm(gamma_hat))
        return (ifelse(is.na(res), 100, res))
      }
      find_next_param_by_ei <- compiler::cmpfun(find_next_param_by_ei)

      find_next_param_by_ucb <- function(x)
      {
        x <- matrix(x, nrow = 1)
        pred_obj <- bayesianrvfl::predict_rvfl(fit_obj2,
                                               newx = x)
        return (-(pred_obj$mean - kappa*pred_obj$sd))
      }
      find_next_param_by_ucb <- compiler::cmpfun(find_next_param_by_ucb)

      find_next_param <- switch(type_acq,
                                "ei" = find_next_param_by_ei,
                                "ucb" = find_next_param_by_ucb)

      if (verbose == FALSE && show_progress == TRUE)
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
                                                                                             parallelType = 0,
                                                                                             itermax = 25))$optim$bestmem)
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

        # NOT CORRECT?
        current_score <- OF(next_param) # consider also replacing this with evaluation of G for s < s_max
        if (verbose == TRUE)
        {
          cat("score", "\n")
          print(current_score)
          cat("\n")
        }

        parameters <- rbind(parameters, next_param)
        scores <- c(scores, current_score)

        fit_obj <- bayesianrvfl::update_params(fit_obj = fit_obj,
                                               newx = next_param,
                                               newy = current_score,
                                               method = "polyak")

        fit_obj2 <- fit_obj
        mat_coefs[ , (iter + 1)] <- fit_obj$coef
        #current_coefs <- current_coefs + (fit_obj$coef - current_coefs)/(iter + 1)
        #print(drop(current_coefs))
        #print(rowSums(mat_coefs)/(iter + 1))
        fit_obj2$coef <- rowSums(mat_coefs)/(iter + 1)
        #fit_obj2$coef <- current_coefs

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

        if (verbose == FALSE && show_progress == TRUE) setTxtProgressBar(pb, iter)
      }
      if (verbose == FALSE && show_progress == TRUE) close(pb)
    }

    if (surrogate_model == "matern52")
    {stop("not implemented")}
  }

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
