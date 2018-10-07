# find regularization parameter and number of nodes with GCV
find_lam_nbhidden <- function(x, y, vec_nb_hidden = 1:100, # was 1:100
                              lams = 10^seq(-2, 10, length.out = 100),
                              activ = c("relu", "sigmoid", "tanh"))
{
  activ <- match.arg(activ)

  mat_GCV <- sapply(vec_nb_hidden,
                    function(i) bayesianrvfl::fit_rvfl(x = x, y = y,
                                                       nb_hidden = i, lambda = lams, activ = activ)$GCV)

  best_coords <- which(mat_GCV == min(mat_GCV), arr.ind = TRUE)

  return(list(best_lambda = lams[best_coords[1]],
              best_nb_hidden = vec_nb_hidden[best_coords[2]]))
}

# check if the set of parameter has already been found by the algo
param_is_found <- function(mat, vec)
{
  mat <- as.matrix(mat)
  if (ncol(mat) > 1) # more than 1 parameter
  {
    res <- sum(apply(mat, 1, identical, vec))
    return(ifelse(res >= 1, TRUE, FALSE))
  } else { # 1 parameter
    return(round(vec, 4) %in% round(mat[,1], 4))
  }
}
param_is_found <- compiler::cmpfun(param_is_found)

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

# likelihood minimization
min_loglik <- function(x, y,
                       surrogate_model = c("rvfl", "matern52"),
                       nodes_sim = c("sobol", "halton", "unif"),
                       activ = c("relu", "sigmoid", "tanh",
                                 "leakyrelu","elu", "linear"))
{
  n <- nrow(x)
  stopifnot(n == length(y))
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  surrogate_model <- match.arg(surrogate_model)

  if (surrogate_model == "rvfl")
  {
    OF <- function(xx)
    {
      x_augmented <- create_new_predictors(x = x,
                                           nb_hidden = floor(xx[1]),
                                           nodes_sim = nodes_sim,
                                           activ = activ)$predictors
      Sigma <- tcrossprod(x_augmented, x_augmented) + xx[2]*diag(n)
      res <- try(0.5*(n*log(2*pi) + log(det(Sigma)) +
                        drop(crossprod(y, chol2inv(chol(Sigma)) )%*%y)),
                 silent = TRUE)

      ifelse(class(res) == "try-error", -1e06, res)
    }

    return(msnlminb(objective = OF, nb_iter = 50,
                    lower =  c(1, 0.01),
                    upper = c(100, 1e04)))
  }

  if (surrogate_model == "matern52")
  {
    OF <- function(xx)
    {
      Sigma <- matern52_kxx_cpp(x = x,
                                sigma = xx[1],
                                l = xx[2]) + xx[3]*diag(n)
      res <- try(0.5*(n*log(2*pi) + log(det(Sigma)) +
                        drop(crossprod(y, chol2inv(chol(Sigma)) )%*%y)),
                 silent = TRUE)

      ifelse(class(res) == "try-error", -1e06, res)
    }

    return(msnlminb(objective = OF, nb_iter = 50,
                    lower =  c(1e-04, 1e-04, 1e-04),
                    upper = c(1e05, 1e05, 1e04)))
  }
}

# expand.grid for data frames
expand_grid_df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
expand_grid_df <- compiler::cmpfun(expand_grid_df)

# find integer params in data frame df
is_integer_params <- function(df)
{
  return(sapply(1:ncol(df),
                function (i) is.integer(df[ , i])))
}

# transform data frame df columns to integers when specified in 'is_integer_bool'
transform_cols_integer <- function(df, is_integer_bool)
{
  if (sum(is_integer_bool) > 0){
    # convert specified columns to floor()
    integer_cols <- names(df)[is_integer_bool]
    for (i in 1:length(integer_cols))
    {
      # convert integer columns
      df[ , integer_cols[i]] <- floor(df[ , integer_cols[i]])
    }
  }
  return(df)
}

# initial design of hyperparameters (parameters of the objective)
initial_design <- function(x, y, # x = covariates, y = response
                           # subsets of dataset (x, y), and we have f_approx(x, s = 1) = f(x) (we want to avoid s = 1)
                           method = "svmRadial",
                           metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                           lower_df = data.frame(sigma = 1e-05, C = 1e-05), # hyperparams bounds for 'fit_method'
                           upper_df = data.frame(sigma = 1e06, C = 1e06), # hyperparams bounds for 'fit_method'
                           # method = "xgbTree",
                           # lower_df = data.frame(nrounds = 100L, max_depth = 1L, eta = 0.01,
                           #                       gamma = 0.01, colsample_bytree = 0.5,
                           #                       min_child_weight = 1, subsample = 0.5),
                           # upper_df = data.frame(nrounds = 200L, max_depth = 9L, eta = 0.5,
                           #                       gamma = 0.5, colsample_bytree = 1,
                           #                       min_child_weight = 2, subsample = 1)
                           nb_init = 10, s_max = 0.5, nb_s = 4L,
                           seed = 123, cl = 4)
{
  stopifnot(0 < s_max && s_max <= 1)
  is_integer_lower <- is_integer_params(lower_df)
  is_integer_upper <- is_integer_params(upper_df)
  if(!all.equal(is_integer_lower, is_integer_upper))
  {
    stop("integer upper and lower bounds required")
  }
  dim_xx <- ncol(lower_df)
  stopifnot(dim_xx == ncol(upper_df))

  # initial design for hyperparams search ----------------------------------------------------

  if (s_max < 1)
  {
    s_vec <- cbind.data.frame(s = s_max*randtoolbox::sobol(n = nb_s),
                              join_var = 1)

    rep_1_nb_init <- rep(1, nb_init)
    lower_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(lower_df))
    upper_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(upper_df))

    df_hyper_params <- lower_mat_init + (upper_mat_init - lower_mat_init)*randtoolbox::sobol(n = nb_init, dim = dim_xx,
                                                                                             init = TRUE, scrambling = 2,
                                                                                             seed = seed, normal = FALSE)
    colnames(df_hyper_params) <- names(lower_df)

    # check for integer columns
    if (sum(is_integer_lower) > 0){
      #convert columns to floor()
      integer_cols <- names(lower_df)[is_integer_lower]
      for (i in 1:length(integer_cols))
      {
        # convert integer columns
        df_hyper_params[ , integer_cols[i]] <- floor(df_hyper_params[ , integer_cols[i]])
      }
    }

    parameters <- cbind.data.frame(df_hyper_params,
                                   join_var = 1)

    init_design <- dplyr::select(dplyr::full_join(parameters,
                                                  s_vec, by = "join_var"),
                                 -join_var)
    # hyperparams without 's'
    df_hyper_params_ <- init_design[, 1:(ncol(init_design)-1)]

    # cv error on initial design ----------------------------------------------------

    nrows_init_design <- nrow(init_design)
    allowParallel <- !is.null(cl) && cl > 0
    if (allowParallel)
    {
      cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
      doSNOW::registerDoSNOW(cl_SOCK)

      pb <- txtProgressBar(min = 0, max = nrows_init_design, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      objective <- foreach::foreach(i = 1:nrows_init_design,
                                    .combine = c,
                                    .packages = "doSNOW",
                                    .options.snow = opts,
                                    .export = "f_approx",
                                    .errorhandling = "stop")%dopar%{
                                      hyper_params <- data.frame(df_hyper_params_[i, ])
                                      colnames(hyper_params) <- names(lower_df)
                                      f_approx(x = x, y = y, s = init_design$s[i],
                                               tune_grid = hyper_params,
                                               method = method, metric = metric,
                                               seed = i*100)$results_metric
                                    }
      parallel::stopCluster(cl_SOCK)
    }  else {
      objective <- foreach::foreach(i = 1:nrows_init_design,
                                    .combine = c, .verbose = TRUE,
                                    .export = c("df_hyper_params_", "f_approx", "init_design"),
                                    .errorhandling = "stop")%do%{
                                      hyper_params <- data.frame(df_hyper_params_[i, ])
                                      colnames(hyper_params) <- names(lower_df)
                                      f_approx(x = x, y = y, s = init_design$s[i],
                                               tune_grid = hyper_params,
                                               method = method, metric = metric,
                                               seed = i*100)$results_metric}
    }

    return(list(hyper_params = df_hyper_params,
                s_vec = as.vector(s_vec[,1]),
                init_design = cbind.data.frame(init_design,
                                               cv_error = objective)))

  } else { # s_max == 1

    rep_1_nb_init <- rep(1, nb_init)
    lower_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(lower_df))
    upper_mat_init <- tcrossprod(rep_1_nb_init, as.numeric(upper_df))

    df_hyper_params <- lower_mat_init + (upper_mat_init - lower_mat_init)*randtoolbox::sobol(n = nb_init, dim = dim_xx,
                                                                                             init = TRUE, scrambling = 2,
                                                                                             seed = seed, normal = FALSE)
    colnames(df_hyper_params) <- names(lower_df)
    df_hyper_params <- transform_cols_integer(df = df_hyper_params,
                                              is_integer_bool = is_integer_upper)

    init_design <- cbind.data.frame(df_hyper_params, s = 1)

    # hyperparams without 's'
    df_hyper_params_ <- init_design[, 1:(ncol(init_design)-1)]
    colnames(df_hyper_params_) <- names(lower_df)

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
                                    f_approx(x = x, y = y, s = 1,
                                             tune_grid = df_hyper_params_[i, ],
                                             method = method, metric = metric,
                                             seed = i*100)$results_metric
                                  }
    parallel::stopCluster(cl_SOCK)

    return(list(hyper_params = df_hyper_params,
                s_vec = 1,
                init_design = cbind.data.frame(init_design,
                                               cv_error = objective)))

  }

}

factor_to_matrix <- function(x)
{
  if (!is.factor(x)) x <- as.factor(x)
  levels_x <- unique(x)
  n_levels_x <- length(levels_x)
  n <- length(x)
  res <- matrix(NA, nrow = n,
                ncol = n_levels_x)

    for (j in 1:n_levels_x)
    {
      res[, j] <- as.numeric(x == levels_x[j])
    }

  colnames(res) <- levels(x)
  rownames(res) <- 1:n

  return(res)
}