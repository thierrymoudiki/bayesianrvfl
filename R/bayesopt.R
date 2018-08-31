
# 0 - utils ---------------------------------------------------------------

  # find regularization parameter and number of nodes with GCV
  find_lam_nbhidden <- function(x, y, vec_nb_hidden = 1:100, # was 1:100
                                lams = 10^seq(-2, 10, length.out = 100),
                                activ = c("relu", "sigmoid", "tanh"))
  {
    activ <- match.arg(activ)

    mat_GCV <- sapply(vec_nb_hidden,
                      function(i) bayesianrvfl::fit_rvfl(x = x, y = y,
                                           nb_hidden = i, lambda = lams, activ = activ)$GCV)
    #colnames(mat_GCV) <- paste0('nb_hidden=', vec_nb_hidden)

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
  braninsc <- compiler::cmpfun(braninsc)

  # Hartmann 6
  hart6sc <- function(xx)
  {
    ##########################################################################
    #
    # HARTMANN 6-DIMENSIONAL FUNCTION, RESCALED
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2, x3, x4, x5, x6)
    #
    ##########################################################################

    alpha <- c(1.0, 1.2, 3.0, 3.2)
    A <- c(10, 3, 17, 3.5, 1.7, 8,
           0.05, 10, 17, 0.1, 8, 14,
           3, 3.5, 1.7, 10, 17, 8,
           17, 8, 0.05, 10, 0.1, 14)
    A <- matrix(A, 4, 6, byrow=TRUE)
    P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                     2329, 4135, 8307, 3736, 1004, 9991,
                     2348, 1451, 3522, 2883, 3047, 6650,
                     4047, 8828, 8732, 5743, 1091, 381)
    P <- matrix(P, 4, 6, byrow=TRUE)

    xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
    inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
    outer <- sum(alpha * exp(-inner))

    y <- -outer
    return(y)
  }
  hart6sc <- compiler::cmpfun(hart6sc)

  # Alpine 01
  alpine01 <- function(x)
  {
    sum(abs(x * sin(x) + 0.1 * x))
  }
  alpine01 <- compiler::cmpfun(alpine01)

  # 1D function
  test1Dfunction <- function(x)
  {
    -(exp(-(x - 2)^2) + exp(-(x - 6)^2/10) + 1/(x^2 + 1))
  }
  test1Dfunction <- compiler::cmpfun(test1Dfunction)

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

# 1 - optimization functions ---------------------------------------------------------------

# 1 - 1 bayesian optimization ---------------------------------------------------------------

  bayes_opt <- function(objective, lower, upper,
                             type_acq = c("ei", "ucb"),
                             kappa = 1.96,
                             nb_init = 10, nb_iter = 25,
                             method = c("standard", "direct_online",
                                        "polyak_online"), # '*_online' available for rvfl only
                             surrogate_model = c("rvfl", "matern52", "rf"),
                             optim_surr = c("GCV", "loglik"),
                             activation_function = c("relu", "tanh"),
                             type_optim = c("nlminb", "DEoptim",
                                            "msnlminb"),
                             seed = 123,
                             verbose = TRUE,
                             show_progress = TRUE,
                             record_points = FALSE, ...)
  {
    OF <- function(y) objective(y, ...)
    nb_is_found <- 0
    dim_xx <- length(lower)
    stopifnot(dim_xx == length(upper))
    type_acq <- match.arg(type_acq)
    type_optim <- match.arg(type_optim)
    method <- match.arg(method)
    activation_function <- match.arg(activation_function)
    optim_surr <- match.arg(optim_surr)
    surrogate_model <- match.arg(surrogate_model)
    if (surrogate_model == "matern52") optim_surr <- "loglik" # to me MODIFIED

    rep_1_nb_init <- rep(1, nb_init)
    lower_mat_init <- tcrossprod(rep_1_nb_init, lower)
    upper_mat_init <- tcrossprod(rep_1_nb_init, upper)

    set.seed(seed)
    parameters <- lower_mat_init + (upper_mat_init - lower_mat_init)*matrix(runif(nb_init*dim_xx),
                                                                            nrow = nb_init, ncol = dim_xx)
    scores <- apply(parameters, 1, OF)

    parameters <- as.matrix(parameters[!is.na(scores),])
    scores <- scores[!is.na(scores)]

    if (verbose == TRUE)
    {
      cat("initial design", "\n")
      print(cbind.data.frame(parameters, scores))
      cat("\n")
    }

    if (optim_surr == "GCV")
    {
      best_params <- bayesianrvfl::find_lam_nbhidden(parameters, scores)
      best_lam <- best_params$best_lambda
      best_nb_hidden <- best_params$best_nb_hidden
    }

    if (optim_surr == "loglik")
    {
      if (surrogate_model == "rvfl")
      {
        best_params <- bayesianrvfl::min_loglik(x = parameters, y = scores,
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

    # method == "standard" ----------------------------------------------------

    if (method == "standard")
    {
      if (surrogate_model == "rvfl")
      {
        find_next_param_by_ei <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                                                        activ = activation_function,
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
        find_next_param_by_ei <- compiler::cmpfun(find_next_param_by_ei)

        find_next_param_by_ucb <- function(x)
        {
          x <- matrix(x, nrow = 1)
          pred_obj <- bayesianrvfl::predict_rvfl(bayesianrvfl::fit_rvfl(x = parameters, y = scores,
                                                                        activ = activation_function,
                                                                        nb_hidden = best_nb_hidden,
                                                                        lambda = best_lam,
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

          current_score <- OF(next_param)
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

          current_score <- OF(next_param)
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

                                OF(searched_points[i, ], ...)
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
                                 OF(searched_points[i, ], ...)
                              })
      snow::stopCluster(cl_SOCK)
    }

    index_opt <- which.min(res)

    return(list(par = searched_points[index_opt, ],
                objective = res[index_opt]))
  }

  # Example with mlrBO -----------------------------------------------------------------

  # fn <- makeSingleObjectiveFunction(
  #   name = "hart6sc",
  #   fn = function(xx) {
  #     alpha <- c(1.0, 1.2, 3.0, 3.2)
  #     A <- c(10, 3, 17, 3.5, 1.7, 8,
  #            0.05, 10, 17, 0.1, 8, 14,
  #            3, 3.5, 1.7, 10, 17, 8,
  #            17, 8, 0.05, 10, 0.1, 14)
  #     A <- matrix(A, 4, 6, byrow=TRUE)
  #     P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
  #                      2329, 4135, 8307, 3736, 1004, 9991,
  #                      2348, 1451, 3522, 2883, 3047, 6650,
  #                      4047, 8828, 8732, 5743, 1091, 381)
  #     P <- matrix(P, 4, 6, byrow=TRUE)
  #
  #     xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
  #     inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  #     outer <- sum(alpha * exp(-inner))
  #
  #     y <- -outer
  #     return(y)
  #   },
  #   par.set = makeParamSet(
  #     makeNumericParam("x1", lower = 0, upper = 1),
  #     makeNumericParam("x2", lower = 0, upper = 1)
  #   )
  # )
  #
  # # Create initial random Latin Hypercube Design of 10 points
  # library(lhs) # for randomLHS
  # library(DiceKriging)
  # des <- generateDesign(n = 25L, getParamSet(fn), fun = randomLHS)
  #
  # # Specify kriging model with standard error estimation
  # surrogate <- makeLearner("regr.km", predict.type = "se",
  #                         covtype = "matern5_2")
  #
  # # Set general controls
  # ctrl <- makeMBOControl()
  # ctrl <- setMBOControlTermination(ctrl, iters = 200L)
  # ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
  # # Start optimization
  # set.seed(1)
  # obj_mbo <- mbo(fn, des, surrogate, ctrl,
  #                show.info = getOption("mlrMBO.show.info", FALSE))
  #
