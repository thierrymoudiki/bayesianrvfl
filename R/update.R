# 0 - data -----

# n <- 100 ; p <- 2
# set.seed(123)
# X <- matrix(rnorm(n * p), n, p) # no intercept!
# y <- rnorm(n)
#
# fit_obj <- fit_rvfl(x = X[1:5,], y = y[1:5],
#                     lambda = 0.01, method = "chol",
#                     compute_Sigma = TRUE)
#
# fit_obj2 <- update_params(fit_obj, newx = X[6, ], newy = y[6])
# fit_obj3 <- update_params(fit_obj2, newx = X[7, ], newy = y[7])
# fit_obj4 <- update_params(fit_obj3, newx = X[8, ], newy = y[8])
# fit_obj5 <- update_params(fit_obj4, newx = X[9, ], newy = y[9])
# fit_obj6 <- update_params(fit_obj5, newx = X[10, ], newy = y[10])
#
# bayesianrvfl::predict_rvfl(fit_obj, newx = X[6:10, ])
# bayesianrvfl::predict_rvfl(fit_obj2, newx = X[7:10, ])
# bayesianrvfl::predict_rvfl(fit_obj3, newx = X[8:10, ])
# bayesianrvfl::predict_rvfl(fit_obj4, newx = X[9:10, ])
#
# 1 - plots -----
#
# par(mfrow=c(2, 3))
#
# plot(fit_obj$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj$fitted_values, type = 'l')
# lines(fit_obj$fitted_values, col = 'red')
#
# plot(fit_obj2$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj2$fitted_values, type = 'l')
# lines(fit_obj2$fitted_values, col = 'red')
#
# plot(fit_obj3$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj3$fitted_values, type = 'l')
# lines(fit_obj3$fitted_values, col = 'red')
#
# plot(fit_obj4$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj4$fitted_values, type = 'l')
# lines(fit_obj4$fitted_values, col = 'red')
#
# plot(fit_obj5$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj5$fitted_values, type = 'l')
# lines(fit_obj5$fitted_values, col = 'red')
#
# plot(fit_obj6$y, type = 'l', xlim = c(1, 10),
#      ylim = c(min(y), max(y)))
# lines(fit_obj6$fitted_values, type = 'l')
# lines(fit_obj6$fitted_values, col = 'red')

# 2 - update function -----

update_params <- function(fit_obj,
                          newx,
                          newy,
                          re_clust = TRUE,
                          method = c("direct", "polyak"),
                          alpha = 0.5)
{
  stopifnot(length(newx) >= 1) # newx is a vector
  stopifnot(is.null(dim(newy)) &&
              length(newy) == 1) # newy is a scalar
  if (is.null(fit_obj$Dn))
    stop("for argument 'fit_obj', you should have 'method == chol' in function 'fit_rvfl'")
  if (is.null(fit_obj$Sigma))
    stop("for argument 'fit_obj', you should compute 'Sigma' in 'fit_rvfl'")
  method <- match.arg(method)

  # new information arriving in the system
  mat_newx <- t(newx) # along with newy
  # initial number of covariates
  p <- ncol(fit_obj$x)
  # number of observations at step n
  n <- nrow(fit_obj$x)
  # number of regularization parameters
  nlambda <- length(fit_obj$lambda)

  # parameters at step n
  Dn <- fit_obj$Dn # Cn^{-1}
  ym <- fit_obj$ym
  xm <- as.vector(fit_obj$xm)
  scales <- as.vector(fit_obj$scales)

  # regression parameters for step n + 1
  centered_newy <- newy - ym

  if (is.null(fit_obj$clust_obj))
  {
    augmented_newx <- create_new_predictors(
      x = mat_newx,
      nb_hidden = fit_obj$nb_hidden,
      nodes_sim = fit_obj$nodes_sim,
      activ = fit_obj$activ,
      nn_xm = fit_obj$nn_xm,
      nn_scales = fit_obj$nn_scales
    )$predictors
  } else {
    pred_clusters <- predict(
      fit_obj$clust_obj,
      bayesianrvfl::my_scale(
        mat_newx,
        xm = fit_obj$clusters_scales$means,
        xsd = fit_obj$clusters_scales$sds
      )
    )

    newX_clust <-
      bayesianrvfl::one_hot_encode(pred_clusters$cluster,
                                   fit_obj$n_clusters)

    augmented_newx <-
      create_new_predictors(
        x = cbind(mat_newx, newX_clust),
        nb_hidden = fit_obj$nb_hidden,
        nodes_sim = fit_obj$nodes_sim,
        activ = fit_obj$activ,
        nn_xm = fit_obj$nn_xm,
        nn_scales = fit_obj$nn_scales
      )$predictors
  }

  scaled_augmented_newx <- my_scale(x = augmented_newx, xm = xm,
                                    xsd = scales)

  # ----- update parameters

  ncol_Sigma <- ifelse(nlambda > 1,
                       ncol(fit_obj$Sigma[[1]]),
                       ncol(fit_obj$Sigma))

  if (method == "direct")
  {
    if (nlambda > 1)
    {
      for (i in 1:nlambda) {
        # update factor
        temp <-
          tcrossprod(Dn[[i]], scaled_augmented_newx) # Cn^{-1}%*%t(mat_newx)
        update_factor <-
          Dn[[i]] - tcrossprod(temp) / (1 + drop(scaled_augmented_newx %*% temp))
        resids <-
          centered_newy - scaled_augmented_newx %*% fit_obj$coef

        # update regression coefficients and covariance with update factor
        gradients <- drop(sapply(1:length(resids),
                                 function (i)
                                   scaled_augmented_newx * resids[i]))
        if (is.null(dim(gradients)))
          gradients <- matrix(gradients, ncol = 1,
                              byrow = FALSE)
        fit_obj$coef[, i] <-
          fit_obj$coef[, i] + update_factor %*% gradients[, i]
        fit_obj$Sigma[[i]] <-
          (diag(ncol_Sigma) - update_factor %*% crossprod(scaled_augmented_newx)) %*%
          fit_obj$Sigma[[i]]
      }
    } else {
      # nlambda == 1
      # update factor
      temp <-
        tcrossprod(Dn, scaled_augmented_newx) # Cn^{-1}%*%(t(mat_newx))
      update_factor <-
        Dn - tcrossprod(temp) / (1 + drop(scaled_augmented_newx %*% temp))
      resids <-
        drop(centered_newy - scaled_augmented_newx %*% fit_obj$coef)

      # update regression coefficients and covariance with update factor
      gradients <- as.vector(scaled_augmented_newx * resids)
      fit_obj$coef <- fit_obj$coef + update_factor %*% gradients
      fit_obj$Sigma <-
        (diag(ncol_Sigma) - update_factor %*% crossprod(scaled_augmented_newx)) %*%
        fit_obj$Sigma
    }
  } else {
    # else method == polyak
    # update factor
    update_factor <- n ^ (-alpha)

    if (nlambda > 1)
    {
      for (i in 1:nlambda) {
        # update factor
        temp <-
          tcrossprod(Dn[[i]], scaled_augmented_newx) # Cn^{-1}%*%t(mat_newx)
        resids <-
          centered_newy - scaled_augmented_newx %*% fit_obj$coef

        # update regression coefficients and covariance with update factor
        gradients <- drop(sapply(1:length(resids),
                                 function (i)
                                   scaled_augmented_newx * resids[i]))
        if (is.null(dim(gradients)))
          gradients <- matrix(gradients, ncol = 1,
                              byrow = FALSE)
        fit_obj$coef[, i] <-
          fit_obj$coef[, i] + update_factor %*% gradients[, i]
        fit_obj$Sigma[[i]] <-
          (diag(ncol_Sigma) - update_factor %*% crossprod(scaled_augmented_newx)) %*%
          fit_obj$Sigma[[i]]
      }
    } else {
      # nlambda == 1
      # update factor
      temp <-
        tcrossprod(Dn, scaled_augmented_newx) # Cn^{-1}%*%(t(mat_newx))
      resids <-
        drop(centered_newy - scaled_augmented_newx %*% fit_obj$coef)

      # update regression coefficients and covariance with update factor
      gradients <- as.vector(scaled_augmented_newx * resids)
      fit_obj$coef <- fit_obj$coef + update_factor %*% gradients
      fit_obj$Sigma <-
        (diag(ncol_Sigma) - update_factor %*% crossprod(scaled_augmented_newx)) %*%
        fit_obj$Sigma
    }
  }

  # update response and covariates with new observations
  fit_obj$x <- rbind(as.matrix(fit_obj$x), newx)
  fit_obj$y <- c(fit_obj$y, newy)

  # update ym
  fit_obj$ym <- mean(fit_obj$y)

  # update xm
  if (is.null(fit_obj$clust_obj))
  {
    list_xreg <- create_new_predictors(
      x = fit_obj$x,
      nb_hidden = fit_obj$nb_hidden,
      nodes_sim = fit_obj$nodes_sim,
      activ = fit_obj$activ
    )
  } else {
    if (re_clust)
    {
      fit_obj$clust_obj <- cclust::cclust(fit_obj$x, fit_obj$n_clusters)
    }

    X_clust_obj  <- predict(
      fit_obj$clust_obj,
      bayesianrvfl::my_scale(
        fit_obj$x,
        xm = fit_obj$clusters_scales$means,
        xsd = fit_obj$clusters_scales$sds
      )
    )

    X_clust <- bayesianrvfl::one_hot_encode(X_clust_obj$cluster,
                                            fit_obj$n_clusters)

    list_xreg <-
      create_new_predictors(
        x = cbind(fit_obj$x, X_clust),
        nb_hidden = fit_obj$nb_hidden,
        nodes_sim = fit_obj$nodes_sim,
        activ = fit_obj$activ
      )
  }
  x_scaled <- my_scale(list_xreg$predictors)
  fit_obj$xm <- x_scaled$xm

  # update scales
  fit_obj$scales <- x_scaled$xsd

  # update nn xm and nn scales
  fit_obj$nn_xm <- list_xreg$nn_xm
  fit_obj$nn_scales <- list_xreg$nn_scales

  # update fitted values
  X <- x_scaled$res
  fit_obj$fitted_values <- drop(fit_obj$ym +  X %*% fit_obj$coef)

  # update Dn
  if (method == "direct")
  {
    XTX <- crossprod(X)
    if (nlambda > 1)
    {
      Dn <- vector("list", length = nlambda)
      names(Dn) <- fit_obj$lambda
      Dn <- lapply(1:nlambda, function (i)
        chol2inv(chol(
          XTX + diag(x = fit_obj$lambda[i],
                     nrow = ncol(XTX))
        )))
      fit_obj$Dn <- Dn
    } else {
      Dn <- chol2inv(chol(XTX + diag(
        x = fit_obj$lambda,
        nrow = ncol(XTX)
      )))
      fit_obj$Dn <- Dn
    }
  }

  # serves for variance prediction
  fit_obj$compute_Sigma <- TRUE

  # return a new fit_obj with updated data.
  return (fit_obj)
}

# function for updating with multiple newx arriving
# - fit_obj should specify the method: fit_obj$method == 'direct' or 'pol' (for the next 'newx')
# - fit_obj should specify the alpha: for fit_obj$method == 'pol' (for the next 'newx')
