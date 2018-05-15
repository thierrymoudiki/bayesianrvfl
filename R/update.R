update_params <- function(fit_obj, newx, newy,
                          method = c("direct", "polyak"),
                          alpha = 0.5)
{
  stopifnot(is.null(nrow(newx)) && length(newx) > 1) # newx is a vector
  stopifnot(is.null(dim(newy)) && length(newy) == 1) # newy is a scalar
  if(is.null(fit_obj$Dn)) stop("for argument 'fit_obj', you should have 'method == ginv' in function 'fit_rvfl'")
  if(is.null(fit_obj$Sigma)) stop("for argument 'fit_obj', you should compute 'Sigma' in 'fit_rvfl'")
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
      augmented_newx <- create_new_predictors(x = mat_newx,
                                    nb_hidden = fit_obj$nb_hidden,
                                    nodes_sim = fit_obj$nodes_sim,
                                    activ = fit_obj$activ,
                                    nn_xm = fit_obj$nn_xm,
                                    nn_scales = fit_obj$nn_scales)$predictors
      scaled_augmented_newx <- my_scale(x = augmented_newx, xm = xm,
                              xsd = scales)

    # ----- update parameters

        ncol_Sigma <- ncol(fit_obj$Sigma[[1]])

        if (method == "direct")
        {

          for(i in 1:nlambda){

            # update factor
              temp <- tcrossprod(Dn[[i]], scaled_augmented_newx) # Cn^{-1}%*%(t(mat_newx))
              update_factor <- Dn[[i]] - tcrossprod(temp)/(1 + drop(scaled_augmented_newx%*%temp))
              resids <- centered_newy - scaled_augmented_newx%*%fit_obj$coef

            # update regression coefficients and covariance with update factor
               gradients <- drop(sapply(1:length(resids),
                                        function (i) scaled_augmented_newx*resids[i]))
               fit_obj$coef[, i] <- fit_obj$coef[, i] + update_factor%*%gradients[, i]
               fit_obj$Sigma[[i]] <- fit_obj$Sigma[[i]] + update_factor%*%crossprod(scaled_augmented_newx)%*%(fit_obj$Sigma[[i]] -
                                                                                        2*diag(ncol_Sigma))
          }

        } else {

          # update factor
            update_factor <- n^(-alpha)

          for(i in 1:nlambda){

            # update regression coefficient with update factor
              resids <- centered_newy - scaled_augmented_newx%*%fit_obj$coef
              gradients <- drop(sapply(1:length(resids),
                                       function (i) scaled_augmented_newx*resids[i]))
              fit_obj$coef[, i] <- fit_obj$coef[, i] + update_factor*gradients[, i]

            # update regression coefficients and covariance with update factor
              fit_obj$Sigma[[i]] <- fit_obj$Sigma[[i]] + update_factor*crossprod(scaled_augmented_newx)%*%(fit_obj$Sigma[[i]] -
                                                                                                             2*diag(ncol_Sigma))
          }

        }

          # update response and covariates with new observations
            fit_obj$x <- rbind(as.matrix(fit_obj$x), newx)
            fit_obj$y <- c(fit_obj$y, newy)

          # update ym
            fit_obj$ym <- mean(fit_obj$y)

          # update xm
            list_xreg <- create_new_predictors(x = fit_obj$x,
                                               nb_hidden = fit_obj$nb_hidden,
                                               nodes_sim = fit_obj$nodes_sim,
                                               activ = fit_obj$activ)
            x_scaled <- my_scale(list_xreg$predictors)
            fit_obj$xm <- x_scaled$xm

          # update scales
            fit_obj$scales <- x_scaled$xsd

          # update nn xm and nn scales
            fit_obj$nn_xm <- list_xreg$nn_xm
            fit_obj$nn_scales <- list_xreg$nn_scales

          # update fitted values
            X <- x_scaled$res
            fit_obj$fitted_values <- drop(fit_obj$ym +  X%*%fit_obj$coef)

          # update Dn
            if (method == "direct")
            {
              XTX <- crossprod(X)
              Dn <- vector("list", length = nlambda)
              names(Dn) <- fit_obj$lambda
              Dn <- lapply(1:nlambda, function (i)
                MASS::ginv(XTX + diag(x = fit_obj$lambda[i],
                                               nrow = ncol(XTX))))
              fit_obj$Dn <- Dn
            }

  # return a new fit_obj with updated data.
  return (fit_obj)
}
