# predicting from an rvfl ----
predict_rvfl <- function(fit_obj, newx, ci = NULL, graph = FALSE)
{

  if (is.vector(newx)) newx <- t(newx)

  newx <- create_new_predictors(x = newx,
                                nb_hidden = fit_obj$nb_hidden,
                                nn_xm = fit_obj$nn_xm,
                                nn_scales = fit_obj$nn_scales,
                                activ = fit_obj$activ,
                                nodes_sim = fit_obj$nodes_sim)$predictors
  xm <- as.vector(fit_obj$xm)
  scales <- as.vector(fit_obj$scales)
  scaled_newx <- my_scale(x = as.matrix(newx), xm = xm,
                          xsd = scales)
  n <- nrow(scaled_newx)

  res <- drop(scaled_newx%*%as.matrix(fit_obj$coef) + fit_obj$ym)

  lambda <- fit_obj$lambda
  nlambda <- length(lambda)
  compute_Sigma <- fit_obj$compute_Sigma

  if(is.matrix(res))
    colnames(res) <- lambda

  if(nlambda == 1)
  {
    if (compute_Sigma == TRUE)
    {
      cat("in predict: sigma_hat", "\n")
      print(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma, scaled_newx) + lambda*diag(n)))
      cat("\n")
      #stop("in predict")

      return (list(mean = res,
                   sd = sqrt(pmax(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma, scaled_newx) +
                                    lambda*diag(n)), 1e-02))
                   ))
    } else {
      return (res)
    }
  } else { # nlambda > 1
    if (compute_Sigma == TRUE)
    {
      i <- NULL
      `%op%` <-  foreach::`%do%`
      sd <- foreach::foreach(i = 1:nlambda, .combine = cbind)%op%{

        cat("in predict: sigma_hat", "\n")
        print(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma[[i]], scaled_newx) + lambda[i]*diag(n)))
        cat("\n")

        sqrt(pmax(diag(scaled_newx%*%tcrossprod(fit_obj$Sigma[[i]], scaled_newx) +
                    lambda[i]*diag(n)), 1e-02))
      }
      colnames(sd) <- lambda
      return (list(mean = res,
                   sd = sd))
    } else {
      return (res)
    }
  }
}

# predicting from Matérn 5/2 model ----
predict_matern52 <- function(fit_obj, newx, ci = NULL)
{
   if (is.vector(newx)) newx <- t(newx)

    y <- fit_obj$y
    xm <- fit_obj$xm
    ym <- fit_obj$ym
    xsd <- fit_obj$xsd
    X <- fit_obj$scaled_x
    sigma <- fit_obj$sigma
    l <- fit_obj$l
    lambda_krls <- fit_obj$lambda_krls
    K <- fit_obj$K
    mat_coefs <- fit_obj$mat_coefs
    compute_Sigma <- fit_obj$compute_Sigma
    scaled_newx <- my_scale(newx, xm = xm, xsd = xsd)
    inv_method <- fit_obj$inv_method
    n_newx <- nrow(newx)

    if (is.vector(scaled_newx)){
      K_star <- matern52_kxy_cpp(x = X, y = scaled_newx,
                                 sigma = sigma, l = l)
    } else {
      K_star <- sapply(1:n_newx, function (i) matern52_kxy_cpp(x = X, y = scaled_newx[i, ],
                                 sigma = sigma, l = l))
    }

    preds <- drop(crossprod(K_star, mat_coefs)) + ym

    if (compute_Sigma == TRUE)
    {
      K_star2 <- matern52_kxx_cpp(x = scaled_newx,
                                  sigma = sigma, l = l)
      shrinked_K <- switch(
        inv_method,
        "chol" = chol2inv(chol(K + lambda_krls * diag(nrow(K)))),
        "ginv" = bayesianrvfl::my_ginv(K + lambda_krls * diag(nrow(K))))
      Sigma <- sqrt(diag(K_star2 - crossprod(K_star, shrinked_K)%*%K_star))
      return (list(mean = preds,
                   sd = Sigma))
    } else {
      return(preds)
    }
}

# predicting from elastic net ----
predict_glmnet <- function(fit_obj, newx, s = 0.1)
{
  if (is.vector(newx)) newx <- t(newx)

  if (fit_obj$compute_Sigma == FALSE)
  {
    newx <- create_new_predictors(x = newx,
                                  nb_hidden = fit_obj$nb_hidden,
                                  nn_xm = fit_obj$nn_xm,
                                  nn_scales = fit_obj$nn_scales,
                                  activ = fit_obj$activ,
                                  nodes_sim = fit_obj$nodes_sim)$predictors
    xm <- as.vector(fit_obj$xm)
    scales <- as.vector(fit_obj$scales)
    scaled_newx <- my_scale(x = as.matrix(newx), xm = xm,
                            xsd = scales)

      return(predict(fit_obj$fit_obj, newx = scaled_newx, s = s))

    } else {

      n_resamples <- length(fit_obj$fit_obj_list)

      boot_preds <- foreach::foreach(i = 1:n_resamples, .combine = cbind)%do%{
        newxpreds <- create_new_predictors(x = newx,
                                      nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$nodes_sim,
                                      activ = fit_obj$activ
                                      )$predictors

        scaled_newx <- my_scale(x = as.matrix(newxpreds),
                                xm = as.vector(fit_obj$xm[[i]]),
                                xsd = as.vector(fit_obj$scales[[i]]))
         predict(fit_obj$fit_obj_list[[i]],
                                    newx = scaled_newx, s = s)
      }

      return(list(mean = rowMeans(boot_preds),
                  sd = apply(boot_preds, 1, sd)))
    }

}

par(mfrow = c(2, 2))

fit_obj <- bayesianrvfl::fit_glmnet(x = mtcars[,-1], y = mtcars[,1],
                                    nb_hidden = 50, compute_Sigma = TRUE)
preds <- predict_glmnet(fit_obj, newx = mtcars[,-1], s = 0.05)

plot(preds$mean, ylim = c(0, 40), type = 'l')
lines(mtcars[,1], col = "blue")
lines(preds$mean - 1.96*preds$sd, col = "red")
lines(preds$mean + 1.96*preds$sd, col = "red")



X <- longley[,-1]
y <- longley[,1]
fit_obj <- bayesianrvfl::fit_glmnet(x = X, y = y,
                                     nb_hidden = 100, compute_Sigma = TRUE)
preds <- predict_glmnet(fit_obj, newx = X, s = 0.05)

 plot(preds$mean, ylim = c(50, 150), type = 'l')
 lines(y, col = "blue")
 lines(preds$mean - 1.96*preds$sd, col = "red")
 lines(preds$mean + 1.96*preds$sd, col = "red")

set.seed(123)
n <- 25 ; p <- 2
X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
fit_obj <- bayesianrvfl::fit_glmnet(x = X, y = y,
                                     nb_hidden = 200, compute_Sigma = TRUE)
preds <- predict_glmnet(fit_obj, newx = X, s = 0.1)

 plot(preds$mean, ylim = c(-2, 3), type = 'l')
 lines(y, col = "blue")
 lines(preds$mean - 1.96*preds$sd, col = "red")
 lines(preds$mean + 1.96*preds$sd, col = "red")


set.seed(225)
n <- 25 ; p <- 2
X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
fit_obj <- bayesianrvfl::fit_glmnet(x = X, y = y,
                                     nb_hidden = 200, compute_Sigma = TRUE)
preds <- predict_glmnet(fit_obj, newx = X, s = 0.1)

 plot(preds$mean, ylim = c(-2, 3), type = 'l')
 lines(y, col = "blue")
 lines(preds$mean - 1.96*preds$sd, col = "red")
 lines(preds$mean + 1.96*preds$sd, col = "red")
