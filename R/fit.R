# fitting and rvfl
fit_rvfl <- function(x, y, nb_hidden = 5,
                     nodes_sim = c("sobol", "halton", "unif"),
                     activ = c("relu", "sigmoid", "tanh",
                               "leakyrelu", "elu", "linear"),
                     lambda = 10^seq(from = -5, to = 2, length.out = 100),
                     method = c("svd", "chol"),
                     compute_Sigma = FALSE)
{
  if (!is.null(dim(y))) stop("'y' must be a vector") # otherwise y - ym is not working

    if (!is.factor(y))
    {
        ## regression
        stopifnot(nb_hidden > 0)
        x <- as.matrix(x)
        y <- as.vector(y)
        nlambda <- length(lambda)
        method <- match.arg(method)

        ym <- mean(y)
        centered_y <- y - ym

        nodes_sim <- match.arg(nodes_sim)
        activ <- match.arg(activ)
        list_xreg <- create_new_predictors(x = x, nb_hidden = nb_hidden,
                                           nodes_sim = nodes_sim,
                                           activ = activ)

        x_scaled <- my_scale(list_xreg$predictors)
        X <- x_scaled$res

        if (method == "svd")
        {
          # inspired from MASS::lm.ridge
          Xs <- La.svd(X)
          rhs <- crossprod(Xs$u, centered_y)
          d <- Xs$d
          nb_di <- length(d)
          div <- d^2 + rep(lambda, rep(nb_di, nlambda))
          a <- drop(d * rhs)/div
          dim(a) <- c(nb_di, nlambda)
          n <- nrow(X)

          if (nlambda == 1)
          {
            vt <- Xs$vt
            coef <- crossprod(vt, a)
            centered_y_hat <- X %*% coef
            GCV <- colSums((centered_y - centered_y_hat)^2)/(n - colSums(matrix(d^2/div,
                                                                                nb_di)))^2

            if (compute_Sigma == TRUE)
            {
              rhsX <- crossprod(Xs$u, X)
              aX <- drop(d * rhsX)/div
              Sigma <- diag(ncol(X)) - crossprod(vt, aX)
              rownames(Sigma) <- colnames(X)
              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          Sigma = Sigma, lambda = lambda, ym = ym,
                          xm = x_scaled$xm, nb_hidden = nb_hidden,
                          nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  centered_y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            } else { #else: compute_Sigma == FALSE && nlambda == 1
              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          lambda = lambda, ym = ym, xm = x_scaled$xm,
                          nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm,
                          nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  centered_y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            }
          } else { #else: nlambda > 1
            coef <- crossprod(Xs$vt, a)
            colnames(coef) <- lambda
            centered_y_hat <- X %*% coef
            fitted_values <- drop(ym +  centered_y_hat)
            colnames(fitted_values) <- lambda
            GCV <- colSums((centered_y - centered_y_hat)^2)/(n - colSums(matrix(d^2/div,
                                                                                nb_di)))^2

            if (compute_Sigma == TRUE) #compute_Sigma == TRUE && nlambda > 1
            {
              rhsX <- crossprod(Xs$u, X)
              `%op%` <-  foreach::`%do%`
              i <- NULL
              Sigma <- foreach::foreach(i = 1:nlambda)%op%{
                div_i <- d^2 + rep(lambda[i], rep(nb_di, 1))
                aX_i <- drop(d * rhsX)/div_i
                Sigma <- diag(ncol(X)) - crossprod(Xs$vt, aX_i)
                rownames(Sigma) <- colnames(X)
                Sigma
              }
              names(Sigma) <- lambda
              return(list(coef = drop(coef), scales = x_scaled$xsd, Sigma = Sigma,
                          lambda = lambda, ym = ym, xm = x_scaled$xm, nb_hidden = nb_hidden,
                          nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = fitted_values,
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            } else { #else: compute_Sigma == FALSE && length(lambda) == 1

              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          lambda = lambda, ym = ym, xm = x_scaled$xm,
                          nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  centered_y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            }
          }
        }

        if (method == "chol")
        {
          XTX <- crossprod(X)

          if (nlambda > 1)
          {
            Dn <- vector("list", length = nlambda)
            names(Dn) <- lambda
            Sigma <- vector("list", length = nlambda)
            names(Sigma) <- lambda
            coef <- foreach::foreach(i = 1:nlambda, .combine = cbind)%do%{
              Dn[[i]] <- chol2inv(chol(XTX + diag(x = lambda[i],
                                                  nrow = nrow(XTX)))) # Cn^{-1}
              Sigma[[i]] <- diag(ncol(X)) - Dn[[i]]%*%XTX # Sigma_n
              Dn[[i]]%*%crossprod(X, centered_y) # beta_n
            }
            colnames(coef) <- lambda
            rownames(coef) <- colnames(x_scaled$res)
          } else {
            Dn <- chol2inv(chol(XTX + diag(x = lambda,
                                           nrow = nrow(XTX)))) # Cn^{-1}
            Sigma <- diag(ncol(X)) - Dn%*%XTX # Sigma_n
            coef <- Dn%*%crossprod(X, centered_y) # beta_n
          }

          return(list(coef = coef, Dn = Dn, Sigma = Sigma,
                      scales = x_scaled$xsd, lambda = lambda,
                      ym = ym, xm = x_scaled$xm,
                      nb_hidden = nb_hidden,
                      nodes_sim = nodes_sim,
                      activ = activ,
                      nn_xm = list_xreg$nn_xm,
                      nn_scales = list_xreg$nn_scales,
                      fitted_values = drop(ym + X %*% coef),
                      compute_Sigma = compute_Sigma,
                      x = x, y = y))
        }

    } else { # if is.factor(y)

        ## classification
        stopifnot(nb_hidden > 0)
        x <- as.matrix(x)
        y <- factor_to_matrix(y)
        nlambda <- length(lambda)

        nodes_sim <- match.arg(nodes_sim)
        activ <- match.arg(activ)
        list_xreg <- create_new_predictors(x = x, nb_hidden = nb_hidden,
                                           nodes_sim = nodes_sim,
                                           activ = activ)

        x_scaled <- my_scale(list_xreg$predictors)
        X <- x_scaled$res

        if (method == "svd")
        {
          # inspired from MASS::lm.ridge
          Xs <- La.svd(X)
          rhs <- crossprod(Xs$u, y)
          d <- Xs$d
          nb_di <- length(d)
          div <- d^2 + rep(lambda, rep(nb_di, nlambda))
          a <- drop(d * rhs)/div
          dim(a) <- c(nb_di, nlambda)
          n <- nrow(X)

          if (nlambda == 1)
          {
            vt <- Xs$vt
            coef <- crossprod(vt, a)
            y_hat <- X %*% coef
            GCV <- colSums((y - y_hat)^2)/(n - colSums(matrix(d^2/div,
                                                                                nb_di)))^2

            if (compute_Sigma == TRUE)
            {
              rhsX <- crossprod(Xs$u, X)
              aX <- drop(d * rhsX)/div
              Sigma <- diag(ncol(X)) - crossprod(vt, aX)
              rownames(Sigma) <- colnames(X)
              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          Sigma = Sigma, lambda = lambda, ym = ym,
                          xm = x_scaled$xm, nb_hidden = nb_hidden,
                          nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            } else { #else: compute_Sigma == FALSE && nlambda == 1
              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          lambda = lambda, ym = ym, xm = x_scaled$xm,
                          nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm,
                          nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            }
          } else { #else: nlambda > 1
            coef <- crossprod(Xs$vt, a)
            colnames(coef) <- lambda
            y_hat <- X %*% coef
            fitted_values <- drop(ym +  y_hat)
            colnames(fitted_values) <- lambda
            GCV <- colSums((y - y_hat)^2)/(n - colSums(matrix(d^2/div,
                                                                                nb_di)))^2

            if (compute_Sigma == TRUE) #compute_Sigma == TRUE && nlambda > 1
            {
              rhsX <- crossprod(Xs$u, X)
              `%op%` <-  foreach::`%do%`
              i <- NULL
              Sigma <- foreach::foreach(i = 1:nlambda)%op%{
                div_i <- d^2 + rep(lambda[i], rep(nb_di, 1))
                aX_i <- drop(d * rhsX)/div_i
                Sigma <- diag(ncol(X)) - crossprod(Xs$vt, aX_i)
                rownames(Sigma) <- colnames(X)
                Sigma
              }
              names(Sigma) <- lambda
              return(list(coef = drop(coef), scales = x_scaled$xsd, Sigma = Sigma,
                          lambda = lambda, ym = ym, xm = x_scaled$xm, nb_hidden = nb_hidden,
                          nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = fitted_values,
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            } else { #else: compute_Sigma == FALSE && length(lambda) == 1

              return(list(coef = drop(coef), scales = x_scaled$xsd,
                          lambda = lambda, ym = ym, xm = x_scaled$xm,
                          nb_hidden = nb_hidden, nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
                          nodes_sim = nodes_sim, activ = activ,
                          fitted_values = drop(ym +  y_hat),
                          GCV = GCV,
                          compute_Sigma = compute_Sigma))
            }
          }
        }

        if (method == "chol")
        {
          XTX <- crossprod(X)

          if (nlambda > 1)
          {
            Dn <- vector("list", length = nlambda)
            names(Dn) <- lambda
            Sigma <- vector("list", length = nlambda)
            names(Sigma) <- lambda
            coef <- foreach::foreach(i = 1:nlambda, .combine = cbind)%do%{
              Dn[[i]] <- chol2inv(chol(XTX + diag(x = lambda[i],
                                                  nrow = nrow(XTX)))) # Cn^{-1}
              Sigma[[i]] <- diag(ncol(X)) - Dn[[i]]%*%XTX # Sigma_n
              Dn[[i]]%*%crossprod(X, y) # beta_n
            }
            colnames(coef) <- lambda
            rownames(coef) <- colnames(x_scaled$res)
          } else {
            Dn <- chol2inv(chol(XTX + diag(x = lambda,
                                           nrow = nrow(XTX)))) # Cn^{-1}
            Sigma <- diag(ncol(X)) - Dn%*%XTX # Sigma_n
            coef <- Dn%*%crossprod(X, y) # beta_n
          }

          return(list(coef = coef, Dn = Dn, Sigma = Sigma,
                      scales = x_scaled$xsd, lambda = lambda,
                      ym = ym, xm = x_scaled$xm,
                      nb_hidden = nb_hidden,
                      nodes_sim = nodes_sim,
                      activ = activ,
                      nn_xm = list_xreg$nn_xm,
                      nn_scales = list_xreg$nn_scales,
                      fitted_values = drop(ym + X %*% coef),
                      compute_Sigma = compute_Sigma,
                      x = x, y = y))
        }

    }
}

# Fitting Matérn 5/2 model
fit_matern52 <- function(x, y,
                   sigma = 2, l = 0.1, lambda_krls = 0.1,
                   inv_method = c("chol", "ginv"),
                   compute_Sigma = FALSE)
{
  if (!is.vector(y)) stop("'y' must be a vector") # otherwise y - ym is not working
  x <- as.matrix(x)
  y <- as.vector(y)
  inv_method <- match.arg(inv_method)

  ym <- mean(y)
  centered_y <- y - ym

  x_scaled <- my_scale(x)
  X <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  K <- bayesianrvfl::matern52_kxx_cpp(x = X,
                        sigma = sigma,
                        l = l)
  mat_coefs <- switch(
      inv_method,
      "chol" = chol2inv(chol(K + lambda_krls * diag(nrow(K)))) %*% centered_y,
      "ginv" = my_ginv(K + lambda_krls * diag(nrow(K))) %*% centered_y
    )

    lsfit <- drop(crossprod(K, mat_coefs))
    fitted_values <- lsfit  + rep(ym, length(lsfit))
    resid <- y - fitted_values

  return(
    list(
      y = y,
      sigma = sigma,
      l = l,
      K = K,
      lambda_krls = lambda_krls,
      inv_method = inv_method,
      mat_coefs = mat_coefs,
      fitted_values = fitted_values,
      ym = ym,
      xm = xm,
      xsd = xsd,
      scaled_x = X,
      resid = resid,
      compute_Sigma = compute_Sigma
    )
  )

}

# Fitting elastic net
fit_glmnet <- function(x, y, nb_hidden = 5,
                       nodes_sim = c("sobol", "halton", "unif"),
                       activ = c("relu", "sigmoid", "tanh",
                                 "leakyrelu", "elu", "linear"),
                       ...)
{
  if (!is.vector(y)) stop("'y' must be a vector") # otherwise y - ym is not working
  x <- as.matrix(x)
  y <- as.vector(y)

  ym <- mean(y)
  centered_y <- y - ym

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  list_xreg <- bayesianrvfl::create_new_predictors(x = x, nb_hidden = nb_hidden,
                                     nodes_sim = nodes_sim,
                                     activ = activ)

  x_scaled <- my_scale(list_xreg$predictors)
  X <- x_scaled$res

  fit_obj <- glmnet::glmnet(x = X, y = y, ...)

  # obtain fitted values

  fitted_values <- predict.elnet(fit_obj, newx = X)

  return(
    list(fit_obj = fit_obj,
         scales = x_scaled$xsd,
         ym = ym, xm = x_scaled$xm,
         nb_hidden = nb_hidden,
         nn_xm = list_xreg$nn_xm, nn_scales = list_xreg$nn_scales,
         nodes_sim = nodes_sim, activ = activ,
         y = y, fitted_values = fitted_values,
         resid = y - fitted_values
    )
  )
}


sigmoid <- function(x) exp(x)/(1 + exp(x))
n_iris <- nrow(iris)
index_train <- sort(sample(1:n_iris, size = floor(0.8*n_iris)))
x_train <- as.matrix(iris[index_train, 1:4])
X_train <- cbind(1, x_train)
y_train <- factor_to_matrix(iris$Species[index_train])
fit_lm <- .lm.fit(x = X_train, y = y_train)

x_test <- as.matrix(iris[-index_train, 1:4])
X_test <- cbind(1, x_test)
yhat <- X_test%*%fit_lm$coefficients
preds <- sapply(1:3, function (i) as.numeric(yhat[,i] > 0.5))
colnames(preds) <- colnames(y_test)
rownames(preds) <- 1:nrow(preds)
probs <- sigmoid(yhat)
y_test <- factor_to_matrix(iris$Species[-index_train])


colMeans(y_test == preds)


