# 0 - data -----

library(MASS)
data("Boston")
X <- as.matrix(Boston[, -ncol(Boston)])
y <- Boston[, ncol(Boston)]
n <- dim(X)[1]; p <- dim(X)[2]

#set.seed(123)
#X <- matrix(rnorm(n * p), n, p) #no intercept!
#y <- rnorm(n)

(idx_train <- caret::createDataPartition(y, p = 0.6)$Resample1)
(idx_comp <- generics::setdiff(1:n, idx_train))
(idx_validation <- base::sample(x = idx_comp,
                               size = floor(length(idx_comp)/2),
                               replace = FALSE))
(idx_test <- generics::setdiff(idx_comp, idx_validation))
generics::intersect(idx_train, idx_validation)
generics::intersect(c(idx_train, idx_validation), idx_test)

X_train <- X[idx_train,]
y_train <- y[idx_train]
X_validation <- X[idx_validation,]
y_validation <- y[idx_validation]
X_test <- X[idx_test,]
y_test <- y[idx_test]

obj_GCV <- bayesianrvfl::fit_rvfl(x = X_train, y = y_train)
(best_lambda <- obj_GCV$lambda[which.min(obj_GCV$GCV)])

(fit_obj <- bayesianrvfl::fit_rvfl(x = X_train,
                                   y = y_train,
                                   method = "solve",
                                   lambda = best_lambda,
                                   compute_Sigma = TRUE))

preds_validation <- bayesianrvfl::predict_rvfl(fit_obj,
                                               newx = X_validation)
level <- 99
summary(preds_validation$mean - y_validation)
preds_upper <- preds_validation$mean + qnorm(1 - (100 - level)/200)*preds_validation$sd
preds_lower <- preds_validation$mean - qnorm(1 - (100 - level)/200)*preds_validation$sd

# coverage
mean((preds_upper >= y_validation)*(preds_lower <= y_validation))

par(mfrow=c(3, 2))
plot(x = log(obj_GCV$lambda), y = obj_GCV$GCV, type='l')
plot(y_validation, type='l', col="red")
lines(preds_validation$mean, col="blue")
lines(preds_upper, col="gray")
lines(preds_lower, col="gray")
plot(x = y_validation,
     y = preds_validation$mean,
     ylim = c(min(c(y_validation, preds_validation$mean)),
              max(c(y_validation, preds_validation$mean))))
abline(a = 0, b = 1, col="green", lwd=2)

fit_obj2 <- bayesianrvfl::update_params(fit_obj, newx = X[idx_test[1], ],
                                        newy = y[idx_test[1]])
fit_obj3 <- bayesianrvfl::update_params(fit_obj2, newx = X[idx_test[2], ],
                                        newy = y[idx_test[2]])
fit_obj4 <- bayesianrvfl::update_params(fit_obj3, newx = X[idx_test[3], ],
                                        newy = y[idx_test[3]])
fit_obj5 <- bayesianrvfl::update_params(fit_obj4, newx = X[idx_test[4], ],
                                        newy = y[idx_test[4]])
fit_obj6 <- bayesianrvfl::update_params(fit_obj5, newx = X[idx_test[5], ],
                                        newy = y[idx_test[5]])

mat_coefs <- cbind(fit_obj$coef, fit_obj2$coef,
                    fit_obj3$coef, fit_obj4$coef,
                    fit_obj5$coef, fit_obj6$coef)

matplot(t(mat_coefs), type = 'l')


preds_validation2 <- bayesianrvfl::predict_rvfl(fit_obj5,
                                               newx = X_validation)
plot(x = y_validation,
     y = preds_validation2$mean,
     ylim = c(min(c(y_validation, preds_validation2$mean)),
              max(c(y_validation, preds_validation2$mean))))
abline(a = 0, b = 1, col="green", lwd=2)

preds_validation3 <- bayesianrvfl::predict_rvfl(fit_obj6,
                                                newx = X_validation)
plot(x = y_validation,
     y = preds_validation3$mean,
     ylim = c(min(c(y_validation, preds_validation3$mean)),
              max(c(y_validation, preds_validation3$mean))))
abline(a = 0, b = 1, col="green", lwd=2)
