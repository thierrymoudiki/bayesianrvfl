# # 1 - Exemple longley ---------------------------------------------------------
#
# lams <- 10^seq(1, 2, length.out = 100)
# vec_nb_hidden <- seq(600, 850, by = 10)
# x <- data.matrix(longley[, 1:6])
# y <- longley[, "Employed"]
# train_index <- createDataPartition(y, p = 0.7)[[1]]
#
# res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
#                   k = 3, repeats = 10,
#                   vec_nb_hidden = vec_nb_hidden,
#                   lams = lams,
#                   seed = 1, cl = 4)
#
#   coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE)
#   res_cv[coord_min[1], coord_min[2]]
#   (best_nb_hidden <- vec_nb_hidden[coord_min[1]])
#   (best_lam <- lams[coord_min[2]])
#
#   fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
#                       nb_hidden = best_nb_hidden, lambda = best_lam,
#                       compute_Sigma = TRUE)
#   (preds <- predict_rvfl(fit_obj, newx = x[-train_index, ]))
#
#   sqrt(mean((preds$mean - y[-train_index])^2))
#
#
# qt95 <- qnorm(1-0.05/2)
# qt80 <- qnorm(1-0.1/2)
# upper <- preds$mean + qt95*preds$sd
# lower <- preds$mean - qt95*preds$sd
# upper80 <- preds$mean + qt80*preds$sd
# lower80 <- preds$mean - qt80*preds$sd
#
# yy <- c(lower, rev(upper))
# yy80 <- c(lower80, rev(upper80))
# nbyy <- length(upper)
# xx <- c(1:nbyy, nbyy:1)
#
# plot(1:nbyy, y[-train_index], type = 'l',
#      ylim = c(min(lower), max(upper)), lwd = 2)
# polygon(xx, yy, col = "gray80", border = FALSE)
# polygon(xx, yy80, col = "gray60", border = FALSE)
# lines(1:nbyy, y[-train_index], lwd = 2)
# lines(preds$mean, col = "blue", lwd = 2)

# nodes_sim <- "halton"
# activ <- "tanh"
# lams <- 10^seq(-2, 2, length.out = 100)
# vec_nb_hidden <- 1:10
# res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
#                   nodes_sim = nodes_sim, activ = activ,
#                   k = 3, repeats = 10,
#                   vec_nb_hidden = vec_nb_hidden,
#                   lams = lams,
#                   seed = 1, cl = 4)
#
#   coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE)
#   res_cv[coord_min[1], coord_min[2]]
#   (best_nb_hidden <- vec_nb_hidden[coord_min[1]])
#   (best_lam <- lams[coord_min[2]])

# halton relu 3.673985
# sobol relu 0.7814588
# unif relu  2.500851
# unif elu  2.502318
# halton elu 3.673985
# sobol elu 0.7846292
# unif leakyrelu 2.505171
# sobol leakyrelu 0.7832614
# halton leakyrelu 3.673985
# halton sigmoid 1.148266
# sobol sigmoid 0.8781506
# unif sigmoid 0.9351309
# unif tanh 0.7586961
# sobol tanh 0.6969438

#
# # 2 - Exemple SLC14 ---------------------------------------------------------
#
library(caret)
set.seed(123)
train_dat <- SLC14_1(250)

head(train_dat)
str(train_dat)

x <- train_dat[, -ncol(train_dat)]
y <- train_dat[, ncol(train_dat)]
train_index <- createDataPartition(y, p = 0.7)[[1]]

lams <- seq(10, 20, length.out = 100)
vec_nb_hidden <- seq(700, 900, by = 25)
res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
                  k = 10, repeats = 5,
                  vec_nb_hidden = vec_nb_hidden,
                  lams = lams,
                  seed = 1, cl = 4)

coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE)
res_cv[coord_min[1], coord_min[2]]
(best_nb_hidden <- vec_nb_hidden[coord_min[1]])
(best_lam <- lams[coord_min[2]])
summary(y)

fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
                    nb_hidden = best_nb_hidden, lambda = best_lam,
                    compute_Sigma = TRUE)
(preds <- predict_rvfl(fit_obj, newx = x[-train_index, ]))

sqrt(mean((preds$mean - y[-train_index])^2))


qt95 <- qnorm(1-0.05/2)
qt80 <- qnorm(1-0.1/2)
upper <- preds$mean + qt95*preds$sd
lower <- preds$mean - qt95*preds$sd
upper80 <- preds$mean + qt80*preds$sd
lower80 <- preds$mean - qt80*preds$sd

yy <- c(lower, rev(upper))
yy80 <- c(lower80, rev(upper80))
nbyy <- length(upper)
xx <- c(1:nbyy, nbyy:1)

plot(1:nbyy, y[-train_index], type = 'l',
     ylim = c(min(lower), max(upper)), lwd = 2)
polygon(xx, yy, col = "gray80", border = FALSE)
polygon(xx, yy80, col = "gray60", border = FALSE)
lines(1:nbyy, y[-train_index], lwd = 2)
lines(preds$mean, col = "blue", lwd = 2)



lams <- 10^seq(-10, 2, length.out = 100)
vec_nb_hidden <- seq(2000, 2500, by = 10)
res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
                  nodes_sim = "sobol", activ = "leakyrelu",
                  k = 10, repeats = 5,
                  vec_nb_hidden = vec_nb_hidden,
                  lams = lams,
                  seed = 1, cl = 4)
coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE)
res_cv[coord_min[1], coord_min[2]]
(best_nb_hidden <- vec_nb_hidden[coord_min[1]])
(best_lam <- lams[coord_min[2]])


#
#   # 3 - Exemple mtcars ---------------------------------------------------------
#
x <- as.matrix(mtcars[, -1])
y <- as.vector(mtcars[, 1])

train_index <- createDataPartition(y, p = 0.7)[[1]]

lams <- seq(10, 20, length.out = 100)
vec_nb_hidden <- 1:10

nodes_sim <- "halton"
activ <- "relu"
lams <- 10^seq(0, 2, length.out = 200)
vec_nb_hidden <- 1:10

res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
                  nodes_sim = nodes_sim, activ = activ,
                  k = 5, repeats = 10,
                  vec_nb_hidden = vec_nb_hidden,
                  lams = lams,
                  seed = 1, cl = 4)

(coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE))
res_cv[coord_min[1], coord_min[2]]
(best_nb_hidden <- vec_nb_hidden[coord_min[1]])
(best_lam <- lams[coord_min[2]])

fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
                    nodes_sim = nodes_sim, activ = activ,
                    nb_hidden = best_nb_hidden, lambda = best_lam,
                    compute_Sigma = TRUE)
(preds <- predict_rvfl(fit_obj, newx = x[-train_index, ]))

sqrt(mean((preds$mean - y[-train_index])^2))

qt95 <- qnorm(1-0.05/2)
qt80 <- qnorm(1-0.1/2)
upper <- preds$mean + qt95*preds$sd
lower <- preds$mean - qt95*preds$sd
upper80 <- preds$mean + qt80*preds$sd
lower80 <- preds$mean - qt80*preds$sd

yy <- c(lower, rev(upper))
yy80 <- c(lower80, rev(upper80))
nbyy <- length(upper)
xx <- c(1:nbyy, nbyy:1)

plot(1:nbyy, y[-train_index], type = 'l',
     ylim = c(min(lower), max(upper)), lwd = 2)
polygon(xx, yy, col = "gray80", border = FALSE)
polygon(xx, yy80, col = "gray60", border = FALSE)
lines(1:nbyy, y[-train_index], lwd = 2)
lines(preds$mean, col = "blue", lwd = 2)
