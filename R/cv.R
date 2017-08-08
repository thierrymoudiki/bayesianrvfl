# inspired from caret::createFolds
create_folds <- function(y, k = 10)
{
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)

  out <- split(seq(along = y), foldVector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                      sep = "")


  return(out)
}
create_folds <- compiler::cmpfun(create_folds)

compute_accuracy <- function(x, y, nb_hidden = 5,
                             nodes_sim = c("sobol", "halton", "unif"),
                             activ = c("relu", "sigmoid", "tanh",
                                       "leakyrelu", "elu", "linear"),
                             lambda = 10^seq(-2, 10, length.out = 100),
                             k = 5, repeats = 1, seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  if(round(nrow(x)/k) < 3)
    warnings("Risk of having empty folds, pick a higher k")

  `%op%` <-  foreach::`%do%`

  if(is.wholenumber(repeats) && repeats > 1)
  {
    set.seed(seed)
    list_folds <- lapply(1:repeats,
                         function (i) create_folds(y = y, k = k))
    i <- j <- NULL
    res <- foreach::foreach(j = 1:repeats, .combine = 'rbind', .errorhandling = "stop")%op%{
      temp <- foreach::foreach(i = 1:k, .combine = 'rbind', .errorhandling = "stop")%op%{
        train_index <- -list_folds[[j]][[i]]
        test_index <- -train_index
        fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
                            nodes_sim = nodes_sim, activ = activ,
                            nb_hidden = nb_hidden, lambda = lambda,
                            compute_Sigma = FALSE)
        predict_rvfl(fit_obj, newx = x[test_index, ]) - y[test_index]
      }
    }

    #return(res)
    return(sqrt(colMeans(res^2)))
  } else {

    set.seed(seed)
    folds <- create_folds(y = y, k = k)

    res <- foreach::foreach(i = 1:k, .combine = rbind)%op%{
      train_index <- -folds[[i]]
      test_index <- -train_index
      fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
                          nodes_sim = nodes_sim, activ = activ,
                          nb_hidden = nb_hidden, lambda = lambda,
                          compute_Sigma = FALSE)
      predict_rvfl(fit_obj, newx = x[test_index, ]) - y[test_index]
    }

    #return(res)
    return(sqrt(colMeans(res^2)))
  }
}
compute_accuracy <- compiler::cmpfun(compute_accuracy)

cv_rvfl <- function(x, y, k = 5, repeats = 10,
                    nodes_sim = c("sobol", "halton", "unif"),
                    activ = c("relu", "sigmoid", "tanh",
                              "leakyrelu", "elu", "linear"),
                    vec_nb_hidden = seq(from = 100, to = 1000, by = 100),
                    lams = 10^seq(-2, 10, length.out = 100), seed = 1,
                    cl = NULL)
{

  nb_iter <- length(vec_nb_hidden)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  allowParallel <- !is.null(cl) && cl > 0
  if(allowParallel)
  {
    cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
    doSNOW::registerDoSNOW(cl_SOCK)
    `%op%` <-  foreach::`%dopar%`
  }  else {
    `%op%` <-  foreach::`%do%`
  }

  pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  i <- NULL
  res <- foreach::foreach(i = 1:nb_iter, .packages = "doSNOW",
                 .combine = rbind,
                 .options.snow = opts, .verbose = FALSE,
                 .export = c("compute_accuracy",
                             "is.wholenumber",
                             "create_folds",
                             "fit_rvfl",
                             "predict_rvfl",
                             "create_new_predictors",
                             "my_scale",
                             "remove_zero_cols",
                             "my_sd"))%op%
  {
   as.vector(compute_accuracy(x = x, y = y,
                              nb_hidden = vec_nb_hidden[i],
                              nodes_sim = nodes_sim, activ = activ,
                              k = k,
                              repeats = repeats,
                              lambda = lams, seed = seed))
  }
  close(pb)
  if(allowParallel) snow::stopCluster(cl_SOCK)

  colnames(res) <- lams
  rownames(res) <- vec_nb_hidden
  return(res)
}

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
#
# # 2 - Exemple SLC14 ---------------------------------------------------------
#
#  library(caret)
#  set.seed(7210)
#  train_dat <- SLC14_1(250)
#  #large_dat <- SLC14_1(10000)
#
#  head(train_dat)
#  str(train_dat)
#
#  x <- train_dat[, -ncol(train_dat)]
#  y <- train_dat[, ncol(train_dat)]
# train_index <- createDataPartition(y, p = 0.7)[[1]]
#
# #lams <- 10^seq(0, 2, length.out = 100)
# lams <- seq(10, 20, length.out = 100)
# vec_nb_hidden <- seq(700, 900, by = 25)
# res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
#                   k = 10, repeats = 5,
#                   vec_nb_hidden = vec_nb_hidden,
#                   lams = lams,
#                   seed = 1, cl = 4)
#
#   coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE)
#   res_cv[coord_min[1], coord_min[2]]
#   (best_nb_hidden <- vec_nb_hidden[coord_min[1]])
#   (best_lam <- lams[coord_min[2]])
#   summary(y)
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
#
#   # 3 - Exemple mtcars ---------------------------------------------------------
#
  # x <- as.matrix(mtcars[, -1])
  # y <- as.vector(mtcars[, 1])
  #
  # train_index <- createDataPartition(y, p = 0.7)[[1]]
  #
  # lams <- seq(10, 20, length.out = 100)
  # vec_nb_hidden <- 1:10

  # best_lam: 12.22222, cv: 2.58
  #nodes_sim <- "sobol"
  #activ <- "elu"
  # best_lam: 14.34343, cv: 2.505841
  #nodes_sim <- "halton"
  #activ <- "elu"
  # best_lam: , cv: 2.505841
  # best_lam: 8.030303, cv: 2.78051
  #nodes_sim <- "unif"
  #activ <- "elu"

  # best_lam: 12.22222 cv: 2.587276
  #nodes_sim <- "sobol"
  #activ <- "relu"
  # best_lam: 14.44444 cv: 2.505732
  #nodes_sim <- "halton"
  #activ <- "relu"
  # best_lam: 8 cv: 2.780729
  #nodes_sim <- "unif"
  #activ <- "relu"

  # best_lam: 12.22222 cv: 2.589279
  #nodes_sim <- "sobol"
  #activ <- "leakyrelu"
  # best_lam: 14.34343 cv: 2.504902
  #nodes_sim <- "halton"
  #activ <- "leakyrelu"
  # best_lam: 8 cv: 2.781534
  #nodes_sim <- "unif"
  #activ <- "leakyrelu"

  # best_lam: 20.25502 cv: 2.573572
  #nodes_sim <- "sobol"
  #activ <- "tanh"
  # best_lam: 11.62322 cv: 2.634311
  #nodes_sim <- "halton"
  #activ <- "tanh"
  # best_lam: 4.824109  cv: 2.682918
  #nodes_sim <- "unif"
  #activ <- "tanh"

  # best_lam: 79.34097 cv: 2.555117
  #nodes_sim <- "sobol"
  #activ <- "sigmoid"
  # best_lam: 17.22586 cv: 2.785476
  #nodes_sim <- "halton"
  #activ <- "sigmoid"
  # best_lam: 3.255089 cv: 2.774128
  #nodes_sim <- "unif"
  #activ <- "sigmoid"

  # nodes_sim <- "halton"
  # activ <- "relu"
  # lams <- 10^seq(0, 2, length.out = 200)
  # vec_nb_hidden <- 1:10
  #
  # res_cv <- cv_rvfl(x = x[train_index, ], y = y[train_index],
  #                   nodes_sim = nodes_sim, activ = activ,
  #                   k = 5, repeats = 10,
  #                   vec_nb_hidden = vec_nb_hidden,
  #                   lams = lams,
  #                   seed = 1, cl = 4)
  #
  # (coord_min <- which(res_cv == min(res_cv), arr.ind = TRUE))
  # res_cv[coord_min[1], coord_min[2]]
  # (best_nb_hidden <- vec_nb_hidden[coord_min[1]])
  # (best_lam <- lams[coord_min[2]])
  # summary(y)
  #
  # fit_obj <- fit_rvfl(x = x[train_index, ], y = y[train_index],
  #                     nodes_sim = nodes_sim, activ = activ,
  #                     nb_hidden = best_nb_hidden, lambda = best_lam,
  #                     compute_Sigma = TRUE)
  # (preds <- predict_rvfl(fit_obj, newx = x[-train_index, ]))
  #
  # sqrt(mean((preds$mean - y[-train_index])^2))
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
