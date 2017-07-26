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
                 lambda = 10^seq(-2, 10, length.out = 100),
                 k = 5, repeats = 1, seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))

  if(round(nrow(x)/k) < 3)
    stop("Risk of having empty folds, pick a higher k")

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
                            nb_hidden = nb_hidden, lambda = lambda)
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
                          nb_hidden = nb_hidden, lambda = lambda)
      predict_rvfl(fit_obj, newx = x[test_index, ]) - y[test_index]
    }

    #return(res)
    return(sqrt(colMeans(res^2)))
  }
}
compute_accuracy <- compiler::cmpfun(compute_accuracy)

cv_rvfl <- function(x, y, k = 5, repeats = 10,
                    vec_nb_hidden = seq(from = 100, to = 1000, by = 100),
                    lams = 10^seq(-2, 10, length.out = 100), seed = 1,
                    cl = NULL)
{

  nb_iter <- length(vec_nb_hidden)

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
                              nb_hidden = vec_nb_hidden[i], k = k,
                              repeats = repeats,
                              lambda = lams, seed = seed))
  }
  close(pb)
  if(allowParallel) snow::stopCluster(cl_SOCK)

  colnames(res) <- lams
  rownames(res) <- vec_nb_hidden
  return(res)
}

# 1 - Exemple longley ---------------------------------------------------------

# lams <- 10^seq(-2, 10, length.out = 1000)
# x <- longley[,-1]
# y <- longley[,1]
# (res <- compute_accuracy(x = x, y = y, nb_hidden = 100, k = 5,
#                          repeats = 10, lambda = lams, seed = 2))
# plot(log(lams), res, type = 'l')
# res[which.min(res)]
# lams[which.min(res)]

# 2 - Exemple SLC14 ---------------------------------------------------------

 library(caret)
 set.seed(7210)
 train_dat <- SLC14_1(250)
 large_dat <- SLC14_1(10000)
#
 head(train_dat)
 str(train_dat)
#
 x <- train_dat[, -ncol(train_dat)]
 y <- train_dat[, ncol(train_dat)]
# ncol(x)
# lams <- 10^seq(-5, 10, length.out = 500)
# (res <- compute_accuracy(x = x, y = y, nb_hidden = 500, k = 10,
#                          repeats = 5, lambda = lams, seed = 2))
# plot(log(lams), res, type = 'l')
# res[which.min(res)]
# lams[which.min(res)]
#
# # parallel
#  vec_nb_hidden <- seq(from = 100, to = 2000, by = 100)
#  lams <- seq(0, 0.0001, length.out = 100)
# #
#  ptm <- proc.time()
#  cv_res <- cv_rvfl(x = x, y = y, vec_nb_hidden = vec_nb_hidden,
#                    lams = lams, k = 10, repeats = 5, seed = 308, cl = 4)
#  proc.time() - ptm
# #
#  coord_min <- which(cv_res == min(cv_res), arr.ind = TRUE)
#  cv_res[coord_min[1], coord_min[2]]
#  (best_nb_hidden <- vec_nb_hidden[coord_min[1]])
#  (best_lam <- lams[coord_min[2]])
# #
#  filled.contour(x = vec_nb_hidden,
#                 y = lams, z = cv_res,
#                 color.palette = terrain.colors,
#                 xlab = "nb hidden", ylab = "lambda")
#
#  filled.contour(x = vec_nb_hidden,
#                 y = log(lams), z = cv_res,
#                 color.palette = terrain.colors,
#                 xlab = "nb hidden", ylab = "log(lambda)")
#
#
folds <- createDataPartition(y, p = 0.9)
fit_obj <- fit_rvfl(x = as.matrix(x[folds$Resample1, ]), y = y[folds$Resample1], nb_hidden = 100, lambda = 1)
preds <- predict_rvfl(fit_obj, newx = as.matrix(x[-folds$Resample1, ]))

upper <- preds$mean + 1.96*sqrt(diag(preds$cov))
lower <- preds$mean - 1.96*sqrt(diag(preds$cov))

plot(y[-folds$Resample1], type = 'l')


plot(preds$mean, type = "l", col = "blue")
lines(upper, type = "l", col = "red")
lines(lower, type = "l", col = "red")
