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
  x <- as.matrix(x)
  y <- as.vector(y)

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

