# 1 - multistart nlminb ---------------------------------------------------------------

#' @title Multistart nlminb
#' @description Multistart nlminb
#' @param objective function to be minimized
#' @param nb_iter number of iterations
#' @param lower lower bounds
#' @param upper upper bounds
#' @param cl number of cores to be used
#' @param max_fails maximum consecutive failures before adjusting parameters
#' @param ... other arguments to be passed to nlminb
#' @return list with the best solution
#' @export
#' @examples
msnlminb <- function(objective, nb_iter = 100,
                     lower, upper, cl = NULL,
                     max_fails = 3, # Maximum consecutive failures before adjusting parameters
                     ...)
{
  ldots <- list(...)
  
  # Robust objective function wrapper
  OF <- function(u, ...) {
    tryCatch({
      val <- objective(u, ...)
      # Check if result is valid
      if (is.null(val) || !is.finite(val) || is.na(val)) {
        return(1e10)  # Return high value for invalid results
      }
      return(val)
    }, error = function(e) {
      return(1e10)  # Return high value on error
    })
  }

  rep_1_nb_iter <- rep(1, nb_iter)
  lower_mat <- tcrossprod(rep_1_nb_iter, lower)
  upper_mat <- tcrossprod(rep_1_nb_iter, upper)

  # Generate starting points using Sobol sequences
  starting_points <- lower_mat + (upper_mat - lower_mat)*randtoolbox::sobol(n = nb_iter,
                                                                           dim = length(lower))
  nb_iter <- nrow(starting_points)
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)
  
  if (is.null(cl)) {
    res <- vector("list", nb_iter)
    consecutive_fails <- 0
    
    for(i in 1:nb_iter) {
      setTxtProgressBar(pb, i)
      
      current_result <- tryCatch({
        opt <- stats::nlminb(
          start = starting_points[i, ],
          objective = OF,
          lower = lower,
          upper = upper,
          control = list(rel.tol = 1e-6, x.tol = 1e-6),
          ...
        )
        
        if (!is.null(opt) && is.finite(opt$objective)) {
          consecutive_fails <- 0  # Reset counter on success
          opt
        } else {
          consecutive_fails <- consecutive_fails + 1
          list(objective = 1e10)
        }
      }, error = function(e) {
        consecutive_fails <- consecutive_fails + 1
        list(objective = 1e10)
      })
      
      res[[i]] <- current_result
      
      # If too many consecutive failures, try different starting point
      if (consecutive_fails >= max_fails) {
        starting_points[i, ] <- lower + (upper - lower) * stats::runif(length(lower))
        consecutive_fails <- 0
      }
    }
    close(pb)
  } else {
    # ... rest of parallel implementation ...
  }

  # Filter out failed optimizations with more lenient criteria
  valid_results <- !sapply(res, is.null) & 
                  sapply(res, function(x) x$objective < 1e10)
  
  if(sum(valid_results) == 0) {
    # If all attempts failed, return best of the worst
    warning("All optimization attempts had high objective values")
    index_opt <- which.min(sapply(res, function(x) x$objective))
    return(res[[index_opt]])
  }
  
  # Filter valid results and return best
  res <- res[valid_results]
  index_opt <- which.min(sapply(res, function(x) x$objective))
  
  return(res[[index_opt]])
}

# 2 - random search ---------------------------------------------------------------

random_search_opt <- function(objective, nb_iter = 100,
                              lower, upper, sim = c("sobol", "unif"),
                              seed = 123, cl = NULL, ...)
{
  OF <- function(y, ...) {return(objective(y, ...))}
  rep_1_nb_iter <- rep(1, nb_iter)
  lower_mat <- tcrossprod(rep_1_nb_iter, lower)
  upper_mat <- tcrossprod(rep_1_nb_iter, upper)
  sim <- match.arg(sim)

  if (sim == "sobol")
  {
    searched_points <- lower_mat + (upper_mat -
                                      lower_mat) * randtoolbox::sobol(n = nb_iter,
                                                                      dim = length(lower))
  }

  if (sim == "unif")
  {
    set.seed(seed)
    searched_points <- lower_mat + (upper_mat -
                                      lower_mat) * matrix(runif(nb_iter *
                                                                  length(lower)),
                                                          nrow = nrow(lower_mat))
  }

  if (is.null(cl)) {
    pb <- txtProgressBar(min = 1,
                         max = nb_iter,
                         style = 3)

    `%op%` <-  foreach::`%do%`

    res <- foreach::foreach(
      i = 1:nb_iter,
      .combine = c,
      .verbose = FALSE,
      .errorhandling = "remove"
    ) %op% {
      setTxtProgressBar(pb, i)
      res <- try(OF(searched_points[i, ], ...),
                 silent = TRUE)
      ifelse(class(res) == "try-error", 1e06, res)
    }

    close(pb)

  } else {
    cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
    doSNOW::registerDoSNOW(cl_SOCK)

    pb <- txtProgressBar(min = 0,
                         max = nb_iter,
                         style = 3)
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    i <- j <- NULL
    res <-
      suppressWarnings(
        foreach::foreach(
          i = 1:nb_iter,
          .combine = c,
          .packages = "doSNOW",
          .options.snow = opts,
          .export = ...,
          .verbose = FALSE,
          .errorhandling = "remove"
        ) %dopar% {
          res <- try(OF(searched_points[i, ], ...),
                     silent = TRUE)
          ifelse(class(res) == "try-error", 1e06, res)
        }
      )
    snow::stopCluster(cl_SOCK)
  }

  index_opt <- which.min(res)

  return(list(par = searched_points[index_opt, ],
              objective = res[index_opt]))
}
random_search_opt <- memoise::memoise(f = random_search_opt)



# 3 - multistart nmkb ---------------------------------------------------------------

#' @title Multistart nmkb
#' @description Multistart nmkb
#' @param objective function to be minimized
#' @param nb_iter number of iterations
#' @param lower lower bounds
#' @param upper upper bounds
#' @param cl number of cores to be used
#' @param ... other arguments to be passed to nmkb
#' @return list with the best solution
#' @export
#' @examples
msnmkb <- function(objective, nb_iter = 100,
                     lower, upper, cl = NULL,
                     ...)
{
  OF <- function(u, ...)
  {return(objective(u, ...))}

  rep_1_nb_iter <- rep(1, nb_iter)
  lower_mat <- tcrossprod(rep_1_nb_iter, lower)
  upper_mat <- tcrossprod(rep_1_nb_iter, upper)

  starting_points <- lower_mat + (upper_mat - lower_mat)*randtoolbox::sobol(n = nb_iter,
                                                                            dim = length(lower))
  nb_iter <- nrow(starting_points)

  pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)

  if (is.null(cl)) {

    res <- foreach::foreach(i = 1:nb_iter,
                            .export = "ldots",
                            .verbose = FALSE,
                            .packages = "dfoptim",
                            .errorhandling = "remove")%do%{
                              setTxtProgressBar(pb, i)
                              dfoptim::nmkb(par = starting_points[i, ],
                                            fn = OF,
                                            lower = lower,
                                            upper = upper,
                                            ...)
                            }
    close(pb)

  } else {

    cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
    doSNOW::registerDoSNOW(cl_SOCK)

    pb <- txtProgressBar(min = 0, max = nb_iter,
                         style = 3)

    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    i <- j <- NULL
    res <- foreach::foreach(i = 1:nb_iter,
                            .packages = c("doSNOW", "dfoptim"),
                            .options.snow = opts,
                            .verbose = FALSE,
                            .errorhandling = "remove")%dopar%{ # fix this

                              dfoptim::nmkb(par = starting_points[i, ],
                                            fn = OF,
                                            lower = lower,
                                            upper = upper,
                                            ...)

                            }

    parallel::stopCluster(cl_SOCK)
  }

  index_opt <- which.min(sapply(1:length(res),
                                function (i)
                                  res[[i]]$value))

  return(res[[index_opt]])
}
