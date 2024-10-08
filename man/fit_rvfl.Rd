% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit_rvfl}
\alias{fit_rvfl}
\title{Fitting base rvfl}
\usage{
fit_rvfl(
  x,
  y,
  nb_hidden = 5,
  n_clusters = 0,
  nodes_sim = c("sobol", "halton", "unif"),
  activ = c("relu", "sigmoid", "tanh", "leakyrelu", "elu", "linear"),
  lambda = 10^seq(from = -10, to = 10, length.out = 100),
  method = c("svd", "solve", "chol"),
  compute_Sigma = FALSE,
  seed = 123
)
}
\arguments{
\item{seed}{}
}
\description{
Fitting base rvfl
}
