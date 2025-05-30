% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_functions.R
\name{msnlminb}
\alias{msnlminb}
\title{Multistart nlminb}
\usage{
msnlminb(objective, nb_iter = 100, lower, upper, cl = NULL, max_fails = 3, ...)
}
\arguments{
\item{objective}{function to be minimized}

\item{nb_iter}{number of iterations}

\item{lower}{lower bounds}

\item{upper}{upper bounds}

\item{cl}{number of cores to be used}

\item{max_fails}{maximum consecutive failures before adjusting parameters}

\item{...}{other arguments to be passed to nlminb}
}
\value{
list with the best solution
}
\description{
Multistart nlminb
}
