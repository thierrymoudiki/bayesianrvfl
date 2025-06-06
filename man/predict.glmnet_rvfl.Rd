% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.glmnet_rvfl}
\alias{predict.glmnet_rvfl}
\title{Predict Method for Elastic Net Models}
\usage{
\method{predict}{glmnet_rvfl}(object, newdata, s = 0.1, ...)
}
\arguments{
\item{object}{Fitted model object}

\item{newdata}{New data for prediction}

\item{s}{Value of the penalty parameter lambda at which predictions are required}

\item{...}{Additional arguments}
}
\value{
Predictions for new data
}
\description{
Predict Method for Elastic Net Models
}
