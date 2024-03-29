% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_estimParam.R
\name{saemixPredictNewdata}
\alias{saemixPredictNewdata}
\alias{estimateMeanParametersNewdata}
\alias{estimateIndividualParametersNewdata}
\title{Predictions for a new dataset}
\usage{
saemixPredictNewdata(
  saemixObject,
  newdata,
  type = c("ipred", "ypred", "ppred", "icpred"),
  nsamp = 1
)
}
\arguments{
\item{saemixObject}{an SaemixObject from a fitted run}

\item{newdata}{a dataframe containing the new data. The dataframe must contain the same information as the original dataset (colunm names, etc...)}

\item{type}{one or several of "ipred" (individual predictions using the MAP estimates), "ppred" (population predictions obtained using the population parameters f(E(theta))), "ypred" (mean of the population predictions (E(f(theta)))), "icpred"  (individual predictions using the conditional mean estimates). Defaults to "ppred".}

\item{nsamp}{an integer, ignored for other types than icpred; if icpred, returns both the mean of the conditional distribution and nsamp samples, with the corresponding predictions. Defaults to 1.}
}
\value{
a list with three components (five if type includes "icpred" and nsamp>1)
\describe{
\item{param}{a dataframe with the estimated parameters. The columns in the dataframe depend on which type of predictions were requested (argument type)}
\item{predictions}{a dataframe with the predictions. The columns in the dataframe depend on which type of predictions were requested (argument type)}
\item{saemixObject}{the SaemixObject with the data slot replaced by the new data. The elements of the results slot pertaining to individual (including population individual parameters) predictions and likelihood will have been removed, and only the elements computed within the function will have been replaced (eg individual estimated parameters and predictions for the new data)}
\item{parSample}{a dataframe with parameters sampled from the conditional distribution of the individual parameters (only present if type includes 'icpred' and nsamp>1)}
\item{predSample}{a dataframe with the predictions corresponding to the parameters sampled from the conditional distribution of the individual parameters (only present if type includes 'icpred' and nsamp>1)}
}
}
\description{
Predictions for a new dataset
}
\details{
This function is the workhorse behind the predict method for SaemixObject. It computes predictions for a new dataframe
based on the results of an saemix fit.
Since the predict function only returns predicted values, this function is provided so that users can access other elements,
for example the different types of parameter estimates (such as individual or population parameters) associated with the predictions
For other purposes such as simply obtaining model predictions, we suggest using the predict() method with or without newdata.

The function uses estimateMeanParametersNewdata() to set the population estimates for the individual parameters
taking into account the individual covariates and doses, and estimateIndividualParametersNewdata() to derive individual estimates
by computing the mean of the conditional distributions (type="icpred") or the MAP estimate (type="ipred")
}
\examples{
# TODO
}
