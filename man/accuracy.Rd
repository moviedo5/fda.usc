% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accuracy.R
\name{accuracy}
\alias{accuracy}
\alias{cat2meas}
\alias{tab2meas}
\alias{pred2meas.}
\alias{pred.MSE}
\alias{pred.RMSE}
\alias{pred.MAE}
\alias{pred2meas}
\title{Performance measures for regression and classification models}
\usage{
cat2meas(yobs, ypred, measure = "accuracy", cost = rep(1, nlevels(yobs)))

tab2meas(tab, measure = "accuracy", cost = rep(1, nrow(tab)))

pred.MSE(yobs, ypred)

pred.RMSE(yobs, ypred)

pred.MAE(yobs, ypred)

pred2meas(yobs, ypred, measure = "RMSE")
}
\arguments{
\item{yobs}{A vector of the labels, true class or observed response. Can be \code{numeric}, \code{character}, or \code{factor}.}

\item{ypred}{A vector of the predicted labels, predicted class or predicted response. Can be \code{numeric, character, or factor}.}

\item{measure}{Type of measure, see \code{details} section.}

\item{cost}{Cost value by class (only for input factors).}

\item{tab}{Confusion matrix (Contingency table: observed class by rows, predicted class by columns).}
}
\description{
\code{\link{cat2meas}} and \code{\link{tab2meas}} calculate the measures for a multiclass classification model.\cr
 \code{\link{pred2meas}} calculates the measures for a regression model.
}
\details{
\itemize{
\item \code{\link{cat2meas}} compute \eqn{tab=table(yobs,ypred)} and calls \code{\link{tab2meas}} function.
\item \code{\link{tab2meas}} function computes the following measures (see \code{measure} argument) for a binary classification model:
\itemize{
\item \code{accuracy}: Proportion of correct predictions.
  \eqn{\frac{TP + TN}{TP + TN + FP + FN}}
\item \code{sensitivity, TPrate, recall}: True Positive Rate or recall.
  \eqn{\frac{TP}{TP + FN}}
\item \code{precision}: Positive Predictive Value.
  \eqn{\frac{TP}{TP + FP}}
\item \code{specificity, TNrate}: True Negative Rate.
  \eqn{\frac{TN}{TN + FP}}
\item \code{FPrate}: False Positive Rate.
  \eqn{\frac{FP}{TN + FP}}
\item \code{FNrate}: False Negative Rate.
  \eqn{\frac{FN}{TP + FN}}
\item \code{Fmeasure}: Harmonic mean of precision and recall.
  \eqn{\frac{2}{\frac{1}{\text{recall}} + \frac{1}{\text{precision}}}}
\item \code{Gmean}: Geometric Mean of recall and specificity.
  \eqn{\sqrt{\left(\frac{TP}{TP + FN}\right) \cdot \left(\frac{TN}{TN + FP}\right)}}
\item \code{kappa}: Cohen's Kappa index.
  \eqn{Kappa = \frac{P_o - P_e}{1 - P_e}} where \eqn{P_o} is the proportion of observed agreement, 
  \eqn{\frac{TP + TN}{TP + TN + FP + FN}}, and \eqn{P_e} is the proportion of agreement expected by chance, 
  \eqn{\frac{(TP + FP)(TP + FN) + (TN + FN)(TN + FP)}{(TP + TN + FP + FN)^2}}.
\item \code{cost}: Weighted accuracy, calculated as 
  \eqn{\frac{\sum (\text{diag(tab)} / \text{rowSums(tab)} \cdot \text{cost})}{\sum(\text{cost})}}
\item \code{IOU}: Mean Intersection over Union.
  \eqn{\frac{TP}{TP + FN + FP}}
\item \code{IOU4class}: Intersection over Union by class level.
  \eqn{\frac{TP}{TP + FN + FP}}#' }
\item \code{\link{pred2meas}} function computes the following  measures of error, usign the \code{measure} argument, for observed and predicted vectors:
\itemize{
\item \code{MSE}: Mean squared error, \eqn{\frac{\sum{(ypred- yobs)^2}}{n} }{\sum (ypred- yobs)^2 /n }.
\item \code{RMSE}: Root mean squared error \eqn{\sqrt{\frac{\sum{(ypred- yobs)^2}}{n} }}{\sqrt(\sum (ypred- yobs)^2 /n )}.
\item \code{MAE}: Mean Absolute Error, \eqn{\frac{\sum |yobs - ypred|}{n}}{\sum |yobs - ypred| /n}.
}
}
}
\seealso{
Other performance: 
\code{\link{weights4class}()}
}
\concept{performance}
