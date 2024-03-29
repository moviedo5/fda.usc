% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.classif.R
\name{predict.classif}
\alias{predict.classif}
\title{Predicts from a fitted classif object.}
\usage{
\method{predict}{classif}(object, new.fdataobj = NULL, type = "class", ...)
}
\arguments{
\item{object}{Object \code{object} estimated by: k nearest neighbors method
\code{classif.knn}, kernel method \code{classif.kernel}.}

\item{new.fdataobj}{New functional explanatory data of \code{fdata} class.}

\item{type}{Type of prediction ("class or probability of each group
membership").}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
If type="class", produces a vector of predictions.
If type="probs", a list with the following components is returned: 
\itemize{
\item \code{group.pred} the vector of predictions. 
\item \code{prob.group} the matrix of predicted probability by factor level. 
}
}
\description{
Classifier of functional data by kernel method using functional data object
of class \code{classif}.  Returns the predicted classes using a previously trained model.
}
\examples{
\dontrun{
data(phoneme)
mlearn <- phoneme[["learn"]][1:100]
glearn <- phoneme[["classlearn"]][1:100]

#	ESTIMATION knn
out1 <- classif.knn(glearn, mlearn, knn = 3)
summary(out1)

#	PREDICTION knn
mtest <- phoneme[["test"]][1:100]
gtest <- phoneme[["classtest"]][1:100]
pred1 <- predict(out1, mtest)
table(pred1, gtest)
 
#	ESTIMATION kernel 
h <- 2^(0:5)
# using metric distances computed in classif.knn
out2 <- classif.kernel(glearn, mlearn, h = h, metric = out1$mdist)
summary(out2)
#	PREDICTION kernel
pred2 <- predict(out2,mtest)
table(pred2,gtest)
}

}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametricc functional
data analysis.} Springer Series in Statistics, New York.

Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data
Analysis}, 2nd ed., Springer, New York.
}
\seealso{
See also \code{\link{classif.np}} \code{\link{classif.glm}},
\code{\link{classif.gsam}} and \code{\link{classif.gkam}} .
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@udc.es}
}
\keyword{classif}
