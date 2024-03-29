% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classif.np.r
\name{classif.np}
\alias{classif.np}
\alias{classif.kernel}
\alias{classif.knn}
\title{Kernel Classifier from Functional Data}
\usage{
classif.np(
  group,
  fdataobj,
  h = NULL,
  Ker = AKer.norm,
  metric,
  weights = "equal",
  type.S = S.NW,
  par.S = list(),
  ...
)

classif.knn(
  group,
  fdataobj,
  knn = NULL,
  metric,
  weights = "equal",
  par.S = list(),
  ...
)

classif.kernel(
  group,
  fdataobj,
  h = NULL,
  Ker = AKer.norm,
  metric,
  weights = "equal",
  par.S = list(),
  ...
)
}
\arguments{
\item{group}{Factor of length \emph{n}}

\item{fdataobj}{\code{\link{fdata}} class object.}

\item{h}{Vector of smoothing parameter or bandwidth.}

\item{Ker}{Type of kernel used.}

\item{metric}{Metric function, by default \code{\link{metric.lp}}.}

\item{weights}{weights.}

\item{type.S}{Type of smothing matrix \code{S}. By default \code{S} is
calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}

\item{par.S}{List of parameters for \code{type.S}: \code{w}, the weights.}

\item{\dots}{Arguments to be passed for \code{\link{metric.lp}} o other
metric function and \code{\link{Kernel}} function.}

\item{knn}{Vector of number of nearest neighbors considered.}
}
\value{
\itemize{
\item {fdataobj}{ \code{\link{fdata}} class object.} 
\item {group}{ Factor of length \code{n}.} 
\item {group.est}{ Estimated vector groups}
\item {prob.group}{ Matrix of predicted class probabilities. For each
functional point shows the probability of each possible group membership.}
\item {max.prob}{ Highest probability of correct classification.}
\item {h.opt}{ Optimal smoothing parameter or bandwidht estimated.} 
\item {D}{ Matrix of distances of the optimal quantile distance \code{hh.opt}.}
\item {prob.classification}{ Probability of correct classification by group.}
\item {misclassification}{ Vector of probability of misclassification by
number of neighbors \code{knn}.} 
\item {h}{ Vector of smoothing parameter or bandwidht.} 
\item {C}{ A call of function \code{classif.kernel}.}
}
}
\description{
Fits Nonparametric Supervised Classification for Functional Data.
}
\details{
Make the group classification of a training dataset using kernel or KNN
estimation: \code{\link{Kernel}}.\cr Different types of metric funtions can
be used.
}
\note{
If \code{fdataobj} is a data.frame the function considers the case of
multivariate covariates. \cr \code{\link{metric.dist}} function is used to
compute the distances between the rows of a data matrix (as
\code{\link{dist}} function.
}
\examples{
\dontrun{
data(phoneme)
mlearn <- phoneme[["learn"]]
glearn <- phoneme[["classlearn"]]
h <- 9:19
out <- classif.np(glearn,mlearn,h=h)
summary(out)
head(round(out$prob.group,4))
}

}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
data analysis.} Springer Series in Statistics, New York.

Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on
line at \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/}
}
\seealso{
See Also as \code{\link{predict.classif}}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{classif}
