#' @name summary.classif
#' @title Summarizes information from kernel classification methods.
#' 
#' @description Summary function for \code{\link{classif.knn}} or \code{\link{classif.kernel}}.
#' 
#' @details \code{object} from \code{\link{classif.knn}} or \code{\link{classif.kernel}}
#' 
#' @aliases summary.classif print.classif
#' @param object Estimated by kernel classification.
#' @param x Estimated by kernel classification.
#' @param digits how many significant digits are to be used for numeric and complex x.
#' @param \dots Further arguments passed to or from other methods.
#' @return Shows:
#' \itemize{ 
#' \item -Probability of correct classification by group \code{prob.classification}.
#' \item -Confusion matrix between the theoretical groups and estimated groups.
#' \item -Highest probability of correct classification \code{max.prob}. 
#' }
#'  If the object is returned from the function \code{\link{classif.knn}}
#' \itemize{ 
#' \item -Vector of probability of correct classification by number of neighbors \code{knn}.
#' \item  -Optimal number of neighbors: \code{knn.opt}. 
#'  }
#'  If the object is returned from the function: \code{\link{classif.kernel}}
#' \itemize{ 
#' \item -Vector of probability of correct classification by banwidth \code{h}.
#' \item  -Functional measure of closeness (optimal distance, \code{h.opt}). 
#  \item {object}{ Estimated by kernel classification.}
#'  } 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{classif.knn}},
#' \code{\link{classif.kernel}} \cr and \code{\link{summary.classif}}
#' @keywords print
#' @examples
#' \dontrun{ 
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-phoneme[["classlearn"]]
#' out=classif.knn(glearn,mlearn,knn=c(3,5,7))
#' summary(out)
#' out2=classif.kernel(glearn,mlearn,h=2^(0:5))
#' summary(out2)
#' }
#' @rdname summary.classif
#' @export
summary.classif<-function (object, ...)
{
 cat("     - SUMMARY - \n")
 cat("\n-Probability of correct classification by group (prob.classification):\n")
 print(object$prob.classification)
 cat("\n-Confusion matrix between the theoretical groups (by rows)
  and estimated groups (by column) \n")
 print(table(object$group,object$group.est))
 if (object$C[[1]]=='classif.np'){
    if (object$ty=="S.NW") {
   cat("\n-Vector of probability of correct classification
    by banwidth (h):\n")
    print(round(1-object$gcv,4))
#   cat("\n-Functional measure of closeness (optimal distance, h.opt):\n")
#   print(round(object$h.opt,4))

cat("\n-Optimal bandwidth: h.opt=",object$h.opt,"with highest probability of
correct classification: max.prob=",object$max.prob,"\n")
  }  }
 if (object$C[[1]]=='classif.np'){
   if (object$ty=="S.KNN") {
   cat("\n-Vector of probability of correct classification
   by number of neighbors (knn):\n")
    print(round(1-object$gcv,4))
    cat("\n-Optimal number of neighbors: knn.opt=",object$h.opt,
    "\nwith highest probability of correct classification max.prob= ",
    object$max.prob,"\n")
    }}
    # if (object$C[[1]]=='classif.gkam' | object$C[[1]]=='classif.gkam2boost' |
    #  object$C[[1]]=='classif.gsam2boost' | object$C[[1]]=='classif.gsam'
    #  | object$C[[1]]=='classif.glm'| object$C[[1]]=='classif.glm2boost'
    #  | object$C[[1]]=='classif.bootstrap')
  if (!is.null(object$max.prob))  {
   cat("\n-Probability of correct classification: ",round(object$max.prob,4),"\n")
    }
    if (object$C[[1]]=='classif.rpart'|object$C[[1]]=='classif.rpart2boost'){
   cat("\n-Probability of correct classification: ",round(object$max.prob,4),"\n")}
   if (object$C[[1]]=='classif.DD'){
   cat("\n-Probability of correct classification: ",round(1-object$misclassification,4),"\n")}
cat("\n")
output<-object
}

#' @rdname summary.classif
#' @export 
print.classif<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\n-Call:\n", deparse(x$C), "\n", sep = "")
  #   cat("\n-Probability of correct classification by group (prob.classification):\n")
  #   print(x$prob.classification)
  if (x$C[[1]]=='classif.knn'){
    cat("\n-Optimal number of neighbors: knn.opt=",x$knn.opt,
        "\nwith highest probability of correct classification max.prob=",
        x$max.prob,"\n")}
  else {if (x$C[[1]]=='classif.kernel'| x$C[[1]]=='classif.np') {
    cat("\n-Optimal bandwidth: h.opt=",x$h.opt,"with highest probability of
   correct classification: max.prob=",x$max.prob,"\n")
  }
    else {
      # if  (x$C[[1]]=='classif.gsam' | x$C[[1]]=='classif.gsam2boost'|
      #          x$C[[1]]=='classif.glm' |x$C[[1]]== 'classif.glm2boost' |
      #          x$C[[1]]=='classif.gkam' | x$C[[1]]=='classif.gkam2boost'|
      #          x$C[[1]]=='classif.rpart' |   x$C[[1]]=='classif.rpart2boost')
      if (!is.null(x$max.prob))  {
      cat("\n-Probability of correct classification: ",round(x$max.prob,4),"\n")
    }
      if (x$C[[1]]=='classif.DD'){
        cat("\n-Probability of correct classification: ",round(1-x$misclassification,4),"\n")}
    }
  }
  cat("\n")
  invisible(x)
}

