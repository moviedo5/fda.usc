#' @rdname fanova.RPm
#' @export
summary.fanova.RPm<-function (object, ndec=NULL,...) {
 if (is.null(ndec)) ndec=5
 if (is(object,"fanova.RPm")) {
    cat("     - SUMMARY fanova.RPm - \n")
    cat("\n p-value for Bonferroni method \n" )
    print(round(object$p.Bonf,ndec))
    cat("\n  p-value for False Discovery Rate method \n")
     print(round(object$p.FDR,ndec))
    if (!is.null(object$p.Boot)){
           cat("\n p-value for Bootstrap method \n" )
           print(round(object$p.Boot,ndec))
     }
 }
 if (class(object)=="fanova.hetero") {
     cat("\n  - SUMMARY ANOVA HETEROCEDASTHIC - \n")
          print(object)
  }
cat("\n")
output<-object
}


