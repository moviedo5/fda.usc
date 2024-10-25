#' @title Weighting tools
#' @family performance
#' @description computes inverse probability weighting.
#' @param x A vector of the labels, true class or observed response. Can be \code{numeric, character, or factor}
#' @param type Type of weights.
#' 
#' @aliases weights4class
#' @rdname weights4class
#' @export weights4class
weights4class <- function(x,type=c("equal","inverse")){
  output<-  switch(type[1],
                   equal=rep(1,len=length(x)),
                   inverse=weights_inverse(x)
  )
  output
}

weights_inverse<-function(x){
  tab<-table(x)
  ii <- sum(tab)/tab
  ii<-as.numeric(ii[x])
}

select <-function(fila, vardep.summary, ...) {
  if(length(which(fila==max(fila)))>1)
  {predclass <-names(vardep.summary[which(fila==max(fila))])[
    order(vardep.summary[which(fila==max(fila))],decreasing=TRUE)[1]]
  }
  else{
    predclass<- as.character(names(vardep.summary)[(order(fila,decreasing=TRUE)[1])])
  } 
  predclass
}

prob2classif<-function(yobs,ypred){
  tab <- table(ypred, yobs)
  diag(tab)/colSums(tab)
  
}

