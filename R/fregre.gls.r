#' Fit Functional Linear Model Using Generalized Least Squares
#' 
#' @description This function fits a functional linear model using generalized least
#' squares. The errors are allowed to be correlated and/or have unequal
#' variances.
#' 
#' @param formula a two-sided linear formula object describing the model, with
#' the response on the left of a \code{~} operator and the terms, separated by
#' \code{+} operators, on the right.
#' @param data an optional data frame containing the variables named in
#' \code{model}, \code{correlation}, \code{weights}, and \code{subset}. By
#' default the variables are taken from the environment from which \code{gls}
#' is called.
#' @param correlation an optional \code{\link{corStruct}} object describing the
#' within-group correlation structure. See the documentation of
#' \code{\link{corClasses}} for a description of the available \code{corStruct}
#' classes. If a grouping variable is to be used, it must be specified in the
#' \code{form} argument to the \code{corStruct} constructor. Defaults to
#' \code{NULL}, corresponding to uncorrelated errors.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for \eqn{\beta(t)} parameter estimation.
#' @param rn List of Ridge parameter.
#' @param lambda List of Roughness penalty parameter.
#' @param weights an optional \code{\link{varFunc}} object or one-sided formula
#' describing the within-group heteroscedasticity structure. If given as a
#' formula, it is used as the argument to \code{\link{varFixed}}, corresponding
#' to fixed variance weights. See the documentation on \code{\link{varClasses}}
#' for a description of the available \code{\link{varFunc}} classes. Defaults
#' to \code{NULL}, corresponding to homoscedastic errors.
#' @param subset an optional expression indicating which subset of the rows of
#' \code{data} should be used in the fit. This can be a logical vector, or a
#' numeric vector indicating which observation numbers are to be included, or a
#' character vector of the row names to be included.  All observations are
#' included by default.
#' @param method a character string.  If \code{"REML"} the model is fit by
#' maximizing the restricted log-likelihood.  If \code{"ML"} the log-likelihood
#' is maximized.  Defaults to \code{"REML"}.
#' @param control a list of control values for the estimation algorithm to
#' replace the default values returned by the function
#' \code{\link{glsControl}}.  Defaults to an empty list.
#' @param verbose an optional logical value. If \code{TRUE} information on the
#' evolution of the iterative algorithm is printed. Default is \code{FALSE}.
#' @param criteria GCCV criteria, see \code{\link{GCCV.S}}.
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this methodl.
#' @return an object of class \code{"gls"} representing the functional linear
#' model fit. Generic functions such as \code{print}, \code{plot}, and
#' \code{summary} have methods to show the results of the fit.\cr 
#' See \code{\link{glsObject}} for the components of the fit. The functions
#' \code{\link{resid}}, \code{\link{coef}} and \code{\link{fitted}}, can be
#' used to extract some of its components.\cr 
#' Beside, the class(z) is "gls", "lm" and "fregre.lm" with the following
#' objects: 
#' \itemize{
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{Vp}{ Estimated covariance matrix for the parameters.} 
#' \item \code{lambda}{ A roughness penalty.}
#' \item \code{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.}
#' \item \code{basis.b}{ Basis used for beta parameter estimation.} 
#' \item \code{beta.l}{ List of estimated beta parameter of functional covariates.} 
#' \item \code{data}{ List that containing the variables in the model.} 
#' \item \code{formula}{ formula used in ajusted model.} 
#' \item \code{formula.ini}{ formula in call.}
#' \item \code{W}{ inverse of covariance matrix} 
#' \item \code{correlation}{ See glsObject for the components of the fit. }
#' }
#' @references Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and
#' Dominguez, A. Predicting seasonal influenza transmission using Functional
#' Regression Models with Temporal Dependence. arXiv:1610.08718.
#' \url{https://arxiv.org/abs/1610.08718}
#' @keywords regression models
#' @examples
#' \dontrun{ 
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d2<-fdata.deriv(x,nderiv=)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' 
#' # plot the response
#' plot(ts(tecator$y$Fat))
#' 
#' nbasis.x=11;nbasis.b=7
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' basis.x=list("x.d2"=basis1)
#' basis.b=list("x.d2"=basis2)
#' ldata=list("df"=dataf,"x.d2"=x.d2)
#' res.gls=fregre.gls(Fat~x.d2,data=ldata, correlation=corAR1(),
#'                    basis.x=basis.x,basis.b=basis.b)
#' summary(res.gls)                   
#' }
#' @rdname fregre.gls
#' @export
fregre.gls<-function(formula,data,correlation = NULL,
basis.x=NULL,basis.b=NULL,rn,lambda,weights=NULL,subset, 
method = c("REML", "ML"),control = list(),verbose = FALSE, criteria="GCCV1",...){
 tf <- terms.formula(formula)

 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))                                                    
 vnf=intersect(terms,names(data$df))
 vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
 n<-length(data[["df"]][,response])
 #print(response) 
# XX=data.frame(data[["df"]][,response],weights)
 XX=data.frame(data[["df"]][,response])
 iresp<-which(names(data[["df"]])==response)
# namxx=names(XX)=c(response,"weights")
 namxx=names(XX)=c(response)

# if (!is.null(correlation)) {
#        groups <- getGroupsFormula(correlation)
#    } 
#    else groups <- NULL
#    cls1<-getCovariate(correlation,data=data$df)
#     print(cls1)    
#    print(groups[[2]])  
  
 if (length(vnf)>0) {          
#print(paste("Functional covariate:",vnf))
# XX=data.frame(XX,data[["df"]][,c(vnf)])    # si es facotr poner k-1 variables dummies
# names(XX)=c(namxx,vnf)
 for ( i in 1:length(vnf)){
#     print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#     else pf <- paste(pf, terms[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
#  contrasts = NULL
    mf<-as.data.frame(model.matrix(formula(pf),data$df))
    vnf2<-names(mf)[-1]
    pf <- rf <- paste(response, "~", sep = "")
    for ( i in 1:length(vnf2))  pf<-paste(pf, "+", vnf2[i], sep = "")
    XX <- data.frame(XX,mf)
}
else {
    if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
      XX <- data.frame(XX)
     }
    else { 
    XX <- data.frame(XX,model.matrix(formula(paste(pf, "1")),data$df))
#    names(XX)[3]<-"(Intercept)"
    names(XX)[2]<-"(Intercept)"
    }
}
#print(paste("Functional covariate:",vfunc))
if (missing(rn))    {    rn0=FALSE;                    rn=list()}
else rn0<-TRUE
if (missing(lambda))    {    lambda0=FALSE;                    lambda=list()}
else lambda0<-TRUE
mat<-rep(0,len=ncol(XX)-2)
imat2<-ncol(XX)-2
mat2<-diag(0,nrow=imat2)
if (length(vfunc)>0) {
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:length(vfunc)) {
	if(is(data[[vfunc[i]]],"fdata")){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      np<-length(tt)
      fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
#          else int<-1:basis.x[[vfunc[i]]]$nbasis
#                  nam[[vfunc[i]]]<-  basis.x[[vfunc[i]]]$names
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
    	    x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
#          Theta = getbasismatrix(tt,basis.b[[vfunc[i]]])
#          Psi = getbasismatrix(tt,basis.x[[vfunc[i]]])
#          J = t(Psi) %*% Theta
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = data.frame(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
           nbasisb<- basis.b[[vfunc[i]]]$nbasis
         R=diag(0,ncol= nbasisb,nrow=nbasisb)
         R<-eval.penalty(basis.b[[vfunc[i]]],Lfdobj = int2Lfd(2))
        mat2<-diag(0,nrow=imat2+nbasisb)
        if (!is.null(lambda[[vfunc[i]]]))  {
        mat2<-diag(0,nrow=imat2+nbasisb)
        mat2[(imat2+1):(imat2+nbasisb),(imat2+1):(imat2+nbasisb)]<-lambda[[vfunc[i]]]*R
        imat2<-imat2+nbasisb
        }     
				}
      else {
        basis<-basis.x[[vfunc[i]]]
        l<-basis$l
        vs <- t(basis$basis$data)           
        basis$x<-basis$x[,l,drop=FALSE]
        Z<-basis$x
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis$basis$data),sep ="")
        XX =data.frame(XX,Z)
        vs.list[[vfunc[i]]]=basis$basis
        mean.list[[vfunc[i]]]=basis$mean
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        if (!is.null(rn[[vfunc[i]]]))  {
           mat<-c(mat,rep(rn[[vfunc[i]]],len=length(l)))
           }     
        else          mat<-c(mat,rep(0,len=length(l)))
        mat2<-diag(0,nrow=imat2+length(l))
        if (!is.null(lambda[[vfunc[i]]]))  {
        mat2<-diag(0,nrow=imat2+length(l))
        R<-P.penalty(1:length(l))
        mat2[(imat2+1):(imat2+length(l)),(imat2+1):(imat2+length(l))]<-R
        imat2<-imat2+length(l)
        }                    
   }
  }
 	else {
 		if(class(data[[vfunc[i]]])[1]=="fd"){
      fdat<-data[[vfunc[i]]]
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)
         basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
         l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
         rangeval=fdat$basis$rangeval)
      else           if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp2=FALSE
      if (bsp1 & bsp2) {
          r=fdat[[2]][[3]]
#          tt = seq(r[1], r[2], len = length(fdat[[3]]$time))
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          mean.list[[vfunc[i]]]<-mean.fd(x.fd)
          x.fd<-center.fd(x.fd)
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = data.frame(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
      else {
        l<-ncol(basis.x[[vfunc[i]]]$scores)
        vs <- basis.x[[vfunc[i]]]$harmonics$coefs
        Z<-basis.x[[vfunc[i]]]$scores
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
        XX = data.frame(XX,Z)
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
         }
    }
   else stop("Please, enter functional covariate")
   }
  }  }
 if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    y<-XX[,1]    
    ycen = y - mean(y)
    scores<-as.matrix(XX[,-c(1)])     
    
#    scores<-as.matrix(XX)     
#    W<-diag(weights)  
    W<-diag(n)
    if (!rn0 & !lambda0) {
    pf<-formula(pf)
#    print(pf)
#    print(class(pf))
#    print(names(XX)) 
#     print("antes gls")    
#     XX=data.frame(XX,data[["df"]][,-iresp])
     XX=data.frame(XX,data[["df"]])
#     print(names(XX))
# print(correlation)
#      print(XX)
#print(weights)
#print("a")
      #z=glse(model=pf,data=XX,correlation=correlation,method =method,control=control,verbose=verbose,weights=weights,...) 
      z=nlme::gls(model=pf,data=XX,correlation=correlation,method =method,control=control,verbose=verbose,weights=weights,...) 
#print("b")
      e<-z$residuals
      coef<-z$coefficients

       mat0<-diag(0,ncol(scores))
#       if (lambda0)              mat0<-mat2  
#       if (rn0)              mat0<-diag(mat)

      S<-solve(t(scores)%*%W%*%scores)          
#       W<-diag(weights)%*%W
      S<-solve(t(scores)%*%W%*%scores+mat0)
# print(dim(S))      
# print(dim(W))      
#print((H))      
###       S<-diag(coef(summary(z))[,2]) # como conseguir S
#     print("no penalizado")

      
  if (!is.null(correlation)) {  
    par.gls=list(correlation=correlation)
    XX<-data[["df"]]
    XX<-data.frame(XX,e=e)
    par.gls$data<-XX
# print("oo0")
    z$out.gls<-corSigma(par.gls)
# print("oo1")    
    W<-z$out.gls$W
#    par.CV$W<-Wend
#    print("peta ak")
#    gcv.it[i]<- do.call(tcv,par.CV)   #si carmack se hace solve(W)
                                     #devolver df como atributo degcv.it    
#      print(W[1:3,1:4])

   }  
        S<-solve(t(scores)%*%W%*%scores)                   
        Cinv<-S%*%t(scores)%*%W   
        z$H<-scores%*%Cinv    
        z$S<-S
        z$gcv<-GCCV.S(y,z$H,criteria=criteria,W)                                                
# print(GCCV.S(y,z$H,criteria=criteria))
        df<-attributes(z$gcv)$df                        
        rdf<-n-df      
        z$sr2 <- sum(e^2)/ rdf
# print(z$sr2)        
        z$r2 <- 1 - sum(e^2)/sum(ycen^2)
        z$r2.adj<- 1 - (1 - z$r2) * ((n -    1)/ rdf)
        z$rank <- df
        z$df.residual <-  rdf    
 #    print("no penalizado3")
        
      class(z)<-c(class(z),"fregre.lm")
      }      
    else {
#print("penalized version")    
       if (!is.null(correlation)) stop("gls with penalization parameter not implemented yet")
       if (lambda0) { S<-solve(t(scores)%*%W%*%scores+mat2)           }             
       if (rn0) {
             mat<-diag(mat)
             S<-solve(t(scores)%*%W%*%scores+mat)       
       }      
       if (rn0 & lambda0) warning("Only ridge penalization is done by rn argument (lambda argument is ignored)")               
#      ddd<-t(scores)%*%W  
#      S<-solve(t(scores)%*%W%*%scores+mat+mat2)       #incluir pesos solve(W)
      Cinv<-S%*%t(scores)%*%W              
      H<-scores%*%Cinv      
      coefs<-Cinv%*%XX[,1]
      z<-list()      
      e<-z$residuals
      yp<-z$fitted.values
#      z$fitted.values<-yp<-drop(scores%*%coefs)
#      e<-z$residuals<-XX[,1]- z$fitted.values      
################################################################################
################################################################################
################################################################################           
#      H<-scores%*%Cinv
      coefs<-drop(coefs)
#      z$coefficients<-coefs
      z$mean.list<-mean.list
#      z$df.residual<-n-df
 #      df<-n-z$df.residual
#      z$H<-H
      z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)       
      if  (class(basis.x[[vfunc[1]]])=="basisfd") {
        z$call[[1]] = "fregre.basis"
        }
       else {
        if  (basis.x[[vfunc[1]]]$type=="pc")  z$call[[1]] = "fregre.pc"
        if  (basis.x[[vfunc[1]]]$type=="pls")  z$call[[1]] = "fregre.pls"        
        }             
    class(z) <- c("fregre.fd","fregre.lm")
#   print("penalizado")
        z$gcv<-GCCV.S(y,z$H,criteria=criteria,W)   
# print(z$gcv)        
# print(GCCV.S(y,z$H,criteria=criteria,diag(n)   ))
        df<-attributes(z$gcv)$df                          
        rdf<-n-df
        sr2 <- sum(e^2)/ rdf
        r2 <- 1 - sum(e^2)/sum(ycen^2)
        r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
        z$rank <- df
        z$df.residual <-  rdf    
#    GCV <- sum(e^2)/(n - df)^2
        z$residuals <- drop(e)
        z$fitted.values <- yp
        z$y <- y

        Z=cbind(rep(1,len=n),Z)
        colnames(Z)[1] = "(Intercept)"
        std.error = sqrt(diag(S) *sr2)
        t.value = coefs/std.error
        p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
        coefficients <- cbind(coefs, std.error, t.value, p.value)
        colnames(coefficients) <- c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)")
        z$coefs<-coefficients  
        z$sr2<-sum(e^2)/z$df.residual  
        class(z) <- "lm"
}       
#    z$call<-z$call[1:2]
for (i in 1:length(vfunc)) {
 if (bsp1) {
 beta.l[[vfunc[i]]]=fd(z[["coefficients"]][name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
 }
 else{
	if(class(data[[vfunc[i]]])[1]=="fdata"){
#     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est$data<-colSums(beta.est$data)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
           if  (basis.x[[vfunc[i]]]$type=="pls") {
             if (basis.x[[vfunc[i]]]$norm)  {
              sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
            }  
        
     beta.l[[vfunc[i]]]<-beta.est     
     }
 else {
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(vs.list[[vfunc[i]]])
#     beta.est<-apply(beta.est,2,sum)
     beta.est<-colSums(beta.est)
     beta.l[[vfunc[i]]]<-fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis)
      }
}
}
 z$Vp=z$sr2*S
 z$beta.l=beta.l
 z$formula=pf
 z$mean=mean.list
 z$formula.ini=formula
 z$basis.x=basis.x
 z$basis.b=basis.b
 z$JJ<-JJ
 z$XX=XX
 z$W<-W
 z$data<-data
 z$fdataobj<-data[[vfunc[1]]]
 z$rn<-rn0
 z$lambda<-lambda0
 z$vs.list=vs.list   
 z$correlation<-correlation

 class(z)<-c("fregre.gls","gls","lm")
 z
}     

