#' @rdname  influence.fregre.fd
#' @aliases influence.fregre.fd
#' @note influence.fdata deprecated.
#' @title Functional influence measures
#' 
#' @description Once estimated the functional regression model with scalar response,
#' influence.fregre.fd function is used to obtain the functional influence measures.
#' @param  model \code{fregre.pc}, \code{fregre.basis} or \code{fregre.basis.cv} object.
#' @param \dots Further arguments passed to or from other methods.

#' @details  Identify influential observations in the functional linear model in which 
#' the predictor is functional and the response is scalar.
#' Three statistics are introduced for measuring the influence:   Distance Cook Prediction
#'  \code{DCP}, Distance Cook Estimation \code{DCE} and Distance
#'   \eqn{\mbox{pe}\tilde{\mbox{n}}\mbox{a} }{} \code{DP} respectively.
#' @return Return:
#' \itemize{
#'  \item \code{DCP}{ Cook's Distance for Prediction.}
#'  \item \code{DCE}{ Cook's Distance for Estimation.}
#'  \item \code{DP}{ \eqn{\mbox{Pe}\tilde{\mbox{n}}\mbox{a's} }{} Distance.}
#'  }
#' @references
#' Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2010). \emph{Measures of influence for the functional linear model with scalar response}. Journal of Multivariate Analysis 101, 327-339.
#' 
#' Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
#' Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso See Also as:  \code{\link{fregre.pc}}, \code{\link{fregre.basis}}, 
#'  \code{\link{influence_quan}}
#' @examples
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata[1:129]
#' y=tecator$y$Fat[1:129]
#' 
#' res1=fregre.pc(x,y,1:5)  
#' # time consuming
#' res.infl1=influence(res1)  
#' res2=fregre.basis(x,y)  
#' res.infl2=influence(res2)  
#' 
#' res<-res1
#' res.infl<-res.infl1
#' mat=cbind(y,res$fitted.values,res.infl$DCP,res.infl$DCE,res.infl$DP)
#' colnames(mat)=c("Resp.","Pred.","DCP","DCE","DP")
#' pairs(mat)
#' }
#' @keywords outliers 
#' 
#' @export 
influence.fregre.fd<-function(model,...){
if (!is.fdata(model$fdataobj)) fdataobj=fdata(model$fdataobj)
else fdataobj<-model$fdataobj
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names<-fdataobj[["names"]]
y=model$y
n<-nrow(x)
J<-ncol(x)
fitted.values=model$fitted.values
sr2=model$sr2
H=model$H
dist.cook.for <- array(NA,dim=c(n,1))
dist.cook.est2<-dist.cook.est <- array(NA,dim=c(n,1))
dist.pena <- array(NA,dim=c(n,1))
S <- array(NA,dim=c(n,n))
if (model$call[[1]]=="fregre.pc") {
   betas<-beta.est<-model$beta.est #/(ncol(fdata)-1)
   l=model$l
   lambdas=model$fdata.comp$d^2
   for (i in 1:n){
     oo <- fregre.pc(fdataobj[-i,],y[-i],l,lambda=model$lambda,P=model$P,weights=model$weights[-i])
#     G <- oo$svd.fdata$x[,l]
#     I <- diag(1/((n-1)*oo$svd.fdata$lambdas[l]),ncol=kn)
     best <- oo$beta.est #/(ncol(fdata)-1)
     aest <- oo$a.est
#    ypi <-   aest * rep(1,n) + fdata %*% best
     ypi<-predict(oo,fdataobj)
     S[i,] <- t(fitted.values-ypi)
     dist.cook.for[i] <- t(S[i,]) %*% S[i,] / sr2
     dist.cook.est[i] <- sum((beta.est$data-best$data)^2)/(sr2/n*(sum(1/lambdas[l])))
     bb<-beta.est-best
     betas<- c(betas,best)
#    a<-sum((beta.est$data-best$data)^2)
       b<-norm.fdata(fdata(beta.est$data-best$data,tt,rtt))
#    dd<-a/b
     dist.cook.est[i] <- as.numeric(norm.fdata(bb))/(sr2/n*(sum(1/lambdas[l])))
    }
    betas$data<-betas$data[-1,]
}
if (model$call[[1]]=="fregre.basis" || model$call[[1]]=="fregre.basis.cv") {
beta.est<-model$beta.est #/(ncol(fdata)-1)
betas<-b2<-eval.fd(tt,beta.est)#*diff(range(tt))
for (i in 1:n){
    fdata_i<-fdataobj[-i,]
    oo <- fregre.basis(fdata_i,y[-i],basis.x=model$basis.x.opt,basis.b=model$basis.b.opt,
    lambda=model$lambda.opt, Lfdobj=model$Lfdobj,weights=model$weights[-i])
    best <- oo$beta.est
    aest <- oo$a.est
    ypi<-predict(oo,fdataobj)
    S[i,] <- t(fitted.values-ypi)
    dist.cook.for[i] <- t(S[i,]) %*% S[i,] / sr2
    b1<-eval.fd(tt,best)
    bb<-best-beta.est
#    dist.cook.est[i] <- sum((b2-b1)^2)/(sr2/n)
     dist.cook.est[i] <- norm.fd(bb)/(sr2/n)
     betas<- cbind(betas,b1)
    }
   betas<- t(betas)
   betas<-fdata(betas[-1,],tt,rtt,names=list(main="beta CV"))
    }
for (i in 1:n){dist.pena[i] <- t(S[,i]) %*% S[,i] / (sr2 * diag(H)[i])}
return(list("H"=diag(H),"DCP"=dist.cook.for,"DCE"=dist.cook.est,"DP"=dist.pena,
betas=betas))
}
