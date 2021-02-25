#' @rdname fanova.hetero
#' @aliases fanova.hetero anova.hetero
#' @title ANOVA for heteroscedastic data
#' @note anova.hetero deprecated
#' @description Univariate ANOVA for heteroscedastic data.
#' 
#' @details This function fits a univariate analysis of variance model and allows
#' calculate special contrasts defined by the user.  The list of special
#' contrast to be used for some of the factors in the formula.  Each matrix of
#' the list has \code{r} rows and \code{r-1} columns.
#' 
#' The user can also request special predetermined contrasts, for example using
#' \code{\link{contr.helmert}}, \code{\link{contr.sum}} or
#' \code{\link{contr.treatment}} functions.
#' 
#' @param object A data frame with dimension (\code{n} x \code{p+1}).  In the
#' first column contains the \code{n} response values and on the following
#' \code{p} columns the explanatory variables specified in the formula.
#' @param formula as \link[stats]{formula}.
#' @param pr If TRUE, print intermediate results.
#' @param contrast List of special contrast to be used, by default no special
#' contrasts are used (\code{contrast}=\code{NULL}).
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Return: 
#' \itemize{
#' \item{ans}{ A list with components including: the Beta estimation \code{Est},
#'  the factor degrees of freedom \code{df1}, the residual degrees of freedom
#'   \code{df2} and \code{p-value} for each factor. }
#' \item{contrast}{ List of special contrasts.}
#' }
#' 
#' @note It only works with categorical variables.
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fanova.RPm}}
#' @references Brunner, E., Dette, H., Munk, A. \emph{Box-Type Approximations
#' in Nonparametric Factorial Designs.} Journal of the American Statistical
#' Association, Vol. 92, No. 440 (Dec., 1997), pp. 1494-1502.
#' @keywords anova
#' @examples
#' \dontrun{
#' data(phoneme)
#' ind=1 # beetwen 1:150
#' fdataobj=data.frame(phoneme$learn[["data"]][,ind])
#' n=dim(fdataobj)[1]
#' group<-factor(phoneme$classlearn)
#' 
#' #ex 1: real factor and random factor
#' group.rand=as.factor(sample(rep(1:3,n),n))
#' f=data.frame(group,group.rand)
#' mm=data.frame(fdataobj,f)
#' colnames(mm)=c("value","group","group.rand")
#' out1=fanova.hetero(object=mm[,-2],value~group.rand,pr=FALSE)
#' out2=fanova.hetero(object=mm[,-3],value~group,pr=FALSE)
#' out1
#' out2
#' 
#' #ex 2: real factor, random factor and  special contrasts
#' cr5=contr.sum(5)  #each level vs last level
#' cr3=c(1,0,-1)			#first level vs last level
#' out.contrast=fanova.hetero(object=mm[,-3],value~group,pr=FALSE,
#' contrast=list(group=cr5))
#' out.contrast
#'}     
#' @export 
fanova.hetero=function(object=NULL,formula,pr=FALSE,contrast=NULL,...){
if (pr) print("INI fanova.hetero")
  if (is.data.frame(object)) data=data.frame(object)
  else if (is.fdata(object)) data=data.frame(object[["data"]])
  nrow=nrow(data);ncol=ncol(data)
  fdata2=intercambio(data) #para evitrar transponer
  Terms <- if (missing(fdata2)) terms(formula, "Error")
  else terms(formula, "Error", fdata2 = fdata2)
  indError <- attr(Terms, "specials")$Error
   if (length(indError) > 1) stop(sprintf(ngettext(length(indError),
   "there are %d Error terms: only 1 is allowed",
  "there are %d Error terms: only 1 is allowed"), length(indError)), domain = NA)
  mf=model.frame(formula,fdata2)
  mfo=model.frame(formula,data) ####
  mt=model.matrix(formula,mf)
  ff=attr(Terms,"factors")
  ffcol=colnames(ff)
  ffrow=rownames(ff)
  ff= as.matrix(ff[-1,],nrow=nrow(ff[-1,]))
  rownames(ff)=ffrow[-1]
  colnames(ff)=ffcol
  nombres=attr(Terms,"term.labels")
  uniq=apply(mf,2,function(x){length(unique(x))})
  lis=sapply(mf,is.factor)
  fact=as.vector(lis)
#  fact=uniq<10
  nfactors= length(uniq[fact])
  nterms=length(Terms)
  aaa=dim(mf)[2]-1
  dddd=(dim(mf)[2]-1)
#  if (nfactors!=(dim(mf)[2]-1)) stop("Not Yet Implemented the use of
# covariables or number of levels great than 10")
 if (nfactors!=(dim(mf)[2]-1)) stop("Not Yet Implemented the use of
 covariables")
  if (is.null(contrast)) {
     ans=matrix(NA,ncol=4,nrow=length(nombres))
     rownames(ans)=nombres}
  else {
            contrast2=intercambio.l(contrast)
            b=length(contrast)
            cnombres=ncontrast=rep(0,len=b)
            tnombres=rep(FALSE,len=length(ncontrast))
            tgroups=rep(FALSE,len=(length(ffcol)+1))
            bb=length(uniq)
            for (i in 1:b)    {
                       a=which(colnames(mfo)==names(contrast)[i])
                       tgroups[a]=tnombres[a]=TRUE
                       if (is.vector(contrast[[i]])) {
                          ncontrast[a]=1
                          contrast[[i]]=matrix(contrast[[i]],ncol=1)
                          }
                       else ncontrast[a]=ncol(contrast[[i]])
                       names(ncontrast)[a]=names(contrast[i])
                       }
            name.contr=rep(NA,len=sum(ncontrast))
            j=1;ji=1;jk=1
            for (i in 1:length(ncontrast))    {
            if (tgroups[ji])  {
                name.contr[j:(j+ncontrast[i]-1)]=paste("C",j:(j+ncontrast[i]-1),".",
                names(contrast[jk]),sep="")
                colnames(contrast[[jk]])=name.contr[j:(j+ncontrast[i]-1)]
              j=j+ncontrast[i];jk=jk+1
              }
              ji=ji+1
              }
    nombres.esp=c(nombres,name.contr)
    ans=matrix(NA,nrow=length(nombres.esp),ncol=4)
    rownames(ans)=nombres.esp          }
  colnames(ans)=c("Est.","df1","df2","p.value")
  nlev=uniq[fact]
  if (pr) print(paste("Levels:",nlev))
  mf2=intercambio(mf)   #para evitar transponer
  n=table(mf2[,fact])
  m=tapply(mf2[,1],mf2[,fact],mean)
  sig=tapply(mf2[,1],mf2[,fact],var)
  sig=sig*(n-1)/n^2
  max.vecm=length(m)
  vecm=matrix(as.vector((m)),nrow=1)
  ssig=as.vector((sig))
  if (pr) {print("Means");print(m)}
  if (pr) {print("Variances");print(sig)}
  vecDelta=as.vector(1/(n-1))
  Delta=diag(vecDelta)
  SN=sum(n)*diag(ssig)
  if (pr) {print("Delta");print(Delta);print("SN");print(SN)}
  for (i in 1:length(nombres)){
      MM=1
      for (k in 1:nrow(ff)){
          J=matrix(1,ncol=nlev[k],nrow=nlev[k])/nlev[k]
          if (ff[k,i]==1) P=diag(1,nlev[k])-J else P=J
          MM=kronecker(MM,P)      #
          }
      DM=diag(diag(MM))
      if (pr) {print(paste("MM",nombres[i]));print(MM)}
      FNM=sum(n)*vecm%*%MM%*%t(vecm)/fdata.trace(DM%*%SN)
      f1=fdata.trace(DM%*%SN)^2/fdata.trace(MM%*%SN%*%MM%*%SN)
      f0=fdata.trace(DM%*%SN)^2/fdata.trace(DM%*%DM%*%SN%*%SN%*%Delta)
      pF=1-pf(FNM,f1,f0)
      ans[i,]=c(FNM,f1,f0,pF)
  }
  if (!is.null(contrast)){
jj2=ind.f=1
for (i in 1:b)    {
#ind.g=which(colnames(mf)==names(ncontrast[(i+1)]))#ok
ind.g=NA
ind.g=which(colnames(mfo)==names(contrast[i])) #ko
 if (is.vector(contrast[[i]])) vvf=matrix(contrast[[i]],ncol=1)
 else { if (is.list(contrast))    vvf=contrast[[i]]
        else  vvf=contrast}
for (jj in 1:ncol(vvf)) {
                   aa=mfo[,ind.g]
                   n.f=table(aa)
                   m.f=tapply(mfo[,1],mfo[,ind.g],mean)
                   sig.f=tapply(mfo[,1],mfo[,ind.g],var)
                   sig.f=sig.f*(n.f-1)/n.f^2
  max.vecm.f=length(m.f)
  vecm.f=matrix(as.vector(m.f),nrow=1)
  ssig.f=as.vector(sig.f)
#  MM.f=vvf%*%ginv(t(vvf)%*%vvf)%*%t(vvf)
  MM.f=vvf[,jj]%*%ginv(t(vvf[,jj])%*%vvf[,jj])%*%t(vvf[,jj])
  DM.f=diag(diag(MM.f))
  vecDelta.f=as.vector(1/(n.f-1))
  Delta.f=diag(vecDelta.f)
  SN.f=sum(n.f)*diag(ssig.f)
      FNM.f=sum(n.f)*vecm.f%*%MM.f%*%t(vecm.f)/fdata.trace(DM.f%*%SN.f)
      f1.f=fdata.trace(DM.f%*%SN.f)^2/fdata.trace(MM.f%*%SN.f%*%MM.f%*%SN.f)
      f0.f=fdata.trace(DM.f%*%SN.f)^2/fdata.trace(DM.f%*%DM.f%*%SN.f%*%SN.f%*%Delta.f)
      pF.f=1-pf(FNM.f,f1.f,f0.f)
      ans[length(nombres)+jj2,]=c(FNM.f,f1.f,f0.f,pF.f)
      jj2=jj2+1
}}}
res=list("ans"=ans)
if (!is.null(contrast)) res$contrast=contrast
if (pr) print(res)
class(res) <- "fanova.hetero"
return(res)
}

intercambio.l=function(contrast){
   n=length(contrast)
   mdata=contrast
   avance=seq(1,n)
   retro=seq(n,1)
   mdata[avance]=contrast[retro]
   names(mdata)[avance]= names(contrast)[retro]
return((mdata))
}

intercambio=function(fdata2){
   n=ncol(fdata2)
   mdata=fdata2
   avance=seq(2,n)
   retro=seq(n,2)
   mdata[,avance]=fdata2[,retro]
  colnames(mdata)[avance]= colnames(fdata2)[retro]
return(data.frame(mdata))
}


