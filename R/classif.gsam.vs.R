#' @title Variable Selection in Functional Data Classification
#' 
#' @description Computes classification by selecting the functional (and non functional)
#' explanatory variables.  
#' @param data List that containing the variables in the model.  "df" element
#' is a  \code{data.frame} with the response and scalar covariates (numeric and factors
#' variables are allowed). Functional covariates of class \code{fdata} or
#' \code{fd} are introduced in the following items in the \code{data} list.
#' @param y  \code{caracter} string with the name of the scalar response variable
#' @param x  \code{caracter} string vector with the name of the scalar and functional
#' potential covariates.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param weights Weights:   
#' \itemize{
#' \item if \code{character} string \code{='equal'} same weights for each observation (by default) and
#' \code{='inverse'} for inverse-probability of weighting.   
#' \item if \code{numeric} vector of length \code{n}, Weight values of each observation.
#' }
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for functional beta parameter estimation.
#' @param prob probability value used for binary discriminant.
#' @param alpha alpha value to test the null hypothesis for the test of
#' independence among covariate X and residual e. By default is \code{0.05}.
#' @param dcor.min lower threshold for the variable X to be considered. X is
#' discarded if the distance correlation \eqn{R(X,e)< dcor.min} (e is the
#' residual).
#' @param type \code{character}, type of scheme classification. \code{'1vsall'}  (by default) 
#' strategy involves training a single classifier per class, with the samples of that class 
#' as positive samples and all other samples as negatives. Other posibility for K-way multiclass problem
#' is the \code{'majority'} voting scheme (also called one vs one). 
#' The procedure  trains the \eqn{K (K - 1) / 2}{K (K - 1) / 2}  binary classifiers and predicts the final class label as the class
#' label that has been predicted most frequently.
#' @param smooth if \code{TRUE}, a smooth estimate is made for all covariates included
#' in the model (less for factors). The model is adjusted with the estimated
#' variable linearly or smoothly. If the models are equivalent, the model is
#' adjusted with the linearly estimated variable.
#' @param measure measure related with correct classification (by default accuracy).
#' @param xydist list with the matrices of distances of each variable (all
#' potential covariates and the response) with itself.
#' @param \dots Further arguments passed to or from other methods.

# @param trace Interactive Tracing and Debugging of Call.
# @param verbose \code{logical}, if \code{TRUE} print relevant information for each iteration.
# @param CV TRUE, Cross-validation (CV) is done.
# @param ncomp.fix if TRUE, a number of basis element is fixed in
# \code{ncomp}. If FALSE, the funcion selects the number of PC components
# between the \code{ncomp}.
# @param par.model Model parameters.
# @param kbs the dimension of the basis used to represent the smooth term. The
# default depends on the number of variables that the smooth is a function of.

#' @aliases classif.gsam.vs
#' @return Return the final fitted model (same result of the classsification method) plus:\cr
#' \itemize{
#' \item \code{dcor}, \code{matrix} with the values of distance correlation for each
#' pontential covariate  (by column) and the residual of the model in each step (by row).
#' \item \code{i.predictor}, \code{vector} with 1 if the variable is selected, 0 otherwise.
#' \item \code{ipredictor}, \code{vector} with the name of selected variables (in order of selection)
#' }
#' 
#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' @seealso See Also as:  \code{\link{classif.gsam}}.
#' @references Febrero-Bande, M., Gonz\'alez-Manteiga, W. and Oviedo de la
#' Fuente, M. Variable selection in functional additive regression models,
#' (2018).  Computational Statistics, 1-19. DOI:
#' \doi{10.1007/s00180-018-0844-5}
#' @keywords classif
#' @note Adapted version from the original method in repression: \code{\link{fregre.gsam.vs}}.
#' @examples
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata
#' x1<-fdata.deriv(x)
#' x2<-fdata.deriv(x,nderiv=2)
#' y=factor(ifelse(tecator$y$Fat<12,0,1))
#' xcat0<-cut(rnorm(length(y)),4)
#' xcat1<-cut(tecator$y$Protein,4)
#' xcat2<-cut(tecator$y$Water,4)
#' ind <- 1:129
#' dat    <- data.frame("Fat"=y, x1$data, xcat1, xcat2)
#' ldat <- list("df"=dat[ind,],"x"=x[ind,],"x1"=x1[ind,],"x2"=x2[ind,])
#' # 3 functionals (x,x1,x2), 3 factors (xcat0, xcat1, xcat2)
#' # and 100 scalars (impact poitns of x1) 
#' 
#' res.gam<-classif.gsam(Fat~s(x2),data=ldat)
#' summary(res.gam1$model)
#' 
#' # Time consuming
#' res.gam.vs<-classif.gsam.vs("Fat",data=ldat)
#' summary(res.gam.vs)
#' res.gam.vs$i.predictor
#' res.gam.vs$ipredictor
#' 
#' # Prediction 
#' newldat <- list("df"=dat[-ind,],"x"=x[-ind,],
#'                 "x1"=x1[-ind,],"x2"=x2[-ind,])
#' pred.gam<-predict(res.gam,newldat)                
#' pred.gam.vs<-predict(res.gam1,newldat)
#' cat2meas(newldat$df$Fat,pred.gam)
#' cat2meas(newldat$df$Fat,pred.gam.vs)
#' }
#' @export 
classif.gsam.vs = function(data = list(), y, x, family = binomial()
                         , weights = "equal"
                         , basis.x = NULL, basis.b = NULL
                         ,type = "1vsall", prob = 0.5, alpha = 0.05
                         ,dcor.min =0.01, smooth=TRUE
                        ,measure = "accuracy", xydist,...){
  #data,y,x,
  if (missing(y)) {stop("The name of the response must be specified in the 'y' argument")  }
  resp <- y
  n=length(data$df[[resp]])
  nesc=length(names(data$df))-1 #Number of scalar variables
  namesc=names(data$df)[which(names(data$df)!=resp)]
  namfunc=names(data)[which(names(data)!="df")]
  
  if (missing(x)) {  
    x<-c(namesc,namfunc)
    ifunc<-namfunc
    infunc<-namesc
  }  else {
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(y,x),namesc)
  } 
  
  # xdatos<-as.list(data$df[,y,drop=F])
  # xdatos<-c(xdatos,as.list(data$df[,infunc,drop=F]),data[ifunc])
  
  #nvar=length(names(data$df))-1+length(data)-1 #Global number of variates
  nvar =length(x)
  
  dcor = matrix(0,nrow=nvar,ncol=nvar)
  ind = numeric(nvar)
  colnames(dcor) = x #c(namesc,namfunc)
  covar0 <- covar <- x #c(namesc,namfunc)
  names(ind) = colnames(dcor)
  rownames(dcor) = 1:nvar
  form = as.formula(paste0(resp,"~1"))
  
  j=1
  Mset=c()
# print(form)  
  modelo <-classif.gsam(formula=form, family = family, data=data,basis.x=basis.x,...)
  #res<-NULL
  #for (i in 1:length(modelo$fit))
#    res<-cbind(res,modelo$fit[[i]]$residuals)
  res0<-res<-model.matrix(as.formula(paste0("~-1+",resp)),data$df)
 #print("aa")
 #print(head(res0)  )
  #res=data$df[,resp]-mean(data$df[,resp])
  modfinish=FALSE
  len.esc <- length(namesc)
  len.fun <- length(namfunc)
  i.predictor<-NULL
  ipredictor<-numeric(nvar)
  names(ipredictor)<-x
#print("a1")
if (missing(xydist)) {
#  print("entra xydist")
  calc.dist<-TRUE
  xydist<-list()
  xydist[[resp]]<-metric.dist(res)  
  if (len.esc>0){
  #  cat(len.esc,"Escalar j=",j,"\n") 
    for (i in 1:len.esc){
      aux<-NA
      if (is.factor(data$df[[namesc[i]]])) {
        y=model.matrix(as.formula(paste0("~-1+",namesc[i])),data$df)
      } else {
        y=data$df[,namesc[i],drop=F]
      }
#      cat(i,namesc[i],dim(y),"\n")
      xydist[[namesc[i]]]<-metric.dist(y) 
    }
#    print(names(xydist))
    }
  if (len.fun>0){
#    cat("Funcional j=",j,"\n")        
    for (i in 1:len.fun){
      xydist[[namfunc[i]]]<-metric.lp(data[[namfunc[i]]])
    }}
  
}# else  calc.dist<-FALSE
#print("a2")
len.var<-length(covar)
#print(calc.dist)
iiii<-0#borrraar
while(!modfinish){
  iiii<-iiii+1 #borrraar
#  cat("iiiii ");print(iiii)
  iname<-NULL
  xydist[[resp]]<-metric.dist(res)
  len.var<-length(covar)
  for (i in 1:len.var){
      icovar<-covar[i]
      tt=dcor.test(xydist[[resp]],xydist[[icovar]],n=n)
      #dcor[j,i]=tt$estimate*(tt$p.value<alpha)
      dcor[j,icovar]=tt$estimate*(tt$p.value<alpha)
  }
  
#print(max(dcor[j,]))
#print(which.max(dcor[j,]))
 #dcor[j,covar %in% Mset]=-dcor[j,covar %in% Mset]
 jj=which.max(dcor[j,])
 aux<-jj
 iname<-covar0[jj]
#print("iname") ;print(iname)
 if  (iname %in% namfunc) esfuncional<-TRUE
 else esfuncional<-FALSE
#print(esfuncional)
#print(dcor[j,])
#print(Mset) ;print(jj)
#print(colnames(dcor)) 
 
# cat(paste0("Iter:",j," dcor(",covar0[jj],")= ",round(dcor[j,jj],4),"\n"))  
 if (max(dcor[j,],na.rm=TRUE) > dcor.min ) {
		modant=modelo
		formant=form
		Mset=c(Mset,iname)
		covar<-setdiff(covar,iname)
		if (!esfuncional){
 # print("escalar")		  
		  #namesc[jj]
  		esfactor<-is.factor(data$df[,iname])
  		namesc<-setdiff(namesc,iname)
  		len.esc<-len.esc-1
		} 
		else { 
# print("funcional")
		 esfactor<-FALSE
		 namfunc<-setdiff(namfunc,iname)
		 len.fun<-len.fun-1
		}
		a1<- Inf
# cat("esfactor " ,esfactor,"\n")
    if (esfactor | !smooth ){
#print("es factor o !smooth")
		form.lin=update.formula(form, paste0(".~.+",iname))
#old		  modelo=fregre.gsam(formula=form.lin,family = family, data=data,basis.x=basis.x,...)
		modelo=classif.gsam(formula=form.lin,family = family, data=data,basis.x=basis.x,...)
		#a1=do.call(msc,list("model"=modelo))	
		a1 <-  cat2meas(modelo$group,modelo$group.est, measure = measure)
    }
    if (!esfactor){
      # print("smooth entra")      
#print(form)
     	 form.smo = update.formula(form, paste0(".~.+s(",iname,")"))
		   modelo.smo=classif.gsam(formula=form.smo,family = family, data=data,basis.x=basis.x,...)
		   # print("smooth sale")		   
#print(modelo.smo)
#print(form.smo)
	     #a2=do.call(msc,list("model"=modelo.smo))
	     a2 <-  cat2meas(modelo.smo$group,modelo.smo$group.est, measure = measure)
# print("lineal vs smo22")
#cat(a1,a2,(a2/a1)<1,"\n")
		   if ((a2/a1)<1) {
#print("smoo22")
		     modelo <-modelo.smo
		     form <- form.smo
		   }
	     else form <- form.lin
		} else form <- form.lin
		#print(form)
	j=j+1

	if (j>nvar){modfinish=TRUE} # entran todas o ninguna
	} else {modfinish=TRUE}
  #!modfinish
   #
#	#print(names(modelo))
##########	 
# aa=do.call(msc,list("model"=modant))
# ab=do.call(msc,list("model"=modelo))

aa <-  cat2meas(modant$group,modant$group.est, measure = measure)
ab <-  cat2meas(modelo$group,modelo$group.est, measure = measure)

#print(msc)
#cat(aa,ab,(ab/aa)<=1,"\n")
#	#print(names(modelo))
 #print(" #################      modelo anterior")
# print(summary(modant))
# print(" #################      modelo actual")
# print(summary(modelo))
if ((ab/aa) <= 1){
  # print("ni entra ni actualiza")
	form=formant
	modelo=modant	
#	modfinish=TRUE
} else {
  if (!modfinish) {
    # print("Entra y actualiza")    
#print(iname)
    i.predictor<-c(i.predictor,iname)
    ipredictor[aux] <- 1
  }
#  res<-NULL
  #for (i in 1:length(modelo$fit))
   # res<-cbind(res,modelo$fit[[i]]$residuals)
# print(head(res0))  
# print(head(modelo$prob.group))
# si 2 grupos solo cogemos 1 grupo
  if (ncol(res0)==1)  res<-res0-modelo$prob.group[,1]
  else res<-res0-modelo$prob.group
  }
# print(head(res))
# print("aaa")
	 if (ab==1) {modfinish=TRUE} # prob.classif=1!
}
#return(list(form=form, data=data, basis.x=basis.x, model=modelo,
#            dcor=dcor,i.predictor=ipredictor,ipredictor=i.predictor))
modelo$dcor=dcor
modelo$formula.ini <-modelo$formula
modelo$i.predictor=ipredictor
modelo$ipredictor=i.predictor
return(modelo)
}
