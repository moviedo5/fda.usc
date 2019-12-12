#' DD-Classifier Based on DD-plot
#' 
#' @description Fits Nonparametric Classification Procedure Based on DD--plot
#' (depth-versus-depth plot) for G dimensions (\eqn{G=g\times h}{G=g x p}, g
#' levels and p data depth).
#' 
#' Make the group classification of a training dataset using DD-classifier
#' estimation in the following steps.\cr
#' 
#' \enumerate{ 
#' \item The function computes the selected \code{depth} measure of
#' the points in \code{fdataobj} w.r.t. a subsample of each g level group and p
#' data dimension (\eqn{G=g \times p}{G=g x p}).  The user can be specify the
#' parameters for depth function in \code{par.depth}.
#' 
#' (i) Type of depth function from functional data, see \code{\link{Depth}}:
#' \itemize{ 
#' \item \code{"FM"}: Fraiman and Muniz depth.  
#' \item \code{"mode"}: h--modal depth. 
#' \item \code{"RT"}: random Tukey depth. 
#' \item \code{"RP"}: random project depth.
#' \item \code{"RPD"}: double random project depth.
#'   }
#' (ii) Type of depth function from multivariate functional data, see \code{\link{depth.mfdata}}:
#' \itemize{ 
#' \item \code{"FMp"}: Fraiman and Muniz depth with common support.
#' Suppose that all p--fdata objects have the same support (same rangevals),
#' see \code{\link{depth.FMp}}.  
#' \item \code{"modep"}: h--modal depth using a p--dimensional metric, see \code{\link{depth.modep}}.
#' \item \code{"RPp"}: random project depth using a p--variate depth with the
#' projections, see \code{\link{depth.RPp}}.  
#' }
#' 
#' If the procedure requires to compute a distance such as in \code{"knn"} or \code{"np"} classifier or
#' \code{"mode"} depth, the user must use a proper distance function:
#' \code{\link{metric.lp}} for functional data and \code{\link{metric.dist}}
#' for multivariate data.
#' 
#' (iii) Type of depth function from multivariate data, see
#' \code{\link{Depth.Multivariate}}: 
#' \itemize{ 
#' \item \code{"SD"}: Simplicial depth (for bivariate data).  
#' \item \code{"HS"}: Half-space depth.  
#' \item \code{"MhD"}: Mahalanobis depth.
#' \item \code{"RD"}: random projections depth. 
#' \item \code{"LD"}: Likelihood depth.  
#' }
#' 
#' \item The function calculates the misclassification rate based on data depth
#' computed in step (1) using the following classifiers.
#' 
#' \itemize{ 
#' \item \code{"MaxD"}: Maximum depth.  
#' \item \code{"DD1"}: Search the best separating polynomial of degree 1.  
#' \item \code{"DD2"}: Search the best separating polynomial of degree 2.
#' \item \code{"DD3"}: Search the best separating polynomial of degree 3.
#' \item \code{"glm"}: Logistic regression is computed using Generalized Linear Models
#'  \code{\link{classif.glm}}.  
#' \item \code{"gam"}: Logistic regression is computed using Generalized Additive Models
#'  \code{\link{classif.gsam}}.
#' \item \code{"lda"}: Linear Discriminant Analysis is computed using
#' \code{\link{lda}}. 
#' \item \code{"qda"}: Quadratic Discriminant Analysis is computed using \code{\link{qda}}.  
#' \item \code{"knn"}: k-Nearest Neighbour classification is computed using \code{\link{classif.knn}}.  
#' \item \code{"np"}: Non-parametric Kernel classifier is computed using
#' \code{\link{classif.np}}.  
#' }
#' The user can be specify the parameters for classifier function in \code{par.classif} such as the smoothing parameter
#' \code{par.classif[[``h'']]}, if \code{classif="np"} or the k-Nearest
#' Neighbour \code{par.classif[[``knn'']]}, if \code{classif="knn"}.
#' 
#' In the case of polynomial classifier (\code{"DD1"}, \code{"DD2"} and
#' \code{"DD3"}) uses the original procedure proposed by Li et al. (2012), by
#' defalut rotating the DD-plot (to exchange abscise and ordinate) using in
#' \code{par.classif} argument \code{rotate=TRUE}. Notice that the maximum
#' depth classifier can be considered as a particular case of DD1, fixing the
#' slope with a value of 1 (\code{par.classif=list(pol=1)}).
#' 
#' The number of possible different polynomials depends on the sample size
#' \code{n} and increases polynomially with order \eqn{k}. In the case of
#' \eqn{g} groups, so the procedure applies some multiple-start optimization
#' scheme to save time:
#' 
#' \itemize{
#' 
#' \item generate all combinations of the elements of n taken k at a time:
#' \eqn{g \times combn(N,k)}{g x combs(N, k)} candidate solutions, and, when
#' this number is larger than \code{nmax=10000}, a random sample of
#' \code{10000} combinations.
#' 
#' \item smooth the empirical loss with the logistic function
#' \eqn{1/(1+e^{-tx})}{1/(1+e^{- tt x})}. The classification rule is
#' constructed optimizing the best \code{noptim} combinations in this random
#' sample (by default \code{noptim=1} and \code{tt=50/range(depth values)}).
#' Note that Li et al.  found that the optimization results become stable for
#' \eqn{t \in [50, 200]}{t between [50, 200]} when the depth is standardized
#' with upper bound 1.  } The original procedure (Li et al. (2012)) not need to
#' try many initial polynomials (\code{nmax=1000}) and that the procedure
#' optimize the best (\code{noptim=1}), but we recommended to repeat the last
#' step for different solutions, as for example \code{nmax=250} and
#' \code{noptim=25}. User can change the parameters \code{pol}, \code{rotate},
#' \code{nmax}, \code{noptim} and \code{tt} in the argument \code{par.classif}.
#' 
#' The \code{classif.DD} procedure extends to multi-class problems by
#' incorporating the method of \emph{majority voting} in the case of polynomial
#' classifier and the method \emph{One vs the Rest} in the logistic case
#' (\code{"glm"} and \code{"gam"}).
#' }
#' 
#' @param group Factor of length \emph{n} with \emph{g} levels.
#' @param fdataobj \code{\link{data.frame}}, \code{\link{fdata}} or \code{list}
#' with the multivariate, functional or both covariates respectively.
#' @param depth Character vector specifying the type of depth functions to use,
#' see \code{Details}.
#' @param classif Character vector specifying the type of classifier method to
#' use, see \code{Details}.
#' @param w Optional case weights, weights for each value of \code{depth}
#' argument, see \code{Details}.
#' @param par.depth List of parameters for \code{depth} function.
#' @param par.classif List of parameters for \code{classif} procedure.
#' @param control List of parameters for controlling the process.
#' 
#' If \code{verbose=TRUE}, report extra information on progress.
#' 
#' If \code{draw=TRUE} print DD-plot of two samples based on data depth.
#' 
#' \code{col}, the colors for points in DD--plot.
#' 
#' \code{alpha}, the alpha transparency used in the background of DD--plot, a
#' number in [0,1].
#' @return \itemize{
#' \item \code{group.est} {Estimated vector groups by classified method
#' selected.}  
#' \item \code{misclassification} { Probability of misclassification.} 
#' \item \code{prob.classification} { Probability of correct classification by group level.} 
#' \item  \code{dep} { Data frame with the depth of the curves for functional data (or points for multivariate data) in
#' \code{fdataobj} w.r.t. each \code{group} level.} 
#' \item \code{depth} { Character vector specifying the type of depth functions used.} 
#' \item \code{par.depth} { List of parameters for \code{depth} function.} 
#' \item \code{classif} { Type of classifier used.} 
#' \item \code{par.classif}{ List of parameters for \code{classif} procedure.}
#' \item \code{w}{ Optional case weights.} 
#' \item \code{fit}{ Fitted object by \code{classif} method using the depth as covariate.}
#' }
#' @author This version was created by Manuel Oviedo de la Fuente and Manuel
#' Febrero Bande and includes the original version for polynomial classifier
#' created by Jun Li, Juan A. Cuesta-Albertos and Regina Y. Liu.
#' 
#' @seealso See Also as \code{\link{predict.classif.DD}}
#' @references 
#' Cuesta-Albertos, J.A., Febrero-Bande, M. and Oviedo de la Fuente, M.
#' \emph{The DDG-classifier in the functional setting}, (2017). Test, 26(1),
#' 119-142. DOI: \url{https://doi.org/10.1007/s11749-016-0502-6}.
#' @keywords classif
#' @examples
#' \dontrun{
#' # DD-classif for functional data
#' data(tecator)
#' ab=tecator$absorp.fdata
#' ab1=fdata.deriv(ab,nderiv=1)
#' ab2=fdata.deriv(ab,nderiv=2)
#' gfat=factor(as.numeric(tecator$y$Fat>=15))
#' 
#' # DD-classif for p=1 functional  data set
#' out01=classif.DD(gfat,ab,depth="mode",classif="np")
#' out02=classif.DD(gfat,ab2,depth="mode",classif="np")
#' # DD-plot in gray scale
#' ctrl<-list(draw=T,col=gray(c(0,.5)),alpha=.2)
#' out02bis=classif.DD(gfat,ab2,depth="mode",classif="np",control=ctrl)
#' 
#' # 2 depth functions (same curves) 
#' out03=classif.DD(gfat,list(ab2,ab2),depth=c("RP","mode"),classif="np")
#' # DD-classif for p=2 functional data set
#' ldata<-list("ab"=ab,"ab2"=ab2)
#' # Weighted version 
#' out04=classif.DD(gfat,ldata,depth="mode",classif="np",w=c(0.5,0.5))
#' # Model version
#' out05=classif.DD(gfat,ldata,depth="mode",classif="np")
#' # Integrated version (for multivariate functional data)
#' out06=classif.DD(gfat,ldata,depth="modep",classif="np")
#' 
#' # DD-classif for multivariate data
#' data(iris)
#' group<-iris[,5]
#' x<-iris[,1:4]
#' out10=classif.DD(group,x,depth="LD",classif="lda")
#' summary(out10)
#' out11=classif.DD(group,list(x,x),depth=c("MhD","LD"),classif="lda")
#' summary(out11)
#' 
#' # DD-classif for functional data: g levels 
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-as.numeric(phoneme[["classlearn"]])-1
#' out20=classif.DD(glearn,mlearn,depth="FM",classif="glm")
#' out21=classif.DD(glearn,list(mlearn,mlearn),depth=c("FM","RP"),classif="glm")
#' summary(out20)
#' summary(out21)
#' }
#' 
#' @export
classif.DD <- function(group,fdataobj,depth="FM",classif="glm",w,
                       par.classif=list(),par.depth=list(),
                       control=list(verbose=FALSE,draw=TRUE,col=NULL,alpha=.25)){
  ################################################################################
  # classif.DD: Fits Nonparametric Classification Procedure Based on DD-plot
  # File created by Manuel Oviedo de la Fuente  using code from paper:
  # Cuesta-Albertos, J.A., Febrero-Bande, M. and Oviedo de la Fuente, M.
  # The DDG-classifier in the functional setting. 
  ################################################################################
  C<-match.call()
  if (!is.factor(group)) group<-factor(group)
  lev<-levels(group)
  group<-factor(group,levels=lev[which(table(group)>0)])
  if (is.null(control$verbose))  control$verbose<-FALSE
  if (is.null(control$draw))  control$draw<-TRUE
  if (is.null(control$fine))  control$fine<-100
  if (is.null(control$alpha))  control$alpha<-0.25
  if (is.null(control$prob))  control$prob<-0.5
  #  if (is.null(control$bg))  control$bg<-"white"
  # if (is.null(control$gray.scale))  control$gray.scale=FALSE
  classif0<-classif
  if (classif=="MaxD") {
    if (!is.null(par.classif$pol) & control$verbose) print("Maximum depth is done using polynomial argument, pol=1")
    par.classif$pol<-1
    classif<-"DD1" 
  }
  # if (classif=="knn" | classif=="np" | classif=="grm")   control$fine<-min(50,control$fine)
  func.clas<-list()
  draw<-control$draw
  
  par.Df<-list("df"=data.frame(group))
  # group<-ifelse(group==lev[1],0,1)
  ng<-length(lev)
  n<-length(group)#nrow(fdataobj)
  Df<-matrix(NA,ncol=ng,nrow=n)
  ind<-matrix(NA,nrow=n,ncol=ng)
  nvec<-table(group)
  p<-nvec[1]/n
  if ( is.fdata(fdataobj) | is.matrix(fdataobj) | is.data.frame(fdataobj)) {
    var.name<-deparse(substitute(fdataobj))
    ldata<-list(fdataobj)
    
    names(ldata)<-var.name
  }
  else {
    if (is.null(names(fdataobj))) names(fdataobj)<-paste("var",1:length(fdataobj),sep="")
    ldata<-fdataobj  
    var.name<-names(fdataobj)
  }
  lenlista<-length(ldata)
  lendepth<-length(depth)
  model<-FALSE
  par.ldata<-list()
  isfdata<-is.fdata(fdataobj) 
  if (is.null(names(par.depth))) {
    for (il in 1:lenlista)     par.ldata[[il]]<-par.depth     # esto no puede ser
  }
  else { 
    if (isfdata)  par.ldata[[1]]<-par.depth       
    else    par.ldata<-par.depth   }     
  lenl<-1
  ng2<-ng
  nam2<-c()
  nam<-NULL
  integrated<-FALSE
  if (missing(w)) {
    if (depth[1] %in% c("RPp","FMp","modep"))   {
      model=FALSE
      # integrated version modal
      multi<-TRUE
      depth0<-depth
      #   if (is.null(par.depth$scale)) par.depth$scale<-TRUE
      integrated<-TRUE  
      if (depth[1]=="modep"){
        hq<-numeric(ng)
        ismdist<-is.matrix(par.depth$metric)
        #     if (is.null(par.depth$par.metric$dscale)) par.depth$par.metric$dscale=
        if (ismdist)    mdist<-par.depth$metric
        else {#calculo  distancas para pasar el mismo h en todos
          par.depth$mfdata<-par.depth$mfdataref<-ldata
          oo<-do.call("depth.modep",par.depth)
          #           mdist<-oo$mdist
          par.depth$h<-oo$hq 
          #           ismdist<-TRUE
          hq<-rep(oo$hq,len=ng)
        }
        fdataobj<-ldata[[1]]
        par.depth$mfdata<-ldata
        for (i in 1:ng) {       
          ind[,i]<-group==lev[i]
          nam<-c(nam,paste("depth ",lev[i],sep=""))
          par.depth$mfdataref<-c.ldata(ldata,ind[,i]) 
          if (ismdist)  par.depth$metric<-mdist[,ind[,i],]
          oo<-do.call("depth.modep",par.depth)        
          par.depth$par.metric<-oo$par.metric
          Df[,i]<-oo$dep
          hq[i]<-oo$hq      
        }     
        w<- attributes(oo$mdist)$method
        par.depth$h<-hq
      } 
      if (depth[1] %in% c("FMp","RPp")){ 
        par.depth$mfdata<-ldata
        fdataori<-ldata
        nam.depth<-paste("depth.",depth[1],sep="")  
        for (i in 1:ng) {
          ind[,i]<-group==lev[i]
          fdataori<-c.ldata(ldata,ind[,i]) 
          nam<-c(nam,paste("depth ",lev[i],sep=""))
          par.depth$mfdataref<-fdataori 
          oo<-do.call(nam.depth,par.depth)
          Df[,i]<-oo$dep
          if (depth[1]=="RPp")  {par.depth$proj<-oo$proj     }
        }  
        w<-oo$dfunc
      }
      nam2<- paste(paste(var.name,collapse=".",sep=""),".",depth,".",lev,sep="")  
      gest<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth  
      depthl<-depth
      colnames(Df)<-nam2 
      par.ldata<-par.depth      
    }
    else{
      lenl<-lenlista
      w<-rep(1/lenlista,len=lenlista)
      ng2<-ng*lenlista
      model<-TRUE
    }   
  }
  if (!integrated){
    lenpesos<-length(w)
    if (sum(w)!=1) stop("Incorrect w argument, w must sum to 1")
    if (any(w<0))  stop("Incorrect w argument, w must be a positive")
    if (lenlista!=lenpesos) stop("Incorrect w argument")
    if (lendepth==1) depthl<-rep(depth,len=lenlista)
    else depthl<-depth
    depth0<-depth    
    for (idat in 1:lenlista) {
      fdataobj<-ldata[[idat]]
      depth<-depthl[idat]
      par.depth<-par.ldata[[idat]]
      nc<-ncol(fdataobj)
      x<-array(NA,dim=c(n,nc,ng))
      Df<-matrix(NA,ncol=ng,nrow=n)
      ind<-matrix(NA,nrow=n,ncol=ng)
      isfdata<-is.fdata(fdataobj)
      mnames<-ls(pattern="^mdepth.*",envir=as.environment("package:fda.usc"),all.names=TRUE)
      fnames<-ls(pattern="^depth.*",envir=as.environment("package:fda.usc"),all.names=TRUE) 
      mnames2<-ls(pattern="^mdepth.*",envir=.GlobalEnv,all.names=TRUE) 
      fnames2<-ls(pattern="^depth.*",envir=.GlobalEnv,all.names=TRUE)  
      mnames<-c(mnames,mnames2)
      fnames<-c(fnames,fnames2) 
      depth.long<-paste("mdepth.",depth,sep="")  
      if (depth.long %in% mnames & !isfdata) {
        multi<-TRUE
        par.depth$x<-fdataobj
      }
      else    {
        depth.long<-paste("depth.",depth,sep="")
        if (depth.long %in% fnames & isfdata) { 
          par.depth[["fdataobj"]]<-fdataobj
          multi=FALSE         }
        else stop("Incorrect depth function or data class object")
      }        
      if (depth %in% c("RHS","RP","RPD","RT")){
        if (is.null(par.depth$proj)) {
          d <- nc#-1
          u <- matrix(runif(d*25,-1,1),25,d)
          norm <- sqrt(rowSums(u*u))
          arg <- u/norm
          par.depth$proj<-arg                         #####################################hacer rproc2fdata y guardarlas para la prediccion!!!!!
          if (!multi & isfdata)  par.depth$proj<-rproc2fdata(50,fdataobj$argvals,sigma="vexponential")
          else   par.depth$proj<-arg
        }
      }
      ismdist<-is.matrix(par.depth$metric)
      if (ismdist) {   mdist<-par.depth$metric  }  
      dmode<-c(depth.long=="depth.mode" | depth.long=="mdepth.mode")
      if (dmode) {#calculo  distancas para pasar el mismo h en todos
        par.depth$fdataori<-par.depth$fdataobj
        oo<-do.call(depth.long,par.depth)
        mdist<-oo$mdist
        par.depth$h<-oo$hq 
        ismdist<-TRUE
        hq<-rep(oo$hq,len=ng)
      }
      #hq<-numeric(ng) 
      #  if (dmode) hq[i]<-oo$hq      
      for (i in 1:ng) {
        ind[,i]<-group==lev[i]
        nam<-c(nam,paste("depth ",lev[i],sep=""))#,paste("depth ",paste(lev[-i],collapse=",")))
        if (ismdist) {
          par.depth$metric<-mdist[,ind[,i]]   
          par.depth$metric2<-mdist[ind[,i],ind[,i]]
        }   
        if (multi)  par.depth$xx<-fdataobj[ind[,i],]
        else   par.depth$fdataori<-fdataobj[ind[,i],]
        oo<-do.call(depth.long,par.depth)
        if (depth %in% c("RHS","RP","RPD","RT")) par.depth$proj<-oo$proj
        Df[,i]<-oo$dep
        
      }
      if (dmode) par.depth$h<-hq
      for (i in 1:length(var.name)) var.name[i]<-unlist(strsplit(var.name[i], "[$]"))[[1]]
      for (i in 1:length(var.name)) var.name[i]<-unlist(strsplit(var.name[i], "[[]"))[[1]]     
      if (model) {
        if (idat==1) Df2<-Df
        else  Df2<-cbind(Df2,Df)
        par.ldata[[idat]]<-par.depth
        par.Df[[paste(".",idat,sep="")]]<-fdata(Df)
        nam2<-c(nam2, paste(var.name[idat],".",depthl[idat] ,".",lev,sep=""))
      }
      else{
        if (idat==1) {
          Df2<-w[idat]*Df
        }
        else  Df2<-Df2+w[idat]*Df
        nam2<- paste(nam2,paste(var.name[idat],depthl[idat],".",sep=""),sep="")
        par.ldata[[idat]]<-par.depth
      }
    }
    if (!model) nam2<- paste(nam2,lev,sep="")
    Df<-Df2
    gest<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth
    nvecs<-c(nvec,nvec)
    k<-1
    colnames(Df)<-nam2
  }
  if (draw) {
    if (is.null(control$col))     col1<-c(4,2,3,1,5:14)[1:ng]
    else col1<-control$col
    if (ng2>2) {
      if (control$verbose) {
        if (ng>2)  warning("DD-plot for more than 2 levels not implemented yet")
        else       warning("DD-plot for more than 1 depth function not implemented yet")
      }
      pch0<-c(21,24,22,23,21,24,22,23,21,24,22,23)[1:ng]   
      pch1<-pch0[group]   
    }
    else {
      rg<-range(Df)
      dff<-diff(rg)*.03
      # 2019/05/02
      #dff<-diff(rg)*.1
      minsq<-max(0,rg[1]-dff)
      maxsq<-rg[2]+dff 
      sq<-seq(minsq,maxsq,len=control$fine)
      vec<-data.frame(expand.grid(sq,sq))
      names(vec)<-nam2       
      pch0<-c(21,3) 
      pch1<-pch0[group]
      #     fill1<-c(4,2)
      #     mycols <- adjustcolor(palette("default"), alpha.f = control$alpha)
      #     opal <- palette(mycols)                
      fill1<- adjustcolor(col1, alpha.f = control$alpha)
    }                       
    col3<-col1
    col2<-col1[group]    
  }
  switch(classif,
         DD1={ 
           
           if (ng2>2) {
             #  stop("DD-plot for more than 2 levels not available")
             warning("Majority voting classification")
             par.classif$rotate<-FALSE #NOT IMPLEMENTED YET
             #ojo ng es num de grupos y ng2 ng*ndepth
             cvot<-combn(ng2,2)
             nvot<-ncol(cvot)
             votos<-matrix(0,ng,n)
             b0<-list()
             for (ivot in 1:nvot) {      
               #cat("votando",ivot)
               eps<-1e-10
               Df[Df==0]<-eps
               i2a2<-which(group==lev[cvot[1,ivot]] | group==lev[cvot[2,ivot]] )
               Df0<-Df[i2a2,cvot[,ivot]]
               ind0<-ind[i2a2,cvot[,ivot]]
               b<-unique(Df0[,1]/Df0[,2])
               mis <- sapply(b,MCR0,Df0,ind0)
               b0[[ivot]] <- min(b[which.min(mis)])
               group.log<-b0[[ivot]]*Df0[,1]<Df0[,2]
               votos[cvot[1,ivot],i2a2]<-votos[cvot[1,ivot],i2a2]+as.numeric(!group.log)
               votos[cvot[2,ivot],i2a2]<-votos[cvot[2,ivot],i2a2]+as.numeric(group.log)    
             }
             maj.voto<-apply(votos,2,which.max)
             group.est<-maj.voto
             for (ii in 1:n) {
               l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = T)]      
               if (length(l) > 1) {
                 abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = T))  
                 group.est[ii] =group[abc]
               }
               group.est <-  factor(group.est,levels = lev)
               incorrect<-group.est!=group
               mis<-mean(incorrect) 
               par.classif$pol<-b0
             }
             if (draw) { 
               draw=FALSE
               warning("Plot for majority voting classification not implemented yet")
             }
           }
           else {
             # en caso de empate dojo el grupo de los empatados con mayor profundiad: votos[cvot[1,ivot],i2a2]
             #En caso de empate ver cual ha sido la clasificacion entre esas 2 clases.
             if (is.null(par.classif$pol)) {
               b<-unique(Df[Df[,1]!=0&Df[,2]!=0,1]/Df[Df[,2]!=0&Df[,1]!=0,2])
               m <- length(b)
               mis <- rep(0,m)
               mis <- sapply(b,MCR0,Df,ind)
               b0 <- min(b[which.min(mis)])
               #####b1 <- optim(b0,AMCR0,dep=Df,ind=ind,tt=ii)$par
               
             } else {b0<-par.classif$pol}
             group.est<-factor(as.numeric(b0*Df[,1]<Df[,2]))
             levels(group.est)=lev
             
             incorrect<-group.est!=group
             mis<-mean(incorrect) 
             group.est2<-factor(as.numeric(b0*Df[,1]<Df[,2]),levels=lev)
             mis2<-mean(group.est!=group)
             par.classif$pol<-b0
             if (draw & ng2<=2) {
               rg<-range(Df)
               dff<-diff(rg)*.1
               bb2<-ifelse(b0*vec[,1]>vec[,2],0,1)
             }          
           }
         },
         DD2={
           #if (ng2>2) {stop("DD-plot for more than 2 levels not available")}
           if (is.null(par.classif$tt)) ind.tt<-50/min(max(Df),1)
           else ind.tt<-par.classif$tt#c(100,90, 20,80,140,70,160,60,180)
           if (is.null(par.classif$nmax))  { 
             nmax <- 100   
           }    else   {
             nmax<-par.classif$nmax
           }
           if (is.null(par.classif$noptim))  { 
             noptim <- 1
           }    else   {
             noptim<-par.classif$noptim
           }
           if (noptim>nmax) stop("nmax must be greather or equal to noptim")
           
           if (ng2>2) {
             par.classif$rotate<-FALSE #NOT IMPLEMENTED YET
             #  stop("DD-plot for more than 2 levels not available")
             warning("Majority voting classification")
             #ojo ng es num de grupos y ng2 ng*ndepth
             cvot<-combn(ng2,2)
             nvot<-ncol(cvot)
             votos<-matrix(0,ng,n)
             b0<-list()
             for (ivot in 1:nvot) {      
               #       cat("votando",ivot)
               eps<-1e-10
               Df[Df==0]<-eps
               i2a2<-which(group==lev[cvot[1,ivot]] | group==lev[cvot[2,ivot]] )
               Df0<-Df[i2a2,cvot[,ivot]]
               n0<-nrow(Df0)
               ind.0<-ind[i2a2,cvot[,ivot]]
               nsample<-choose(n0,2)
               combs1<-NULL
               if (nsample<nmax |  nmax==0 )  combs1 <- t(combn(n0,2))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n0,2))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }    
               #      b<-unique(Df0[,1]/Df0[,2])
               #      mis <- sapply(b,MCR0,Df0,ind0)
               #      b0 <- min(b[which.min(mis)])
               #      group.log<-b0*Df0[,1]<Df0[,2]
               if (noptim!=nmax) {
                 mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)    
                 ioptim<-order(mcrs)[1:noptim]
               } else {
                 ioptim<-1:nmax    
               }
               mcrs <- apply(combs1[ioptim,,drop=FALSE],1, quad.fit.opt, dep=Df0,ind=ind.0,tt=ind.tt)    
               mout <- matrix(unlist(mcrs),4)
               mcrs0<- mout[1,]
               wmin.mcrs<-which.min(mcrs0)
               tt.opt<-mout[2,wmin.mcrs]
               a0.2<- mcrs[[wmin.mcrs]][[3]]
               mcrs.opt<- mout[1,wmin.mcrs]
               group.log<-sapply(Df0[,1],RR,a=a0.2)<Df0[,2]
               #################
               votos[cvot[1,ivot],i2a2]<-votos[cvot[1,ivot],i2a2]+as.numeric(!group.log)
               votos[cvot[2,ivot],i2a2]<-votos[cvot[2,ivot],i2a2]+as.numeric(group.log)    
               b0[[ivot]]<-a0.2
             }
             par.classif$pol<-b0
             #    n<-nrow(Df)
             maj.voto<-apply(votos,2,which.max)
             group.est<-maj.voto
             for (ii in 1:n) {
               l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = T)]      
               if (length(l) > 1) {
                 abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = T))  
                 group.est[ii] =group[abc]
               }
               group.est <-  factor(group.est,levels = lev)
             }
             incorrect<-group.est!=group
             mis<-mean(incorrect)         
             if (draw) { 
               draw=FALSE
               warning("Plot for majority voting classification not implemented yet")
             }
           }
           else{ 
             # print("  #DD2 con 2 grupos")
             if (is.null(par.classif$pol)) {
               nsample<-choose(n,2)
               combs1<-NULL
               if (nsample<nmax |  nmax==0 )  combs1 <- t(combn(n,2))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n,2))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }
               Df<-Df+1e-15 ############################### sino NaN 
               #usar solo los q Df != 0 
               par.classif$mpol<-par.classif$mpol.opt<-NULL  
               par.classif$mis<-par.classif$mis.opt<-NULL  
               calculando<-TRUE
               if (is.null(par.classif$rotate))   par.classif$rotate<-TRUE
               iwhile=0
               isrotated<-FALSE
               tt.opt<-ind.tt[1]
               mis.new<-mis<-1
               a2<-c(1,1)
               while (iwhile<=2){
                 iwhile<-iwhile+1
                 if (noptim!=nmax) {
                   mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)    
                   ioptim<-order(mcrs)[1:noptim]
                 } else {
                   ioptim<-1:nmax    
                 }
                 mcrs <- apply(combs1[ioptim,,drop=FALSE],1, quad.fit.opt, dep=Df,ind=ind,tt=ind.tt)    
                 mout <- matrix(unlist(mcrs),4)
                 mcrs0<- mout[1,]
                 wmin.mcrs<-which.min(mcrs0)
                 tt.opt<-mout[2,wmin.mcrs]
                 a0.2<- mcrs[[wmin.mcrs]][[3]]
                 mcrs.opt<- mout[1,wmin.mcrs]
                 if (!any(is.na(a0.2))) {
                   input0<-sapply(Df[,1],RR,a=a0.2)
                   input1<-ifelse(input0 < Df[,2],lev[2],lev[1])            
                   group.est.new<-factor(input1,levels=lev) 
                   incorrect.new<-group.est.new!=group
                   mis.new<-mean(incorrect.new)
                   if (mis>mis.new){
                     if( par.classif$rotate & iwhile==3) isrotated<-TRUE
                     a2<-a0.2
                     mis<-mis.new
                     group.est<-group.est.new
                     #  par.classif$tt.opt<-lout$tt.opt      
                   }
                 }
                 if (par.classif$rotate & iwhile==1)   {
                   ####### #new exchange axis (rotate)
                   #      print("            se rotataaaaaaaaaaaaaaaaaaaaaaa  ")
                   Df<-Df[,2:1]
                   ind<-ind[,2:1]
                   iwhile=iwhile+1
                   lev<-lev[2:1]     
                 }    else {
                   iwhile=3
                 }
               }#fin calculando
               if (isrotated & par.classif$rotate){
                 warning("the axis are rotated")
                 #   lev<-lev[2:1]        
               }  else{
                 Df<-Df[,2:1]
                 ind<-ind[,2:1]  
                 par.classif$rotate<-FALSE #PQ EN LA PREDICCIONNO SE ROTA
               } 
               if (control$verbose) cat("pol=",a2," misclassification=",mis,"\n")
             }
             else a2<-par.classif$pol 
             # group.est<-factor(ifelse(sapply(Df[,1],RR,a=a2)<Df[,2],lev[2],lev[1])   )
             # incorrect<-group.est!=group
             # mis<-mean(incorrect) 
             par.classif$tt.opt <- tt.opt
             par.classif$pol<-a2
             if (draw & ng2<=2)  {
               bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],0,1)
               if (par.classif$rotate) bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],1,0)
             }
           }#fin else 2 grupos
         },
         DD3={
           if (is.null(par.classif$tt)) ind.tt<-50/min(max(Df),1)
           else ind.tt<-par.classif$tt#c(100,90, 20,80,140,70,160,60,180)
           # if (ng2>2) {stop("DD-plot for more than 2 levels not available")}
           if (is.null(par.classif$nmax))  { 
             nmax <- 100   
           }    else   {
             nmax<-par.classif$nmax
           }
           if (is.null(par.classif$noptim))  { 
             noptim <- 1
           }    else   {
             noptim<-par.classif$noptim
           }
           if (noptim>nmax) stop("nmax must be greather or equal to noptim")
           if (ng2>2) {
             #  stop("DD-plot for more than 2 levels not available")
             warning("Majority voting classification")
             #ojo ng es num de grupos y ng2 ng*ndepth
             cvot<-combn(ng2,2)
             nvot<-ncol(cvot)
             votos<-matrix(0,ng,n)
             b0<-list()
             for (ivot in 1:nvot) {      
               #cat("votando",ivot)
               eps<-1e-10
               #       Df[Df==0]<-eps
               i2a2<-which(group==lev[cvot[1,ivot]] | group==lev[cvot[2,ivot]] )
               Df0<-Df[i2a2,cvot[,ivot]]
               n0<-nrow(Df0)
               ind.0<-ind[i2a2,cvot[,ivot]]
               nsample<-choose(n0,3)
               combs1<-NULL
               if (nsample<nmax |  nmax==0)  combs1 <- t(combn(n0,3))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n0,3))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }   
               if (noptim!=nmax) {
                 mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)    
                 ioptim<-order(mcrs)[1:noptim]
               } else {
                 ioptim<-1:nmax    
               }
               mcrs <- apply(combs1[ioptim,,drop=FALSE],1, cubic.fit.opt, dep=Df0,ind=ind.0,tt=ind.tt)    
               mout <- matrix(unlist(mcrs),5)
               mcrs0<- mout[1,]
               wmin.mcrs<-which.min(mcrs0)
               tt.opt<-mout[2,wmin.mcrs]
               a0.3<- mcrs[[wmin.mcrs]][[3]]
               mcrs.opt<- mout[1,wmin.mcrs]
               #################
               group.log<-sapply(Df0[,1],RR,a=a0.3)<Df0[,2]
               votos[cvot[1,ivot],i2a2]<-votos[cvot[1,ivot],i2a2]+as.numeric(!group.log)
               votos[cvot[2,ivot],i2a2]<-votos[cvot[2,ivot],i2a2]+as.numeric(group.log)    
               b0[[ivot]]<-a0.3
             }
             par.classif$pol<-b0
             #    n<-nrow(Df)
             maj.voto<-apply(votos,2,which.max)
             group.est<-maj.voto
             for (ii in 1:n) {
               l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = T)]      
               if (length(l) > 1) {
                 abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = T))  
                 group.est[ii] =group[abc]
               }
               group.est <-  factor(group.est,levels = lev)
             }
             incorrect<-group.est!=group
             mis<-mean(incorrect)         
             if (draw) { 
               draw=FALSE
               warning("Plot for majority voting classification not implemented yet")
             }
           }
           else{
             # 3D con 2 grupos 
             # print("  #DD3 con 2 grupos")
             if (is.null(par.classif$pol)) {
               nsample<-choose(n,3)
               combs1<-NULL
               if (nsample<nmax |  nmax==0)  combs1 <- t(combn(n,3))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n,3))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }   
               mcrs <- apply(combs1,1,cubic.fit0, dep=Df,ind=ind)
               Df<-Df+1e-15 ############################### sino NaN 
               #usar solo los q Df != 0 
               par.classif$mpol<-par.classif$mpol.opt<-NULL  
               par.classif$mis<-par.classif$mis.opt<-NULL  
               calculando<-TRUE
               if (is.null(par.classif$rotate))   par.classif$rotate<-TRUE
               iwhile=0
               iwhile=0
               isrotated<-FALSE
               tt.opt<-ind.tt[1]
               mis.new<-mis<-1
               a2<-c(1,1)
               while (iwhile<=2){
                 iwhile<-iwhile+1
                 if (noptim!=nmax) {
                   mcrs <- apply(combs1,1, cubic.fit0, dep=Df,ind=ind)    
                   ioptim<-order(mcrs)[1:noptim]
                 } else {
                   ioptim<-1:nmax    
                 }
                 mcrs <- apply(combs1[ioptim,,drop=FALSE],1, cubic.fit.opt, dep=Df,ind=ind,tt=ind.tt)    
                 
                 mout <- matrix(unlist(mcrs),5)
                 mcrs0<- mout[1,]
                 wmin.mcrs<-which.min(mcrs0)
                 tt.opt<-mout[2,wmin.mcrs]
                 a0.2<- mcrs[[wmin.mcrs]][[3]]
                 mcrs.opt<- mout[1,wmin.mcrs]
                 if (!any(is.na(a0.2))) {
                   input0<-sapply(Df[,1],RR,a=a0.2)
                   input1<-ifelse(input0 < Df[,2],lev[2],lev[1])            
                   group.est.new<-factor(input1,levels=lev) 
                   incorrect.new<-group.est.new!=group
                   mis.new<-mean(incorrect.new)
                   if (mis>mis.new){      
                     if (par.classif$rotate & iwhile==3) isrotated<-TRUE
                     a2<-a0.2
                     mis<-mis.new
                     group.est<-group.est.new
                     #  par.classif$tt.opt<-lout$tt.opt      
                   }
                 }
                 if (par.classif$rotate & iwhile==1)   {
                   ####### #new exchange axis (rotate)
                   #      print("            se rotataaaaaaaaaaaaaaaaaaaaaaa  ")
                   Df<-Df[,2:1]
                   ind<-ind[,2:1]
                   iwhile=iwhile+1
                   lev<-lev[2:1]     
                 }    else {
                   iwhile=3
                 }
               }#fin calculando    
               if (isrotated & par.classif$rotate){
                 warning("the axis are rotated")
                 #   lev<-lev[2:1]        
               }  else{
                 Df<-Df[,2:1]
                 ind<-ind[,2:1]  
                 par.classif$rotate<-FALSE #PQ EN LA PREDICCIONNO SE ROTA
               } 
               if (control$verbose) cat("pol=",a2," misclassification=",mis,"\n")
             }
             else a2<-par.classif$pol 
             # group.est<-factor(ifelse(sapply(Df[,1],RR,a=a2)<Df[,2],lev[2],lev[1])   )
             # incorrect<-group.est!=group
             # mis<-mean(incorrect) 
             par.classif$tt.opt <- tt.opt
             par.classif$pol<-a2
             if (draw & ng2<=2)  {
               bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],0,1)
               if (par.classif$rotate) bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],1,0)
             }
           }},
         DDk={
           if (ng2>2) {
             par.classif$rotate<-FALSE #NOT IMPLEMENTED YET
             #  stop("DD-plot for more than 2 levels not available")
             warning("Majority voting classification")
             cvot<-combn(ng2,2)
             nvot<-ncol(cvot)
             votos<-matrix(0,ng,n)
             if (is.null(par.classif$nmax))  { 
               nmax <- 1000   
             }    else   {
               nmax<-par.classif$nmax
             }
             if (is.null(par.classif$noptim))  { 
               noptim <- 1
             }    else   {
               noptim<-par.classif$noptim
             }
             if (noptim>nmax) stop("nmax must be greather or equal to noptim")
             
             b0<-list()
             for (ivot in 1:nvot) {      
               #       cat("votando",ivot)
               eps<-1e-10
               Df[Df==0]<-eps
               i2a2<-which(group==lev[cvot[1,ivot]] | group==lev[cvot[2,ivot]] )
               Df0<-Df[i2a2,cvot[,ivot]]
               n0<-nrow(Df0)
               ind.0<-ind[i2a2,cvot[,ivot]]
               nsample<-choose(n0,2)
               combs1<-NULL
               if (nsample<nmax |  nmax==0 )  combs1 <- t(combn(n0,2))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n0,2))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }    
               #      b<-unique(Df0[,1]/Df0[,2])
               #      mis <- sapply(b,MCR0,Df0,ind0)
               #      b0 <- min(b[which.min(mis)])
               #      group.log<-b0*Df0[,1]<Df0[,2]
               mcrs <- apply(combs1,1, quad.fit0, dep=Df0,ind=ind.0)
               ind0 <- combs1[which.min(mcrs),] 
               A <- matrix(c( Df0[,1][ind0[1]],( Df0[,1][ind0[1]])^2, Df0[,1][ind0[2]],( Df0[,1][ind0[2]])^2),byrow=TRUE,2,2)
               ww <- c(Df0[,2][ind0[1]],Df0[,2][ind0[2]])
               a0.2 <- solve.ab(A,ww)
               group.log<-sapply(Df0[,1],RR,a=a0.2)<Df0[,2]
               votos[cvot[1,ivot],i2a2]<-votos[cvot[1,ivot],i2a2]+as.numeric(!group.log)
               votos[cvot[2,ivot],i2a2]<-votos[cvot[2,ivot],i2a2]+as.numeric(group.log)    
               b0[[ivot]]<-a0.2
             }
             par.classif$pol<-b0
             #    n<-nrow(Df)
             maj.voto<-apply(votos,2,which.max)
             group.est<-maj.voto
             for (ii in 1:n) {
               l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = T)]      
               if (length(l) > 1) {
                 abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = T))  
                 group.est[ii] =group[abc]
               }
               group.est <-  factor(group.est,levels = lev)
             }
             incorrect<-group.est!=group
             mis<-mean(incorrect)         
             if (draw) { 
               draw=FALSE
               warning("Plot for majority voting classification not implemented yet")
             }
           }
           else{      
             #DD2 con 2 grupos
             if (is.null(par.classif$pol)) {
               if (is.null(par.classif$nmax))   nmax=50000   #0
               else   nmax<-par.classif$nmax
               nsample<-choose(n,2)
               combs1<-NULL
               if (nsample<nmax |  nmax==0 )  combs1 <- t(combn(n,2))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n,2))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }
               mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)
               ind0 <- combs1[which.min(mcrs),] 
               A <- matrix(c( Df[,1][ind0[1]],( Df[,1][ind0[1]])^2, Df[,1][ind0[2]],( Df[,1][ind0[2]])^2),byrow=TRUE,2,2)
               ww <- c(Df[,2][ind0[1]],Df[,2][ind0[2]])
               a0.2 <- solve.ab(A,ww)
               #################################################################
               #      mcrs <- apply(comb1,1, quad.fit, Df=Df,Dg=Dg,x=x,y=y,n1=n1,n2=n2)
               #      ind <- comb1[which.min(mcrs),] 
               #      A <- matrix(c(Df[ind[1]],(Df[ind[1]])^2,Df[ind[2]],(Df[ind[2]])^2),byrow=TRUE,2,2)
               #      w <- c(Dg[ind[1]],Dg[ind[2]]) 
               #      a0.2 <- a2.MD <- solve(A,w)
               n0<-n
               nsample<-choose(n,3)
               combs1<-NULL
               if (nsample<nmax |  nmax==0)  combs1 <- t(combn(n0,3))
               else {
                 for(i in 1:nmax)      combs1<-rbind(combs1,sample(n0,3))
                 if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
               }   
               mcrs <- apply(combs1,1,cubic.fit0, dep=Df,ind=ind)
               ind0<- combs1[which.min(mcrs),]
               A <- matrix(c( Df[,1][ind0[1]],( Df[,1][ind0[1]])^2,( Df[,1][ind0[1]])^3,
                              Df[,1][ind0[2]],( Df[,1][ind0[2]])^2,( Df[,1][ind0[2]])^3,
                              Df[,1][ind0[3]],( Df[,1][ind0[3]])^2,( Df[,1][ind0[3]])^3),byrow=TRUE,3,3)
               ww <- c(Df[,2][ind0[1]],Df[,2][ind0[2]],Df[,2][ind0[3]])
               a0.3 <- solve.ab(A,ww)
             x.g1<-fdataobj[group==lev[1], ]
               x.g2<-fdataobj[group==lev[2], ]
               if (is.fdata(fdataobj))
                 nam.depth2<-paste("depth.",depth0,sep="")
               else
                 nam.depth2<-paste("mdepth.",depth0,sep="")
               mis.cv <- DD.cv.depth(x.g1,x.g2,a0.2,a0.3,nam.depth2)[-1]
               if (which.min(mis.cv)==1) {      
                 group.est<-factor(ifelse(Df[,1]<Df[,2],lev[2],lev[1])   )
               }      
               if (which.min(mis.cv)==2) {      
                 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a0.2)<Df[,2],lev[2],lev[1])   )
               }
               if (which.min(mis.cv)==3) {      
                 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a0.3)<Df[,2],lev[2],lev[1])   )
                 a0.2<-a0.3
               }
             
               incorrect<-group.est!=group
               mis<-mean(incorrect) 
               ####### #new exchange axis (rotate)
               if (is.null(par.classif$rotate))   par.classif$rotate<-FALSE #################
               
               if (par.classif$rotate)   {
                 Df<-Df[,2:1]
                 ind<-ind[,2:1]
                 mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)
                 ind0 <- combs1[which.min(mcrs),] 
                 A <- matrix(c( Df[,1][ind0[1]],( Df[,1][ind0[1]])^2, Df[,1][ind0[2]],( Df[,1][ind0[2]])^2),byrow=TRUE,2,2)
                 ww <- c(Df[,2][ind0[1]],Df[,2][ind0[2]])
                 a0.22 <- solve.ab(A,ww)
                 group.est2<-factor(ifelse(sapply(Df[,1],RR,a=a0.22)<Df[,2],lev[1],lev[2])   )#intercambio
                 incorrect<-group.est2!=group
                 mis2<-mean(incorrect) 
                 if (mis2<mis) {
                   a0.2<-a0.22
                   mis<-mis2
                   warning("the axis are rotated")
                   par.classif$rotate<-TRUE
                   lev<-lev[2:1]
                   group.est<-group.est2
                 }
                 else{
                   Df<-Df[,2:1]
                   ind<-ind[,2:1]  
                   par.classif$rotate<-FALSE
                 }  
               }
               a2<-a0.2
               if (control$verbose) cat("pol=",a2," misclassification=",mis,"\n")
               if (is.null(par.classif$tt)) ind.tt<-50/min(max(Df),1)
               else ind.tt<-par.classif$tt#c(100,90, 20,80,140,70,160,60,180)
               for (ii in ind.tt) {
                 a02 <- optim(a0.2,AMCR0,dep=Df,ind=ind,tt=ii)$par
                 group.est.new<-factor(ifelse(sapply(Df[,1],RR,a=a02)<Df[,2],lev[2],lev[1])   ) 
                 incorrect.new<-group.est.new!=group
                 mis.new<-mean(incorrect.new)
                 if (control$verbose) cat("pol=",a02," tt=",ii," misclassif=",mis.new,"\n")
                 if (mis>=mis.new) {a2<-a02;mis<-mis.new;group.est<-group.est.new}
               }
             }
             else a2<-par.classif$pol 
             group.est<-factor(ifelse(sapply(Df[,1],RR,a=a2)<Df[,2],lev[2],lev[1])   )
             incorrect<-group.est!=group
             mis<-mean(incorrect) 
             
             par.classif$pol<-a2
             if (draw & ng2<=2)  {
               bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],0,1)
               if (par.classif$rotate) bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],1,0)
             }
           }},
         lda={
           dat<-data.frame(as.numeric(group)-1,Df)
           par.classif$x<-Df
           par.classif$grouping<-group
           func.clas<-do.call(classif,par.classif)
           group.est<-predict(func.clas)$class
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {  
             bb2<-as.integer(predict(func.clas,vec)$class)
           }   
         },
         qda={
           dat<-data.frame(as.numeric(group)-1,Df)
           par.classif$x<-Df
           par.classif$grouping<-group    
           func.clas<-do.call(classif,par.classif)
           group.est<-predict(func.clas)$class
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {  
             bb2<-as.integer(predict(func.clas,vec)$class)
           }
         },
         knn={
           dat<-data.frame(as.numeric(group)-1,Df)
           names(dat)<-c("group1",nam2)
           par.classif$fdataobj<-Df
           par.classif$group<-group
           if (is.null(par.classif$knn))  par.classif$knn<-seq(3,19,by=2)
           func.clas<-do.call("classif.knn",par.classif)
           group.est<-func.clas$group.est
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {
             bb2<-as.numeric(predict(func.clas,vec))}
         },
         np={
           dat<-data.frame(as.numeric(group)-1,Df)
           names(dat)<-c("group1",nam2)
           par.classif$fdataobj<-Df
           par.classif$group<-group
           func.clas<-do.call("classif.kernel",par.classif)
           group.est<-func.clas$group.est
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {
             bb2<-as.numeric(predict.classif(func.clas,vec))
             }
         },
         rpart={
           dat<-data.frame(group,Df)
           names(dat)<-c("group1",nam2)
           par.classif$formula<-formula(paste("group1~",names(dat)[-1]))
           par.classif$data<-dat
           func.clas<-do.call("rpart",par.classif)
           group.est<-predict(func.clas,dat,type="class")
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {
             bb2<-predict(func.clas,vec,type="class")
             bb2<-ifelse(bb2==lev[1],1,2)     
           }
           
         },
         glm={
           dat<-data.frame(Df)
           names(dat)<-nam2
           par.classif$fdataobj<-dat
           par.classif$group<-group
           func.clas<-do.call("classif.glm2boost",par.classif)
#print(length(func.clas$fit))           
           group.est<-factor(func.clas$group.est,levels=lev)
           mis<-mean(group.est!=group)
           if (draw & ng2<=2) {
             bb2<-predict.classif(func.clas,list("df"=vec),type="class")
             bb2<-ifelse(bb2==lev[1],1,2)
           }
         },
         gam={
           dat<-data.frame(Df)
           names(dat)<-nam2  
           nam.classif<-paste("classif.",classif,sep="")
           par.classif$fdataobj<-dat
           par.classif$group<-group                    
           #if (is.null(par.classif$family)) par.classif$family<-binomial()
           
           func.clas<-do.call("classif.gsam2boost",par.classif)
           group.est<-func.clas$group.est
           incorrect<-group.est!=group
           mis<-mean(incorrect)
           if (draw & ng2<=2) {
             bb2<-predict.classif(func.clas,list("df"=vec),type="class")
             bb2<-ifelse(bb2==lev[1],1,2)
           }
         }
  )
  if (draw) {
    incorrect<-group.est!=group
    if (is.null(control$main)) {
      if (model) tit<-paste("DD-plot(",paste(depth0,collapse="-"),",",classif0,")",sep="")
      else    tit<-paste("DD-plot(",paste(depth0,w,collapse="-",sep=""),",",classif0,")",sep="")
    }
    else tit<-control$main   
    incorrect<-group!=group.est
    if (ng2>2) {
      bg1<-col1[group]
      bg1[!incorrect] <-0
      #       if (!control$gray.scale)               col2<-col1[group]       
      col2<-col1[group]       
      pairs(Df,col=col2,bg=bg1,main=tit,pch=pch1,diag.panel= function(...) {  
        legend("bottomleft",box.col=0,legend=lev,pch=pch0,col=col1,cex=max(.7,min(1,3/ncol(Df))),pt.bg=0,title="Points",title.adj=0.1,xjust=.1,yjust=.1)
        legend("bottomright",box.col=0,legend=c("good","bad"),pch=c(1,16),col=c(1,1),cex=min(1,3/ncol(Df)),pt.bg=c(0,1),title="Classify",title.adj=0.1,xjust=.1,yjust=.1)      
      } )} 
    else {
      #     mycols <- adjustcolor(palette("default"), alpha.f = control$alpha)
      #     opal <- palette(mycols)     
      #       if (control$gray.scale) {
      #          col2<-rep(1,len=15)[1:ng]
      #          col1<-gray(c(.75,.5),alpha=control$alpha)[1:ng]
      #          col3<-col2 
      #          fill1<-col1
      #      } 
      image(sq,sq,matrix(bb2,control$fine),xlab=nam[1],ylab=nam[2],main=tit,col=fill1,ylim=c(minsq,maxsq),useRaster=TRUE)
      #    palette("default")
      points(Df,col=col2,lwd=1,bg=col2,pch=pch1)
      #    palette(mycols) 
      sinpuntos<-TRUE
      zona<-c("topright","bottomright","topleft","bottomleft")
      izona<-1                                                                                               
      if (is.null(control$box.col)) control$box.col="black"    #borde legenda
      if (is.null(control$bg)) control$bg="white"
      while (sinpuntos) {
        le<- legend(zona[izona],title="Data points",legend=lev,col=col3,pt.bg=col3,border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)            #
        le2<-legend(zona[izona],title="Class Rule ",fill=fill1,legend=lev,border=1,pt.cex=2,pt.bg=col2,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)        
        sinpuntos<-switch(izona,
                          "1"={
                            xleg<-le$rect$left
                            xleg2<-le2$rect$left             
                            yleg2<-le$rect$top-le$rect$h -le2$rect$h
                            yleg<-le$rect$top-le$rect$h
                            sinpuntos= any(Df[,1]>xleg & Df[,2]>yleg2)},
                          "2"={
                            xleg<-le$rect$left
                            xleg2<-le2$rect$left             
                            yleg2<-le$rect$top-le2$rect$h  
                            yleg<-le$rect$top                      
                            sinpuntos= any(Df[,1]>xleg & Df[,2]<yleg)},
                          "3"={
                            xleg<-le$rect$left
                            xleg2<-le2$rect$left             
                            yleg2<-le$rect$top-le$rect$h -le2$rect$h
                            yleg<-le$rect$top-le$rect$h
                            sinpuntos= any(Df[,1]<(xleg+le$rect$w) & Df[,2]>yleg2) },
                          "4"={ 
                            xleg<-le$rect$left
                            xleg2<-le$rect$left             
                            yleg2<-le$rect$top-le2$rect$h
                            yleg<-le$rect$top
                            sinpuntos= any(Df[,1]<(xleg+le$rect$w) & Df[,2]<yleg2)})              
        
        if (!sinpuntos | izona==4) {
          if (!sinpuntos & izona==4 ) {
            control$box.col=0   #borde legenda
            control$bg=0
            le<- legend(zona[1],title="Data points",legend=lev,col=col3,pt.bg=col3,border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)            #
            le2<-legend(zona[1],title="Class Rule ",fill=fill1,legend=lev,border=1,pt.cex=2,pt.bg=col1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)         
            xleg<-le$rect$left
            xleg2<-le2$rect$left             
            yleg2<-le$rect$top-le$rect$h -le2$rect$h
            yleg<-le$rect$top-le$rect$h                     
          }     
          sinpuntos<-FALSE
        }
        else   izona<-izona+1
      }        
      legend(xleg2,yleg2, title="Class Rule ",fill=fill1,legend=lev,border=1, pt.bg=col2,
             horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0,yjust=0,box.col=control$box.col,bg=control$bg)             
      #     palette("default")
      legend(xleg,yleg,title="Data points",legend=lev,pch=pch0,col=col3,pt.bg=col3,
             border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0,yjust=0,box.col=control$box.col,bg=control$bg)
    }
  }
  prob.classification <- diag(table(group,group.est))/table(group)
  output <- list("group.est"=group.est,"misclassification"=mis,
               "prob.classification"=prob.classification,"dep"=Df,"depth"=depthl,
               "par.depth"=par.ldata,"classif"=classif,"par.classif"=par.classif,"w"=w,
               "group"=group,"fdataobj"=fdataobj,C=C,
               "model"=model,multi=multi,ldata=ldata
               ,prob=control$prob)

  output$fit <- func.clas
  class(output) <- "classif"
  return(output)
}