# For multiclass-classification with k levels, k>2, libsvm uses the ‘one-against-one’-approach, in which k(k-1)/2 binary classifiers are trained; the appropriate class is found by a voting scheme.

wlda2<-function (x, grouping, weights = rep(1, nrow(x)), 
                 method = c("unbiased", "ML"), ...) 
{
  if (is.null(dim(x))) 
    stop("'x' is not a matrix")
  x <- as.matrix(x)
  if (any(!is.finite(x))) 
    stop("infinite, NA or NaN values in 'x'")
  n <- nrow(x)
  if (n != length(weights)) 
    stop("nrow(x) and length(weights) are different")
  if (any(weights < 0)) 
    stop("weights have to be larger or equal to zero")
  if (all(weights == 0)) 
    stop("all weights are zero")
  names(weights) <- rownames(x)
  if (n != length(grouping)) 
    stop("'nrow(x)' and 'length(grouping)' are different")
  x <- x[weights > 0, , drop = FALSE]
  grouping <- grouping[weights > 0]
  w <- weights[weights > 0]
  n <- nrow(x)
  if (!is.factor(grouping)) 
    warning("'grouping' was coerced to a factor")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if (any(counts == 0)) {
    empty <- lev[counts == 0]
    warning(sprintf(ngettext(length(empty), "group %s is empty or weights in this group are all zero", 
                             "groups %s are empty or weights in these groups are all zero"), 
                    paste(empty, collapse = ", ")), domain = NA)
    lev1 <- lev[counts > 0]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  if (length(lev1) == 1L) 
    stop("training data from only one group given")
  method <- match.arg(method)
  class.weights <- tapply(w, g, sum)
  prior <- class.weights/sum(w)
  print(w)
  print(prior)
  print(class.weights)
  print("a")
  ng <- length(prior)
  names(counts) <- lev1
  xwt <- w * x
  print(dim(x))
  print(dim(xwt))
  center <- t(matrix(sapply(lev1, function(z) colSums(xwt[g == 
                                                            z, , drop = FALSE])), ncol = ng, dimnames = list(colnames(x), 
                                                                                                             lev1)))/as.numeric(class.weights)
  z <- x - center[g, , drop = FALSE]
  cov <- crossprod(w * z, z)/sum(w)
  
  if (method == "unbiased") {
    norm.weights <- w/class.weights[g]
    cov <- cov/(1 - sum(prior * tapply(norm.weights^2, g, 
                                       sum)))
  }
  cl <- match.call()
  cl[[1L]] <- as.name("wlda")
  return(structure(list(prior = prior, counts = counts, means = center, 
                        cov = cov, lev = lev, N = n, weights = weights, method = method, 
                        call = cl), class = "wlda"))
}
# https://github.com/rforge/locclass/tree/master/locClassMlr
#install.packages("mlbench")
library(mlbench)
data(PimaIndiansDiabetes)

train <- sample(nrow(PimaIndiansDiabetes), 500)

# weighting observations from classes pos and neg according to their
# frequency in the data set:
ws <- as.numeric(1/table(PimaIndiansDiabetes$diabetes)
                 [PimaIndiansDiabetes$diabetes])
sum(ws[train])
fit <- locClass:::wlda.default(as.matrix(PimaIndiansDiabetes[train,-9]),
                               PimaIndiansDiabetes[train,"diabetes"],
            , weights = ws[train])

fit2 <- wlda2(as.matrix(PimaIndiansDiabetes[train,-9]),
              PimaIndiansDiabetes[train,"diabetes"],
              , weights = rev(ws[train]))
pred <- predict(fit, newdata = PimaIndiansDiabetes[-train,-9])
pred2 <- predict(fit2, newdata = PimaIndiansDiabetes[-train,-9])
ytest<-  PimaIndiansDiabetes[-train,9]
mean(pred$class==ytest)
mean(pred2$class==ytest)
pred$posterior[1:11,];pred2$posterior[1:11,]
pred$posterior[1:11,]-pred2$posterior[1:11,]


install.packages("flexmix")
install.packages("locClass", repos="http://R-Forge.R-project.org")
library(locClass)
?wlda
# Hand, D. J., Vinciotti, V. (2003), Local versus global models for classification problems: Fitting models where it matters, The American Statistician, 57(2) 124–130.
?dalda

#' @export classif.gbm
classif.gbm=function(formula, data, basis.x=NULL, 
                              weights = "equal", type = "1vsall", ...) 
{
  rqr <- "randomForest"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'randomForest'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights","type"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  #mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  y <- data$df[[response]]
  lev <- levels(y)
  prob2<-prob1 <- ny <- nlevels(y)
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x,pf,tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  vs.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- NROW(XX)
  
  par.method <- as.list(substitute(list(...)))[-1L]
  if (weights[1] =="equal") 
    class.weights=NULL
  
  if (weights[1] == "inverse")   {
    #print("inverse")
    weights<-weights4class(y,type=weights)
    wtab<-tapply(weights,y,mean)
    class.weights<-wtab/sum(wtab)
    names(class.weights)<- lev
  }
  
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  if (length(weights)==ny){
    #print("weights como clas lev")
    class.weights <- weights
    names(class.weights) <- lev
  } 
  # if (is.null(names(class.weights)))     names(class.weights) <- lev
  if (!is.null(par.method$classwt)) 
    class.weights=classwt
  
  par.method <- c(list(formula=pf, data=XX,classwt=class.weights),par.method)
  #par.method<-c(list(formula=pf, data=XX),par.method)
  if (is.null(par.method$votes)) par.method$votes=TRUE
  
  
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  
  
  if (type == "majority" |  ny==2){
    z=do.call("randomForest",par.method)
    out$fit<-z
    out$group.est = z$predicted
    out$prob.group <- z$votes
    
  }    else { # One vs Other
    # 2019/08/30
    prob.group<-matrix(NA,n,ny)
    colnames(prob.group)<-lev
    z<-list()
    for (i in 1:ny) {
      igroup  <- y==lev[i]
      newy<-ifelse(igroup, 0,1)
      par.method$data[,response]<-factor(newy)
      newx<- par.method$data
      z[[i]] <-  do.call("randomForest",par.method)
      prob.group[,i] <- z[[i]]$votes[,1]
    }
    out$prob.group <- prob.group
    out2glm<-classifKgroups(y,prob.group,lev) # hacer una par prob<0 y >0
    out$group.est = out2glm$yest
    out$fit <- z
  }
  out$prob <- prob
  out$group <- y
  out$max.prob <- mean(y==out$group.est,na.rm=T) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$vs.list=vs.list
  #out$method <- method
  #out$par.method <- par.method
  
  tab <- table(out$group.est,y)
  prob.group2 <- array(NA, dim = c(n, ny))
  prob.group2 <- prob.group2/apply(prob.group2, 1, sum)
  for (i in 1:ny) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  #colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  out$type <- type
  # If norm.votes=TRUE, the fraction is given, which can be taken as predicted probabilities for the classes.
  
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out)<-"classif"
  out
}

# install.packages("gbm")
library(gbm)

?gbm



#############################
ii<-c(1:150,191:210)
dat$df<-data.frame(glearn=phoneme$classlearn[ii])
dat$x<-phoneme$learn[ii,]
# plot(dat$x)
newdat<-list("df"=data.frame("glearn"=phoneme$classtest),"x"=phoneme$test)
tp<-"majority"
tp1<-"1vsall"
wt<-"inverse"
wt2<-"equal"
set.seed(1)
args(classif.multinom)
a5<-classif.multinom(glearn~x, data = dat, weights = wt2)

set.seed(1)
a6<-classif.nnet(glearn~x, data = dat, type=tp,weights = wt2)

a5$max.prob;a6$max.prob
a5$prob.classification;a6$prob.classification

names(a5)
a5$type
tp2<-"class"
p5 <- predict(a5,newdat,type=tp2) #no va probs!!
p5
p6<- predict(a6,newdat,type=tp2) #no va probs!!
a5$max.prob;a6$max.prob
a5$prob.classification;a6$prob.classification
mean(p5==newdat$df$glearn);mean(p6==newdat$df$glearn)


set.seed(1)
a6<-classif.rpart(glearn~x, data = dat, type=tp,weights="inverse")

set.seed(1)
a6<-classif.rpart(glearn~x, data = dat, type=tp,weights=c(1,1,1,1,1))

p5<- predict(a5,newdat);p6<- predict(a6,newdat)
a5$max.prob;a6$max.prob
a5$prob.classification;a6$prob.classification
mean(p5==newdat$df$glearn);mean(p6==newdat$df$glearn)


table(p5,a5$group.est)
mean(p5==dat$df$glearn)

a5<-classif.randomForest(glearn~x, data = dat, weights="equal",ntree=10)
set.seed(1)
p5<- predict(a5,dat)

p5<- predict(a5,dat,type="prob")
names(p5)
table(p5,a5$group.est)
a5$type
mean(p5==dat$df$glearn)

set.seed(1)
a5<-classif.randomForest(glearn~x, data = dat, weights="equal", type="majority")
set.seed(1)
p5<- predict.classif(a5,dat)
table(p5,a5$group.est)

####################################
a4<-classif.randomForest(glearn~x, data = dat, weights="equal")
a4$prob.group
a5<-classif.randomForest(glearn~x, data = dat, weights="equal", type="majority")
a5$prob.group
a4$max.prob;a5$max.prob

p5<- predict(a5$fit, newx = a3$XX[,-1])
table(p5,a5$group.est)


p4<- predict(a4,dat)
p5<- predict(a5,dat)
p6<- predict.classif(a5,dat)
a5$C;class(a5)
table(p5,a5$group.est)
table(p6,a5$group.est)
table(p6,p5)


a<- predict(a3$fit[[1]], newx = ,a3$fit[[1]]$x,type ='linear.predictor')
a<- predict(a3$fit[[1]], newx = a3$fit[[1]]$x,type ='votes')
a
####################################
fdata2model<-fda.usc:::fdata2model
weights4class<-fda.usc:::weights4class

 require(fda.usc)
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
mtest<-phoneme[["test"]]
gtest<-phoneme[["classtest"]]
dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.glm(glearn~x, data = dat)
a2<-classif.svm(glearn~x, data = dat)
#a3<-classif.ksvm(glearn~x, data = dat, weights="equal")
a4<-classif.lda(glearn~x, data = dat)
a5<-classif.lda(glearn~x, data = dat, weights=c(.1,.2,.2,.25,.25))
a4$max.prob
a5$max.prob
names(a4) 
names(a4$fit[[1]])
summary(a4)
a4$fit

class(a4$fit)

# 2019/08/29
# cbind(round(a2$fit$decision.values,3),a2$fit$fitted)[200:250,]
# 1/2    1/3    1/4    1/5    2/3    2/4    2/5    3/4    3/5    4/5  
# svm hace majority voting por defecto
# se incluye "1vsall" para que sea homogeneo con el resto
#################################
classifKgroups<-fda.usc:::classifKgroups

type="inverse" 
type="equal"
a1<-classif.glm(glearn~x, data = dat,weights=type)
a2<-classif.glm(glearn~x, data = dat,weights=type,type="majority")
a3<-classif.svm(glearn~x, data = dat,weights=type)
a4<-classif.svm(glearn~x, data = dat,weights=type,type="majority")

a1<-classif.lda(glearn~x, data = dat)
a2<-classif.lda(glearn~x, data = dat, weights="inverse")
a3<-classif.lda(glearn~x, data = dat, type="majority")

a1<-classif.qda(glearn~x, data = dat)
a2<-classif.qda(glearn~x, data = dat, weights="inverse")
a3<-classif.qda(glearn~x, data = dat, type="majority")
a4<-classif.lda(glearn~x, data = dat)


a4$group.est
a3$prob.classification 
a4$prob.classification 
a3$prob.group[1:3,]
a4$prob.group[1:3,]
a1$max.prob;a2$max.prob;a3$max.prob;a4$max.prob
newdat<-list("x"=phoneme$test)
ytest<-phoneme$classtest
ytrain<-dat$df$glearn
newdat<-dat;ytest<-ytrain

p1<-predict.classif(a1,newdat)
p2<-predict.classif(a2,newdat)
p3<-predict.classif(a3,newdat)
p4<-predict.classif(a4,newdat)

mean(p1==ytest);mean(p2==ytest);mean(p3==ytest);mean(p4==ytest);
mean(p1==ytrain);mean(p2==ytrain);mean(p3==ytrain);


newdat<-list(df=data.frame(glearn=phoneme$classlearn),"x"=phoneme$test)
ytest<-phoneme$classtest
ytrain<-phoneme$classtest

a1<-classif.qda(glearn~x, data = newdat)
a2<-classif.qda(glearn~x, data = newdat, weights=c(.3,.3,.2,.1,.1))

p1<-predict(a1$fit[[1]],a1$XX)$class
p2<-predict(a2$fit[[1]],a2$XX)$class
p3<-predict.classif(a1,newdat)
mean(p1==p2)
mean(p1==ytest);mean(p2==ytest)
mean(p1==ytrain);mean(p2==p3)
table(p1,p2)
table(p2,p3)
#################################
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])
## your answer may differ
##  c  s  v
## 22 23 30
z <- lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
p<-predict(z, Iris[train, ])$class

names(Iris)
Iris2<-Iris[1:100,]
Iris2$Sp2<-as.numeric(Iris2$Sp)-1
z2 <- gbm(formula("Sp2 ~ Sepal.L.+Sepal.W.+Petal.L.+Petal.W."), data=Iris2,distribution="bernoulli")
z2
summary(z2)

predict.gbm(z2, Iris2,n.trees=10, single.tree=F )>0
predict.gbm(z2, Iris2,n.trees=10, single.tree=F,type="response" )>0.5

install.packages("xgboost")
library(xgboost)
?xgboost


install.packages("h2o")
library(h2o)
#http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
#cuenta de usuario (como aemet!)
h2o.init()
#For H2O package documentation, ask for help:
?h2o

h2o.no_progress()
h2o.init(max_mem_size = "5g")

names(Iris2)
y <- "Sp2"
x <- setdiff(names(Iris2), y)

# turn training set into h2o object
train.h2o <- as.h2o(Iris2)

# training basic GBM model with defaults
h2o.fit1 <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train.h2o,
  nfolds = 5
)

# assess model results
(h2o.fit1)


# training basic GBM model with defaults
h2o.fit2 <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train.h2o,
  nfolds = 5,
  ntrees = 5000,
  stopping_rounds = 10,
  stopping_tolerance = 0,
  seed = 123
)

# model stopped after xx trees
h2o.fit2@parameters$ntrees

# cross validated RMSE
h2o.rmse(h2o.fit2, xval = TRUE)




#split <- h2o.splitFrame(train.h2o, ratios = 0.75)
#train <- split[[1]]
#valid <- split[[2]]

train <- train.h2o
valid <- train.h2o

# create hyperparameter grid
hyper_grid <- list(
  max_depth = c(1, 3, 5),
  min_rows = c(1, 5, 10),
  learn_rate = c(0.01, 0.05, 0.1),
  learn_rate_annealing = c(.99, 1),
  sample_rate = c(.5, .75, 1),
  col_sample_rate = c(.8, .9, 1)
)

# perform grid search 
grid <- h2o.grid(
  algorithm = "gbm",
  grid_id = "gbm_grid1",
  x = x, 
  y = y, 
  training_frame = train,
  validation_frame = valid,
  hyper_params = hyper_grid,
  ntrees = 5000,
  stopping_rounds = 10,
  stopping_tolerance = 0,
  seed = 123
)

# collect the results and sort by our model performance metric of choice
grid_perf <- h2o.getGrid(
  grid_id = "gbm_grid1", 
  sort_by = "mse", 
  decreasing = FALSE
)
grid_perf
proc.time()


#' \dontrun{
#' require(fda.usc)
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-phoneme[["classlearn"]]
#' mtest<-phoneme[["test"]]
#' gtest<-phoneme[["classtest"]]
#' dataf<-data.frame(glearn)
#' dat=list("df"=dataf,"x"=mlearn)
#' a1<-classif.glm(glearn~x, data = dat)
#' a2<-classif.svm(glearn~x, data = dat)
#' a3<-classif.svm2(glearn~x, data = dat)
#' newdat<-list("x"=mtest)
#' p1<-predict.classif(a1,newdat)
#' table(gtest,p1)
#' sum(p1==gtest)/250
#' }
#' @export


rp
 
 data(tecator)
 ab=tecator$absorp.fdata
 ab1=fdata.deriv(ab,nderiv=1)
 ab2=fdata.deriv(ab,nderiv=2)
 gfat<-cut(tecator$y$Fat,breaks=3)

 dat1<-data.frame(y=gfat,x=ab$data[,c(1,50,100)])
 
 a0<-rpart(y~.,data=dat1)
 predict(a0,dat1,type="prob")
 
 
 a0<-svm(y~.,data=dat1,probability=T,fitted=T)
 predict(a0,dat1,probability=T,decision.values=T)
 
 
 a4<- weighted.ksvm(y=gfat,x=ab$data[,c(1,50,100)],weights=rep(1,len=length(gfat)))
 names(a4)
 range(predict(a4,ab$data[,c(1,50,100)],type="benefit.score"))
 
 
 a1<-ksvm(y~.,dat1,prob.model=T)
 
 predict(a1,dat1,type = 'probabilities')
 predict(a1,dat1,type = 'votes')
 
 
 a0$decision.values[1:11,]
 
table(a0$fitted,gfat)
