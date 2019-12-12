Ker0.Kim=function(u){630*(4*u^2-1)^2*u^4*(abs(u)<=.5)} # c(-.5,.5)
Ker0.bi1=function(u){2*u^2*exp(-u^2)/sqrt(pi)}  # (-4,4)
Ker0.bi2=function(u){.5*abs(u)*exp(-abs(u))}    # (-10,10)

Qk=function(K,maxx=1){
  a=integrate(function(x){x^2*K(x)},-maxx,maxx)$value
  b=integrate(function(x){K(x)^2},-maxx,maxx)$value
  return(a*b^2)}
  muj=function(i=0,j=0,K,maxx=1){
  integrate(function(x){x^(i+j)*K(x)},-maxx,maxx)$value
}

Smat=function(p=3,K=Ker.epa,maxx=1){
  iaux=0:p
  jaux=0:p
  Smat=matrix(NA,nrow=length(iaux),ncol=length(jaux))
  for (i in 1:length(iaux)){for (j in 1:length(jaux)){Smat[i,j]=muj(iaux[i],jaux[j],K=K,maxx=maxx)}} 
  return(Smat)
}

Ker.star=function(x,K=Ker.epa,p=3,maxx=1){
  S=Smat(p=p,K=K,maxx=maxx)
  if (is.vector(x)) x=matrix(x,nrow=1)
  a=apply(x[1,,drop=FALSE],2,function(v){drop((c(1,rep(0,p))%*%solve(S)%*%v^(0:p)*K(v)))})
  if (nrow(x)>=2){
  for (i in 2:nrow(x)){
  a=rbind(a,apply(x[i,,drop=FALSE],2,function(v){drop((c(1,rep(0,p))%*%solve(S)%*%v^(0:p)*K(v)))}))
  }}
  if (nrow(x)==1) {a=a/sum(a)} else {a=sweep(a,1,rowSums(a),"/")}
  return(a)
}

CpK=function(K=Ker.epa,p=3,maxx=1){
  a=factorial(p+1)^2*integrate(function(x){Ker.star(x,K=K,p=p,maxx=maxx)^2},-maxx,maxx)$value
  b=2*(p+1)*(integrate(function(x){x^{p+1}*Ker.star(x,K=K,p=p,maxx=maxx)},-maxx,maxx)$value)^2
  return((a/b)^(1/(2*p+3)))
}

