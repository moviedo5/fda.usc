########################################################
cumgamma=function(D,w=min(nrow(D),ncol(D))){
n=nrow(D)
m=ncol(D)
DTW=matrix(0,nrow=n+1,ncol=m+1)
DTW[1,]=Inf
DTW[,1]=Inf
DTW[1,1]=0
for (i in 2:(n+1)){for (j in max(2,i-w):min(m+1,i+w)){DTW[i,j]=0}}
for (i in 2:(n+1)){
	for (j in max(2,i-w):min(m+1,i+w)){
	DTW[i,j]=D[i-1,j-1]+min(DTW[i-1,j],DTW[i,j-1],DTW[i-1,j-1])
	}}
return(DTW)
}
############################

############################
findPath=function(D,i=nrow(D),j=ncol(D)){
ca=c(i)
cb=c(j)
while(i>1 & j>1) {
if (i>1 & j>1){
aa=min(D[i-1,j-1],D[i-1,j],D[i,j-1])
if (aa==D[i-1,j-1]) {i=i-1;j=j-1} else if (aa==D[i-1,j]){i=i-1} else {j=j-1}
} else if (j>1){j=j-1} else {i=i-1}
ca=c(ca,i)
cb=c(cb,j)
}
aa=cbind(ca,cb)
aa=aa[-nrow(aa),]-1
return(aa)
}
############################
DTW=function(a,b,w=min(length(a),length(b)),p=2){
D=outer(a,b,"-")
D=abs(D)^p
DTW=cumgamma(D,w=w)
return(list(dist=DTW[length(a)+1,length(b)+1],D=DTW))
}
############################

############################
WDTW=function(a,b,w=length(a),p=2,wmax=1,g=0.05){
#g=0 (constante), 0.05 (linear), 0.25 (sigmoid), 3 (dos pesos)
n=length(b)
wei=outer(1:n,1:n,"-")
wei=wmax/(1+exp(-g*(abs(w)-n/2)))
M=outer(a,b,"-")
D=abs(wei*M)^p
DTW=cumgamma(D,w=w)
return(list(dist=DTW[length(a)+1,length(b)+1],D=DTW))
}
############################

############################
TWED=function(a,b,p=2,lambda=1,nu=1){
n=length(a)
m=length(b)
D=matrix(0,nrow=n+1,ncol=m+1)
D[1,1]=0
D[2,1]=abs(a[1])^p
D[1,2]=abs(b[1])^p
for (i in 3:(n+1)){ D[i,1]=D[i-1,1]+abs(a[i-2]-a[i-1])^p}
for (j in 3:(m+1)){ D[1,j]=D[1,j-1]+abs(b[i-2]-b[i-1])^p}
for (i in 2:(n+1)){
	for (j in 2:(m+1)){
	 if (i>2 & j>2) {d1=D[i-1,j-1]+2*nu*abs(i-j)+abs(a[i-1]-b[j-1])^p+abs(a[i-2]-b[j-2])^p} else {d1=D[i-1,j-1]+nu*abs(i-j)+abs(a[i-1]-b[j-1])^p}
	if (i>2) {d2=D[i-1,j]+abs(a[i-1]-a[i-2])^p+lambda+nu} else {d2=D[i-1,j]+abs(a[i-1])^p+lambda}
	if (j>2) {d3=D[i,j-1]+abs(b[j-1]-b[j-2])^p+lambda+nu} else {d3=D[i,j-1]+abs(b[j-1])^p+lambda}
	D[i,j]=min(d1,d2,d3)
	}}
return(list(dist=D[length(a)+1,length(b)+1],D=D))
}
############################
