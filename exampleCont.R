rm(list=ls())
source("_functionsCont.R") 

maxT=20
n=500
M=maxT+1
dpar=c(4,2,4,2)
P=matrix(c(0,1,1,0),2,2)

X=array(NA,c(n,maxT,2))
for(i in 1:n) {
X[i,,1]=rnorm(maxT)
X[i,,2]=rbinom(maxT,1,0.5)}

d=matrix(c(diff(pgamma(0:M,dpar[1],dpar[2])),
diff(pgamma(0:M,dpar[3],dpar[4]))),M,2)
d[M,]=1-apply(d[1:(M-1),],2,sum)

xi=matrix(-3,2,2)
xi[,2]=3

Y=array(NA,c(n,maxT,2))
for(i in 1:n) {
u=sample(2,1)
dw=sample(M,1,prob=d[,u])
tot=0
while(tot<maxT) {
for(j in (tot+1):min(tot+dw,maxT)) {
mu1=X[i,j,1]+X[i,j,2]+xi[1,u]
mu2=X[i,j,1]+X[i,j,2]+xi[2,u]
Y[i,j,1]=rnorm(1,mu1,1)
Y[i,j,2]=rnorm(1,mu2,1)}
tot=tot+dw
u=c(1,2)[-u]
dw=sample(M,1,prob=d[,u])}
}

cl=makeCluster(max(c(detectCores()-1,2)))

thetaInit=c(0,log(c(dpar)),c(xi,rep(1,4),c(0,0)))
mod=estHSMMpanelCont(Y,X,k=2,cl, init=thetaInit, Pinit=P, se=TRUE)

stopCluster(cl)

