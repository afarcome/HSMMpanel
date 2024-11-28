library(parallel)
library(snipEM)
library(numDeriv)
library(compiler)
library(mhsmm)
library(snow)

lk=function(theta,Y,X,P,cl,k) {
p=dim(Y)[3]
q=dim(X)[3]
init=c(1,exp(theta[1:(k-1)]))
init=init/sum(init)
d=matrix(NA,M,k)
for(j in 1:k) {
d[,j]=diff(pgamma(0:M,exp(theta[k+(j-1)*2]),exp(theta[k+1+(j-1)*2])))
d[M,j]=1-sum(d[1:(M-1),j])}
xi=matrix(theta[k-1+k*2+1:(k*p)],p,k)
be=matrix(theta[k-1+k*2+k*p+1:(q*p)],p,q)
clusterExport(cl,c("Y","X","P","likn","xi","d","be","init"),envir=environment())
-sum(parSapply(cl,1:n,function(i) likn(init,d,xi,be,P,Y[i,,],X[i,,])))}

estHSMMpanelBin=function(Y,X,k,cl,init, Pinit, se=FALSE) {

p=dim(Y)[3]
q=dim(X)[3]
n=dim(X)[1]
maxT=dim(X)[2]
M=maxT+1

if(k==2) {
op=optim(init,lk,Y=Y,X=X,P=Pinit,cl=cl,k=k)
thetaInit=op$par
d=matrix(c(diff(pgamma(0:M,exp(thetaInit[2]),exp(thetaInit[3]))),
diff(pgamma(0:M,exp(thetaInit[4]),exp(thetaInit[5])))),M,2)
d[M,]=1-apply(d[1:(M-1),],2,sum)
init=c(1,exp(op$par[1]))/(1+exp(op$par[1]))
dpar=exp(thetaInit[2:5])
xi=matrix(thetaInit[5+1:(p*2)],p,2)
be=matrix(thetaInit[5+2*p+1:q*p],nrow=p)}

if(k>2) {



lkP=function(ps,th,k) {
P=matrix(0,k,k)
for(j in 1:k) {P[j,(1:k)[-j]]=exp(ps[1:(k-1)+(j-1)*(k-1)])
P[j,(1:k)[-j]]=P[j,(1:k)[-j]]/sum(P[j,(1:k)[-j]])}
lk(th,Y,X,P,cl,k)}

lp=NULL
for(j in 1:k) {lp=c(lp,log(P[j,(1:k)[-j]]))}

lik=lk(init,Y,X,P,cl,k)
likold=lik+1000
while(likold-lik>1e-3) {
op=optim(thetaInit,lk,Y=Y,X=X,P=P,cl=cl,k=k,control=list(maxit=100))
thetaInit=op$par
op2=optim(lp,lkP,th=op$par,k=k,control=list(maxit=100))
lp=op2$par

likold=lik
lik=op2$value
for(j in 1:k) {P[j,(1:k)[-j]]=exp(lp[1:(k-1)+(j-1)*(k-1)])
P[j,(1:k)[-j]]=P[j,(1:k)[-j]]/sum(P[j,(1:k)[-j]])}
}

thetaInit=op$par
init=c(1,exp(thetaInit[1:(k-1)]))
init=init/sum(init)
d=matrix(NA,M,k)
for(j in 1:k) {
d[,j]=diff(pgamma(0:M,exp(thetaInit[k+(j-1)*2]),exp(thetaInit[k+1+(j-1)*2])))
d[M,j]=1-sum(d[1:(M-1),j])}
xi=matrix(thetaInit[k-1+k*2+1:(k*p)],p,k)
be=matrix(thetaInit[k-1+k*2+k*p+1:(q*p)],p,q)

}

score=ses=info=NA 

if(se==TRUE) {
clusterExport(cl,c("forwn","gradn","dpar","hessn","Y","X","P","likn","xi","d","be","init"),envir=environment())
alphas=parSapply(cl,1:n,function(i) forwn(init,d,xi,be,P,Y[i,,],X[i,,]),simplify=FALSE)
clusterExport(cl,c("forwn","gradn","dpar","hessn","Y","X","P","likn","xi","d","be","init","alphas"),envir=environment())
grn=parSapply(cl,1:n,function(i) gradn(init,dpar,xi,be,P,Y[i,,],X[i,,],alphas[[i]]),simplify=FALSE)
clusterExport(cl,c("forwn","gradn","dpar","hessn","Y","X","P","likn","xi","d","be","init","alphas","grn"),envir=environment())
hs=parSapply(cl,1:n,function(i) hessn(init,dpar,xi,be,P,Y[i,,],X[i,,],alphas[[i]],grn[[i]]),simplify=FALSE)

score=sapply(1:n,function(i) apply(grn[[i]][maxT,,,],3,sum),simplify=FALSE)
likelihoods=sapply(1:n,function(i) sum(alphas[[i]][maxT,,]))
Js=sapply(1:n,function(i) apply(hs[[i]][maxT,,,,],3:4,sum),simplify=FALSE)

info=0
for(i in 1:n){
info=info+Js[[i]]/likelihoods[[i]]-(score[[i]]%*%t(score[[i]]))/(likelihoods[[i]]^2)}

ses=sqrt(diag(-solve(info)))}

return(list(pars=op$par,lik=-op$value,info=info,score=score, se=ses))
}

sumlog=function (x, lower = -745, upper = 709) 
{
    if (missing(x)) 
        stop("'x' missing")
    s <- tryCatch(.Call("fast_sumlog", x, lower, upper, length(x)), 
        `std::range_error` = function(e) {
            conditionMessage(e)
        })
    if (!is.finite(s)) return(max(x))
    return(s)
}

## forward probs

likn=function(init,d,xi,be,P,x,cov) {
library(snipEM)
p=ncol(x)
lx=nrow(x)
k=ncol(d)

if(nrow(d)<lx) {
d=rbind(d,matrix(rep(d[nrow(d),],(lx-nrow(d))),byrow=T,ncol=k))}

mu=as.vector(be%*%cov[1,])+xi

probs=array(NA,c(lx,lx,k))
probs[1,1,]=log(init)+apply(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T),2,sum)

for(ti in 2:lx) {

if(ti>2) {dwell=log(1-apply(d[1:(ti-1),],2,sum))} else {dwell=log(1-d[1,])}

lk=sapply(1:ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi
apply(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T),2,sum)})

probs[ti,ti,]=log(init)+dwell+apply(lk,1,sum)

for(s in 1:(ti-1)) {

if(s>1) {lks=apply(lk[,(ti-s+1):ti],1,sum)}
if(s==1) {lks=lk[,ti]}

if(ti-s>1) {jnk=apply(probs[ti-s,1:(ti-s),],2,sumlog)} else {jnk=probs[1,1,]}

for(u in 1:k) {probs[ti,s,u]=sumlog(log(P[-u,u])+jnk[-u])+log(d[s,u])+lks[u]}
}}

sumlog(probs[lx,,])}

forwn=function(init,d,xi,be,P,x,cov) {
library(snipEM)
p=ncol(x)
lx=nrow(x)
k=ncol(d)

if(nrow(d)<lx) {
d=rbind(d,matrix(rep(d[nrow(d),],(lx-nrow(d))),byrow=T,ncol=k))}

mu=as.vector(be%*%cov[1,])+xi

probs=array(NA,c(lx,lx,k))
probs[1,1,]=log(init)+apply(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T),2,sum)

for(ti in 2:lx) {

if(ti>2) {dwell=log(1-apply(d[1:(ti-1),],2,sum))} else {dwell=log(1-d[1,])}

lk=sapply(1:ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi
apply(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T),2,sum)})

probs[ti,ti,]=log(init)+dwell+apply(lk,1,sum)

for(s in 1:(ti-1)) {

if(s>1) {lks=apply(lk[,(ti-s+1):ti],1,sum)}
if(s==1) {lks=lk[,ti]}

if(ti-s>1) {jnk=apply(probs[ti-s,1:(ti-s),],2,sumlog)} else {jnk=probs[1,1,]}

for(u in 1:k) {probs[ti,s,u]=sumlog(log(P[-u,u])+jnk[-u])+log(d[s,u])+lks[u]}
}}

exp(probs)}

gradn=function(init,dpar,xi,be,P,x,cov,alphas) {
# discrete gamma dwelling
library(snipEM)
library(numDeriv)
p=ncol(x)
lx=nrow(x)
k=length(dpar)/2
d=matrix(NA,lx,k)
for(j in 1:k) {
d[,j]=diff(pgamma(0:lx,dpar[1+(j-1)*2],dpar[2+(j-1)*2]))}
d[lx,]=1-apply(d[1:(lx-1),],2,sum)
fd=function(theta) {
r=diff(pgamma(0:lx,theta[1],theta[2]))
r[lx]=1-sum(r[1:(lx-1)])
r
}
np=k-1+length(dpar)+length(xi)+length(be)
lb=length(be)

probs=array(0,c(lx,lx,k,np))

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
mu=as.vector(be%*%cov[1,])+xi[,u]
exp(sum(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T)))
}

# 1:(k-1) iniziali
# k*2 dwelling
# k*p xi
# k*q be 
# first x%*%(sz-x[1,]) <- + xi 
# numerically approximated 

# derivata rispetto a \pi_c
mu=as.vector(be%*%cov[1,])+xi[,k]
sk=exp(sum(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T)))

for(u in 1:(k-1)) {
mu=as.vector(be%*%cov[1,])+xi[,u]
gr=grad(f,c(xi,be),method="simple",u=u)
probs[1,1,u,u]=exp(sum(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T)))
probs[1,1,u,k-1+2*k+1:(k*p+lb)]=init[u]*gr
probs[1,1,k,u]=-sk}
gr=grad(f,c(xi,be),method="simple",u=k)
probs[1,1,k,k-1+2*k+1:(k*p+lb)]=init[k]*gr

for(ti in 2:lx) {

if(ti>2) {dwell=1-apply(d[1:(ti-1),],2,sum)} else {dwell=1-d[1,]}

lk=sapply(1:ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi
exp(apply(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T),2,sum))})

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
prod(sapply(1:ti,function(j) {
mu=as.vector(be%*%cov[j,])+xi[,u]
exp(sum(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T)))}))}

# derivata rispetto a \pi_c 

ap=apply(lk,1,prod)

for(u in 1:(k-1)) {probs[ti,ti,u,u]=ap[u]*dwell[u]
probs[ti,ti,k,u]=-ap[k]*dwell[k]}

for(u in 1:k) {
probs[ti,ti,u,k-1+(u-1)*2+1]=init[u]*ap[u]*grad(method="simple",function(theta) pgamma(ti,theta,dpar[2+(u-1)*2],lower.tail=F),dpar[1+(u-1)*2])
probs[ti,ti,u,k-1+(u-1)*2+2]=init[u]*ap[u]*grad(method="simple",function(theta) pgamma(ti,dpar[1+(u-1)*2],theta,lower.tail=F),dpar[2+(u-1)*2])

# derivata rispetto a xibe

probs[ti,ti,u,k-1+2*k+1:(k*p+lb)]=init[u]*pgamma(ti,dpar[1+(u-1)*2],dpar[2+(u-1)*2],lower.tail=F)*grad(method="simple",f,c(xi,be),u=u)

}

for(s in 1:(ti-1)) {

if(s>1) {lks=apply(lk[,(ti-s+1):ti],1,prod)}
if(s==1) {lks=lk[,ti]}

if(ti-s>1) {jnk=apply(probs[ti-s,1:(ti-s),,],2:3,sum)} else {jnk=probs[1,1,,]}
if(ti-s>1) {jnk2=apply(alphas[ti-s,1:(ti-s),],2,sum)} else {jnk2=alphas[1,1,]}

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
prod(sapply((ti-s+1):ti,function(j) {
mu=as.vector(be%*%cov[j,])+xi[,u]
exp(sum(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T)))}))}

# primo addendo (generale)
for(u in 1:k) {probs[ti,s,u,]=P[-u,u]*jnk[-u,]*d[s,u]*lks[u]}

# terzo addendo

for(u in 1:k) {
probs[ti,s,u,k-1+(u-1)*2+1]=probs[ti,s,u,k-1+(u-1)*2+1]+jnk2[-u]*lks[u]*grad(method="simple",function(theta) pgamma(s,theta,dpar[2+(u-1)*2])-pgamma(s-1,theta,dpar[2+(u-1)*2]),dpar[1+(u-1)*2])
probs[ti,s,u,k-1+(u-1)*2+2]=probs[ti,s,u,k-1+(u-1)*2+2]+jnk2[-u]*lks[u]*grad(method="simple",function(theta) pgamma(s,dpar[1+(u-1)*2],theta)-pgamma(s-1,dpar[1+(u-1)*2],theta),dpar[2+(u-1)*2])

# quarto addendo 

probs[ti,s,u,k-1+2*k+1:(k*p+lb)]=P[-u,u]*jnk2[-u]*d[s,u]*grad(method="simple",f,c(xi,be),u=u)}
}

}

probs}

hessn=function(init,dpar,xi,be,P,x,cov,alphas,grn) {
# discrete gamma dwelling
library(snipEM)
library(numDeriv)
p=ncol(x)
lx=nrow(x)
k=length(dpar)/2
d=matrix(NA,lx,k)
for(j in 1:k) {
d[,j]=diff(pgamma(0:lx,dpar[1+(j-1)*2],dpar[2+(j-1)*2]))}
d[lx,]=1-apply(d[1:(lx-1),],2,sum)
fd=function(theta) {
r=diff(pgamma(0:lx,theta[1],theta[2]))
r[lx]=1-sum(r[1:(lx-1)])
r
}
np=k-1+length(dpar)+length(xi)+length(be)
lb=length(be)

probs=array(0,c(lx,lx,k,np,np))

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
mu=as.vector(be%*%cov[1,])+xi[,u]
exp(sum(dbinom(x[1,],1,exp(mu)/(1+exp(mu)),log=T)))
}

# 1:(k-1) iniziali
# k*2 dwelling
# k*p xi
# k*q be 
# first x%*%(sz-x[1,]) <- + xi 
# second t(x)%*%diag(sz*(1-sz))%*%x 
# numerically approximated 

# derivata rispetto a \pi_c, e rispetto a opportuno \xi | \be
for(u in 1:(k-1)) {
gr=grad(f,c(xi,be),method="simple",u=u)
hs=hessian(f,c(xi,be),u=u)
probs[1,1,u,u,k-1+2*k+1:(k*p+lb)]=probs[1,1,u,k-1+2*k+1:(k*p+lb),u]=gr
# derivate doppie rispetto a \xi|\be
probs[1,1,u,k-1+2*k+1:(k*p+lb),k-1+2*k+1:(k*p+lb)]=init[u]*hs}

probs[1,1,k,k,k-1+2*k+1:(k*p+lb)]=probs[1,1,k,k-1+2*k+1:(k*p+lb),k]=-grad(f,c(xi,be),method="simple",u=k)
probs[1,1,k,k-1+2*k+1:(k*p+lb),k-1+2*k+1:(k*p+lb)]=init[k]*hessian(f,c(xi,be),u=k)

for(ti in 2:lx) {

if(ti>2) {dwell=1-apply(d[1:(ti-1),],2,sum)} else {dwell=1-d[1,]}

lk=sapply(1:ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi
exp(apply(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T),2,sum))})

# derivata rispetto a \pi_c e dwelling corrispondente (primo addendo) 

ap=apply(lk,1,prod)

gr1=grad(method="simple",function(theta) pgamma(ti,theta,dpar[2+(k-1)*2],lower.tail=F),dpar[1+(k-1)*2])
gr2=grad(method="simple",function(theta) pgamma(ti,dpar[1+(k-1)*2],theta,lower.tail=F),dpar[2+(k-1)*2])

probs[ti,ti,k,k-1+1:2+(k-1)*2,k-1+1:2+(k-1)*2]=init[k]*ap[k]*hessian(function(theta) pgamma(ti,theta[1],theta[2],lower.tail=F),dpar[(k-1)*2+1:2])
probs[ti,ti,k,k,k-1+1+(k-1)*2]=probs[ti,ti,k,k-1+1+(k-1)*2,k]=-ap[k]*gr1
probs[ti,ti,k,k,k-1+2+(k-1)*2]=probs[ti,ti,k,k-1+2+(k-1)*2,k]=-ap[k]*gr2

for(u in 1:(k-1)) {

probs[ti,ti,u,u,k-1+1+(u-1)*2]=probs[ti,ti,u,k-1+1+(u-1)*2,u]=ap[u]*grad(method="simple",function(theta) pgamma(ti,theta,dpar[2+(u-1)*2],lower.tail=F),dpar[1+(u-1)*2])
probs[ti,ti,u,u,k-1+2+(u-1)*2]=probs[ti,ti,u,k-1+2+(u-1)*2,u]=ap[u]*grad(method="simple",function(theta) pgamma(ti,dpar[1+(u-1)*2],theta,lower.tail=F),dpar[2+(u-1)*2])
probs[ti,ti,u,u,k-1+(k-1)*2+1]=probs[ti,ti,u,k-1+(k-1)*2+1,u]=-ap[u]*gr1
probs[ti,ti,u,u,k-1+(k-1)*2+2]=probs[ti,ti,u,k-1+(k-1)*2+2,u]=-ap[u]*gr2
probs[ti,ti,u,k-1+1:2+(u-1)*2,k-1+1:2+(u-1)*2]=init[u]*ap[u]*hessian(function(theta) pgamma(ti,theta[1],theta[2],lower.tail=F),dpar[(u-1)*2+1:2])
}

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
sum(sapply(1:ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi[,u]
sum(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T))}))}

# secondo addendo
gr=grad(method="simple",function(theta) exp(f(theta,u=k)),c(xi,be))
probs[ti,ti,k,k,k-1+2*k+1:(k*p+lb)]=probs[ti,ti,k,k-1+2*k+1:(k*p+lb),k]=-pgamma(ti,dpar[1+(k-1)*2],dpar[2+(k-1)*2],lower.tail=F)*gr

for(u in 1:k) {
gr=grad(method="simple",function(theta) exp(f(theta,u=u)),c(xi,be))
hs=hessian(function(theta) exp(f(theta,u=u)),c(xi,be))

if(u<k) {probs[ti,ti,u,u,k-1+2*k+1:(k*p+lb)]=probs[ti,ti,u,k-1+2*k+1:(k*p+lb),u]=pgamma(ti,dpar[1+(u-1)*2],dpar[2+(u-1)*2],lower.tail=F)*gr}
# quarto addendo
probs[ti,ti,u,k-1+1+(u-1)*2,k-1+2*k+1:(k*p+lb)]=probs[ti,ti,u,k-1+2*k+1:(k*p+lb),k-1+1+(u-1)*2]=gr*init[u]*grad(method="simple",function(theta) pgamma(ti,theta,dpar[2+(u-1)*2],lower.tail=F),dpar[1+(u-1)*2])
probs[ti,ti,u,k-1+2+(u-1)*2,k-1+2*k+1:(k*p+lb)]=probs[ti,ti,u,k-1+2*k+1:(k*p+lb),k-1+2+(u-1)*2]=gr*init[u]*grad(method="simple",function(theta) pgamma(ti,dpar[1+(u-1)*2],theta,lower.tail=F),dpar[2+(u-1)*2])

# quinto addendo

probs[ti,ti,u,k-1+2*k+1:(k*p+lb),k-1+2*k+1:(k*p+lb)]=init[u]*pgamma(ti,dpar[1+(u-1)*2],dpar[2+(u-1)*2],lower.tail=F)*hs 
}

for(s in 1:(ti-1)) {

if(s>1) {lks=apply(lk[,(ti-s+1):ti],1,prod)}
if(s==1) {lks=lk[,ti]}

if(ti-s>1) {jnk=apply(probs[ti-s,1:(ti-s),,,],2:4,sum)
grn1=apply(grn[ti-s,1:(ti-s),,],2:3,sum)
a1=apply(alphas[ti-s,1:(ti-s),],2,sum)
} else {jnk=probs[1,1,,,]
grn1=grn[1,1,,]
a1=alphas[1,1,]}

f=function(theta,u) {
xi=matrix(theta[1:(k*p)],p,k)
be=matrix(theta[-c(1:(k*p))],nrow=p)
sum(sapply((ti-s+1):ti,function(j)
{
mu=as.vector(be%*%cov[j,])+xi[,u]
sum(dbinom(x[j,],1,exp(mu)/(1+exp(mu)),log=T))}))}


for(u in 1:k) {

gr=grad(method="simple",function(theta) exp(f(theta,u=u)),c(xi,be))
hs=hessian(function(theta) exp(f(theta,u=u)),c(xi,be))

# primo addendo

probs[ti,s,u,,]=P[-u,u]*jnk[-u,,]*d[s,u]*lks[u]

# dwelling doppia

probs[ti,s,u,k-1+2*(u-1)+1:2,k-1+2*(u-1)+1:2]=probs[ti,s,u,k-1+2*(u-1)+1:2,k-1+2*(u-1)+1:2]+hessian(function(theta) pgamma(s,theta[1],theta[2])-pgamma(s-1,theta[1],theta[2]),dpar[(u-1)*2+1:2])*P[-u,u]*lks[u]*a1[-u]

# terzo (con symmetry)

probs[ti,s,u,,k-1+1+(u-1)*2]=probs[ti,s,u,k-1+1+(u-1)*2,]=probs[ti,s,u,k-1+1+(u-1)*2,]+grn1[-u,]*P[-u,u]*lks[u]*grad(method="simple",function(theta) pgamma(s,theta,dpar[2+(u-1)*2])-pgamma(s-1,theta,dpar[2+(u-1)*2]),dpar[1+(u-1)*2])
probs[ti,s,u,,k-1+2+(u-1)*2]=probs[ti,s,u,k-1+2+(u-1)*2,]=probs[ti,s,u,k-1+1+(u-1)*2,]+grn1[-u,]*P[-u,u]*lks[u]*grad(method="simple",function(theta) pgamma(s,dpar[1+(u-1)*2],theta)-pgamma(s-1,dpar[1+(u-1)*2],theta),dpar[2+(u-1)*2])

# quarto (con symmetry)

probs[ti,s,u,,k-1+2*k+1:(k*p+lb)]=probs[ti,s,u,,k-1+2*k+1:(k*p+lb)]+d[s,u]*P[-u,u]*(grn1[-u,]%*%t(gr))
probs[ti,s,u,k-1+2*k+1:(k*p+lb),]=t(probs[ti,s,u,,k-1+2*k+1:(k*p+lb)])

# decimo (con symmetry) - dwelling & manifesta

probs[ti,s,u,k-1+2*k+1:(k*p+lb),k-1+1+(u-1)*2]=probs[ti,s,u,k-1+1+(u-1)*2,k-1+2*k+1:(k*p+lb)]=probs[ti,s,u,k-1+1+(u-1)*2,k-1+2*k+1:(k*p+lb)]+a1[-u]*P[-u,u]*grad(method="simple",function(theta) pgamma(s,theta,dpar[2+(u-1)*2])-pgamma(s-1,theta,dpar[2+(u-1)*2]),dpar[1+(u-1)*2])*gr

probs[ti,s,u,k-1+2*k+1:(k*p+lb),k-1+2+(u-1)*2]=probs[ti,s,u,k-1+2+(u-1)*2,k-1+2*k+1:(k*p+lb)]=probs[ti,s,u,k-1+2+(u-1)*2,k-1+2*k+1:(k*p+lb)]+a1[-u]*P[-u,u]*grad(method="simple",function(theta) pgamma(s,dpar[1+(u-1)*2],theta)-pgamma(s-1,dpar[1+(u-1)*2],theta),dpar[2+(u-1)*2])*gr

# undicesimo (senza symmetry) - doppia manifesta

probs[ti,s,u,k-1+2*k+1:(k*p+lb),k-1+2*k+1:(k*p+lb)]=probs[ti,s,u,k-1+2*k+1:(k*p+lb),k-1+2*k+1:(k*p+lb)]+a1[-u]*P[-u,u]*d[s,u]*hs

}

}}

probs}

likn=cmpfun(likn)
forwn=cmpfun(forwn)
gradn=cmpfun(gradn)
hessn=cmpfun(hessn)

