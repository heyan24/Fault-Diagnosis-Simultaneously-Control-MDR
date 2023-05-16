library(decon)
library(pheatmap)
library(ggplot2)
library(scatterplot3d)
library(RColorBrewer)
#setwd("~/Dropbox/diag_MDRmatlab/R")

#------------------------------------- data preprocessing ------------------------------------------

y=read.table(file = 'secom_labels.data')
y=y[,1]
x=read.table(file = 'secom.data')# the origin dataset


for (i in 1:length(x)){
  x[is.na(x[,i]),i]=mean(x[,i],na.rm = T)
} # deal with the missing values
count=0
cons=c()
for (i in 1:length(x)){
  if (length(unique(x[,i]))==1){count=count+1
  cons=c(cons,i)}
}
count

x=x[-cons] # find the constant variables and remove them from the dataset
x=x[-c(69,108,189,192,218,288,293,321,389,417)]
x=x[-c(79,147,148,199,209,243,244,297,309,393,403)]
#x=as.matrix(x)
x_ic=x[y==-1,];x_oc=x[y==1,]#get the IC and OC datasets
p=c()
for (i in 1:length(x_ic)){
  p=c(p,shapiro.test(x_ic[,i])$p.value)
}
sum(p<=0.01)#normality test

xx_ic=x_ic
for (i in 1:length(x_ic)){
  Fn=ecdf(x_ic[,i])
  xx_ic[,i]=qnorm(Fn(x_ic[,i]))
}
index=xx_ic==Inf
xx_ic[index]=3.75

pp=c()
for (i in 1:length(xx_ic)){
  pp=c(pp,shapiro.test(xx_ic[is.finite(xx_ic[,i]),i])$p.value)
}
sum(pp<=0.01)#normality test for IC data 

COV=cor(xx_ic)
pdf('cov.pdf',width=6,height=5)
pheatmap(COV[1:100,1:100],cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color='grey')
dev.off()

xx_oc=x_oc
for (i in 1:length(x_oc)){
  Fn=ecdf(x_ic[,i])
  xx_oc[,i]=qnorm(Fn(x_oc[,i]))
}
index=xx_oc==Inf
xx_oc[index]=3.75
index=xx_oc==-Inf
xx_oc[index]=-3.75
# p=c()
# for (i in 1:length(xx_oc)){
#   p=c(p,shapiro.test(xx_ic[,i])$p.value)
# }
# sum(p<=0.01)#normality test
#meanXIC=apply(xx_ic[1:1463,],MARGIN=2,FUN=mean)

#-------------------------------------------- Data-driven procedure -------------------------------------------

N=c(1,2,5,10,15,20,25,30)
num=c()
for (n in N){
  #n=104
  n=104
  meanX=apply(xx_oc[1:n,],MARGIN=2,FUN=mean)
  m=length(meanX)
  #para=ESTNULL(meanX)
  mu=0;sigma=1/sqrt(n)
  #lam=bw.dboot1(meanX,sigma)
  phat=EPSIEST(meanX,mu,sigma); print(phat) # estimate proportions
  
  MU=seq(-5,5,length.out = 500);MU1=MU[251:499];MU2=MU[1:250]
  f=density(as.vector(meanX))
  f<- approxfun(f$x, f$y)
  fxi=f(MU) # f(xi)
  fxi[is.na(fxi)]=0
  plot(MU,fxi,type='l',xlim = c(-0.75,1),ylim=c(0,4),ylab='Density',xlab='Index',lwd=2)
  lam=0.08 # choice of lambda
  GMU=rep(0,length(MU-1))
  for (i in 1:length(GMU)){
    GMU[i]=gmu2((1-phat$p0),sigma,lam,meanX,MU[i]) # estiamte h(mu_i)
  }
  GMU1=GMU[251:499]*(1-phat$p0)/phat$p1 # estimate h1(mu_i) and h2(mu_i)
  GMU2=GMU[1:250]*(1-phat$p0)/phat$p2
  
  pdf("hmu.pdf",width=10,height=6)
  plot(MU,GMU,type='l',xlim = c(-0.75,1),ylim=c(0,4),ylab='Density',xlab=expression(mu),lwd=3)
  #lines(MU,GMU,type='l',lwd=2,lty=2)
  lines(MU[249:499],c(0.4,0.8,GMU1),type='l',col='red',lty=2,lwd=3)
  lines(MU2[-(249:250)],GMU2[-(249:250)],type='l',col='blue',lty=4,lwd=3)
  #lines(MU1,g1,type='l',col='green')
  legend('topright',c(expression(hat(h)(mu)),expression(hat(h)[1](mu)),expression(hat(h)[2](mu))),col=c('black','red','blue'),lty=c(1,2,4))
  dev.off()
  
  g1=rep(0,m);g2=g1 # estimate conditional distribution of X under different thetas
  for (j in 1:m){
    x=meanX[j]
    g1[j]=sum(dnorm(x,MU1,sigma)*GMU1)/50
    g2[j]=sum(dnorm(x,MU2,sigma)*GMU2)/50
  }
  g0=dnorm(meanX,0,sigma)
  ######## 下面用EM算法估计 ############
  Int <- matrix(c(
    0.8,0.1,0.1,
    0.8,0.1,0.1,
    0.7,0.2,0.1,
    0.7,0.1,0.2
  ),nrow = 4, byrow = T)
  em <- EM(m,g0,g1,g2,Int,500) #estimate parameters in HMM
  pi0_e=em$pi0_e; pi1_e=em$pi1_e; pi2_e=em$pi2_e; A_e=em$A_e
  
  c0=rep(0,m)
  Alpha0=rep(0,m); Alpha1=rep(0,m); Alpha2=rep(0,m)
  Beta0=rep(0,m); Beta1=rep(0,m); Beta2=rep(0,m)
  Alpha0[1]=pi0_e*g0[1]; Alpha1[1]=pi1_e*g1[1]; Alpha2[1]=pi2_e*g2[1]
  c0[1]=1/(Alpha0[1]+Alpha1[1]+Alpha2[1])
  Alpha0[1]=Alpha0[1]*c0[1]; Alpha1[1]=Alpha1[1]*c0[1]; Alpha2[1]=Alpha2[1]*c0[1]; 
  Beta0[m]=1; Beta1[m]=1; Beta2[m]=1
  
  for (i in 2:m){
    Alpha0[i] <- g0[i]*(Alpha0[i-1]*A_e[1,1]+Alpha1[i-1]*A_e[2,1]+Alpha2[i-1]*A_e[3,1])
    Alpha1[i] <- g1[i]*(Alpha0[i-1]*A_e[1,2]+Alpha1[i-1]*A_e[2,2]+Alpha2[i-1]*A_e[3,2])
    Alpha2[i] <- g2[i]*(Alpha0[i-1]*A_e[1,3]+Alpha1[i-1]*A_e[2,3]+Alpha2[i-1]*A_e[3,3])
    c0[i] <- 1/(Alpha0[i]+Alpha1[i]+Alpha2[i])
    Alpha0[i]=Alpha0[i]*c0[i]; Alpha1[i]=Alpha1[i]*c0[i]; Alpha2[i]=Alpha2[i]*c0[i]
    
  }
  
  Beta0[m]=Beta0[m]*c0[m]; Beta1[m]=Beta1[m]*c0[m]; Beta2[m]=Beta2[m]*c0[m]
  for (i in 2:m){
    Beta0[m-i+1] <- (A_e[1,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[1,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[1,3]*g2[m-i+2]*Beta2[m-i+2])*c0[m-i+1]
    Beta1[m-i+1] <- (A_e[2,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[2,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[2,3]*g2[m-i+2]*Beta2[m-i+2])*c0[m-i+1]
    Beta2[m-i+1] <- (A_e[3,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[3,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[3,3]*g2[m-i+2]*Beta2[m-i+2])*c0[m-i+1]
   
  }
  H0 <- Alpha0*Beta0
  H1 <- Alpha1*Beta1
  H2 <- Alpha2*Beta2
  H <- H0+H1+H2
  H0 <- H0/H; H1 <- H1/H; H2 <- H2/H
  
  alpha1 <- 0.1; alpha2 <- 0.05 #setting MDR level for both directions
  
  LAMBDA1 <- H1/H0; LAMBDA2 = H2/H0
  LAM1_sort <- sort(LAMBDA1); LAM2_sort <- sort(LAMBDA2)
  Ind1 <- sort(LAMBDA1, index.return=T)$ix
  Ind2 <- sort(LAMBDA2, index.return=T)$ix
  H1_sort <- H1[Ind1]; H2_sort <- H2[Ind2]
  result1 <- rep(0,m); result2 <- rep(0,m)
  for (j in 1:m){
    result1[j] <- sum(H1_sort[j:m]) / sum(H1)
    result2[j] <- sum(H2_sort[j:m]) / sum(H2)
  }
  r1 <- which(result1-1+alpha1<=0)[1] - 1
  r2 <- which(result2-1+alpha2<=0)[1] - 1
  lam1 <- 1/LAM1_sort[r1]; lam2 <- 1/LAM2_sort[r2]
  
  N_lam1 <- Nhat(1, DRule(c(lam1,lam2),H0,H1,H2),H1,alpha1)
  N_lam2 <- Nhat(-1, DRule(c(lam1,lam2),H0,H1,H2),H2,alpha2)
  N_tilde1 <- Nhat(1,DRule(c(1/LAM1_sort[r1-1],lam2),H0,H1,H2),H1,alpha1)
  N_tilde2 <- Nhat(-1,DRule(c(lam1,1/LAM2_sort[r2-1]),H0,H1,H2),H2,alpha2)
  step <- 1
  
  if (N_lam1 >= 0 & N_lam2 >= 0){
    while (!(N_lam1>=0 & N_lam2>=0 & N_tilde1<0 & N_tilde2<0)){
      step <- step + 1
      N_pot1 <- rep(0,r1); N_pot2 <- rep(0,r2)
      for (j in 1:r1){
        N_pot1[j] <- Nhat(1,DRule(c(1/LAM1_sort[j],lam2),H0,H1,H2),H1,alpha1)
      }
      for (j in 1:r2){
        N_pot2[j] <- Nhat(-1,DRule(c(lam1,1/LAM2_sort[j]),H0,H1,H2),H2,alpha2)
      }
      r1 <- which(N_pot1 > 0)[1]; r2 <- which(N_pot2 > 0)[1]
      lam1 <- 1/LAM1_sort[r1]
      lam2 <- 1/LAM2_sort[r2]
      
      if (r1 == 1 | r2 == 1){
        break
      }
      
      N_lam1 <- Nhat(1, DRule(c(lam1,lam2),H0,H1,H2),H1,alpha1)
      N_lam2 <- Nhat(-1, DRule(c(lam1,lam2),H0,H1,H2),H2,alpha2)
      N_tilde1 <- Nhat(1,DRule(c(1/LAM1_sort[r1-1],lam2),H0,H1,H2),H1,alpha1)
      N_tilde2 <- Nhat(-1,DRule(c(lam1,1/LAM2_sort[r2-1]),H0,H1,H2),H2,alpha2)
      
    }
  }
  
  lam_star <- c(lam1, lam2)
  delta <- DRule(lam_star, H0, H1, H2)
  
  
  num=c(num,cutt)
  print(c(n,cutt))
}
col=rep(1,m)
col[which(delta == 1)]=2;col[which(delta == -1)]=3
plot(meanX,cex=0.8,col=col)
pdf("../figures/plot_real2.pdf",width=16,height=10)
df <- data.frame(Index = 1:453, Value = (meanX), Status=factor(col,labels = c('IC','positive shifts','negative shifts')))
ggplot(data = df, mapping = aes(x = Index, y = Value,colour = Status,shape=Status)) +
  geom_point(size = 6) +
  theme(legend.position='top',legend.title = element_text(size = 20),legend.text = element_text(size = 20,face = 'bold'),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15)) +
  scale_color_brewer(type='qual',palette = 2) +
  labs(x = 'Component i', y = expression(paste('Value ' ,(X[i]), sep=' ')))
dev.off()

# alpha1 = alpha2 = 0.1: sum(delta == 1) = 248; sum(delta == -1) = 17

# alpha1 = 0.1, alpha2 = 0.05: sum(delta == 1) = 253; sum(delta == -1) = 74



#------------------------------ Functions ---------------------------------------

gmu2=function(p,sigma,lam,X,mu){
  #sigma=1/sqrt(n)
  Time=seq(-1/lam,1/lam,length.out = 100)
  result=0
  for (t in 1:(length(Time)-1)){
    result=result+mean(cos(Time[t]*(mu-X)))*exp(0.5*sigma^2*Time[t]^2)/lam/50
  }
  result=result-2*(1-p)*sin(mu/lam)/mu
  result=result/(2*pi*p)
  return(max(result,0))
}

ESTNULL=function(x,gamma=0.1){
  n=length(x)
  t=seq(0,5,0.005)
  T=length(t)
  gan=n^(-gamma)
  phiplus=rep(1,T)
  phiminus=rep(1,T)
  dphiplus=rep(1,T)
  dphiminus=rep(1,T)
  phi=rep(1,T)
  for (k in 1:T){
  s=t[k]
  phiplus[k]=mean(cos(s*x))
  phiminus[k]=mean(sin(s*x))
  dphiplus[k]=-mean(x*sin(s*x))
  dphiminus[k]=mean(x*cos(s*x))
  phi[k]=sqrt(phiplus[k]^2 + phiminus[k]^2)
  }
  TT=1:T
  ind=min(TT[phi-gan<=0])
  tt=t[ind]
  a=phiplus[ind]
  b=phiminus[ind]
  da=dphiplus[ind]
  db=dphiminus[ind]
  c=phi[ind]
  
  u=-(da*b-db*a)/(c*c)
  sigma=sqrt(-(a*da+b*db)/(tt*c*c))
  result = list(mu=u,sigma=sigma)
  return(result)
}

EPSIEST <- function(x,u,sigma){
  z=(x-u)/sigma
  xi= seq(0,1,0.01)
  gamma=0.5
  tmax=sqrt(2*gamma*log(length(x)))
  tt=seq(0,tmax,length.out = 100)
  epshat=rep(0,length(tt))
  for (k1 in 1:length(tt)){
    t=tt[k1]
    f=t*xi
    f=exp(f^2/2)
    w=(1-abs(xi))
    co=rep(0,length(xi))
    for (k2 in 1:length(xi)) {
      co[k2]=mean(cos(t*xi[k2]*z))
    }
    epshat[k1]=1-sum(w*f*co)/sum(w)
  }
  p0=1-max(epshat)
  ind=sort.int(abs(x),index.return = TRUE)$ix
  x_sort=x[ind]
  m=453
  x_oc=x_sort[(p0*m):m]
  p1=sum(x_oc>0)/m
  p2=1-p0-p1
  results=list(p0=p0,p1=p1,p2=p2)
  return(results)
}

EM <- function(m,g0,g1,g2,Int,limit){
  pi0_e=Int[1,1]; pi1_e=Int[1,2]; pi2_e=Int[1,3]
  A_e=matrix(c(
    Int[2,1],Int[2,2],Int[2,3],
    Int[3,1],Int[3,2],Int[3,3],
    Int[4,1],Int[4,2],Int[4,3]),
    nrow=3, byrow=T
  )
  l=0
  #limit=50
  while (l<limit){
    c0=rep(0,m)
    d0=rep(0,m)
    Alpha0=rep(0,m); Alpha1=rep(0,m); Alpha2=rep(0,m)
    Beta0=rep(0,m); Beta1=rep(0,m); Beta2=rep(0,m)
    Alpha0[1]=pi0_e*g0[1]; Alpha1[1]=pi1_e*g1[1]; Alpha2[1]=pi2_e*g2[1]
    c0[1]=1/(Alpha0[1]+Alpha1[1]+Alpha2[1])
    Alpha0[1]=Alpha0[1]*c0[1]; Alpha1[1]=Alpha1[1]*c0[1]; Alpha2[1]=Alpha2[1]*c0[1]; 
    Beta0[m]=1; Beta1[m]=1; Beta2[m]=1
    d0[m]=1/(Beta0[m]+Beta1[m]+Beta2[m])
    
    for (i in 2:m){
      Alpha0[i] <- g0[i]*(Alpha0[i-1]*A_e[1,1]+Alpha1[i-1]*A_e[2,1]+Alpha2[i-1]*A_e[3,1])
      Alpha1[i] <- g1[i]*(Alpha0[i-1]*A_e[1,2]+Alpha1[i-1]*A_e[2,2]+Alpha2[i-1]*A_e[3,2])
      Alpha2[i] <- g2[i]*(Alpha0[i-1]*A_e[1,3]+Alpha1[i-1]*A_e[2,3]+Alpha2[i-1]*A_e[3,3])
      c0[i] <- 1/(Alpha0[i]+Alpha1[i]+Alpha2[i])
      Alpha0[i]=Alpha0[i]*c0[i]; Alpha1[i]=Alpha1[i]*c0[i]; Alpha2[i]=Alpha2[i]*c0[i]
      
    }
    
    Beta0[m]=Beta0[m]*d0[m]; Beta1[m]=Beta1[m]*d0[m]; Beta2[m]=Beta2[m]*d0[m]
    for (i in 2:m){
      Beta0[m-i+1] <- A_e[1,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[1,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[1,3]*g2[m-i+2]*Beta2[m-i+2]
      Beta1[m-i+1] <- A_e[2,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[2,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[2,3]*g2[m-i+2]*Beta2[m-i+2]
      Beta2[m-i+1] <- A_e[3,1]*g0[m-i+2]*Beta0[m-i+2]+A_e[3,2]*g1[m-i+2]*Beta1[m-i+2]+A_e[3,3]*g2[m-i+2]*Beta2[m-i+2]
      d0[m-i+1] <- 1/(Beta0[m-i+1]+Beta1[m-i+1]+Beta2[m-i+1])
      Beta0[m-i+1]=Beta0[m-i+1]*d0[m-i+1]; Beta1[m-i+1]=Beta1[m-i+1]*d0[m-i+1]; Beta2[m-i+1]=Beta2[m-i+1]*d0[m-i+1]
    }
    
    r0=Alpha0*Beta0; r1=Alpha1*Beta1; r2=Alpha2*Beta2
    s=r0+r1+r2
    r0=r0/s; r1=r1/s; r2=r2/s
    k00=rep(0,m-1); k01=rep(0,m-1); k02=rep(0,m-1)
    k10=rep(0,m-1); k11=rep(0,m-1); k12=rep(0,m-1)
    k20=rep(0,m-1); k21=rep(0,m-1); k22=rep(0,m-1)
    for (i in 1:(m-1)){
      k00[i] <- A_e[1,1]*r0[i]*g0[i+1]*Beta0[i+1]*d0[i]/Beta0[i]
      k01[i] <- A_e[1,2]*r0[i]*g1[i+1]*Beta1[i+1]*d0[i]/Beta0[i]
      k02[i] <- A_e[1,3]*r0[i]*g2[i+1]*Beta2[i+1]*d0[i]/Beta0[i]
      k10[i] <- A_e[2,1]*r1[i]*g0[i+1]*Beta0[i+1]*d0[i]/Beta1[i]
      k11[i] <- A_e[2,2]*r1[i]*g1[i+1]*Beta1[i+1]*d0[i]/Beta1[i]
      k12[i] <- A_e[2,3]*r1[i]*g2[i+1]*Beta2[i+1]*d0[i]/Beta1[i]
      k20[i] <- A_e[3,1]*r2[i]*g0[i+1]*Beta0[i+1]*d0[i]/Beta2[i]
      k21[i] <- A_e[3,2]*r2[i]*g1[i+1]*Beta1[i+1]*d0[i]/Beta2[i]
      k22[i] <- A_e[3,3]*r2[i]*g2[i+1]*Beta2[i+1]*d0[i]/Beta2[i]
    }
    pi0_new=r0[1]; pi1_new=r1[1]; pi2_new=r2[1]
    A_new=matrix(c(
      sum(k00)/sum(r0[1:m-1]),sum(k01)/sum(r0[1:m-1]),sum(k02)/sum(r0[1:m-1]),
      sum(k10)/sum(r1[1:m-1]),sum(k11)/sum(r1[1:m-1]),sum(k12)/sum(r1[1:m-1]),
      sum(k20)/sum(r2[1:m-1]),sum(k21)/sum(r2[1:m-1]),sum(k22)/sum(r2[1:m-1])
    ),nrow = 3,byrow = T)
    
    pi0_e=pi0_new; pi1_e=pi1_new; pi2_e=pi2_new
    A_e=A_new
    l=l+1
    #result=list(pi0_e=pi0_e, pi1_e=pi1_e, pi2_e=pi2_e, A_e=A_e)
    #return(result)
  }
  result=list(pi0_e=pi0_e, pi1_e=pi1_e, pi2_e=pi2_e, A_e=A_e)
  return(result)
  
}

stable <- function(A) {
  Q <- A
  Q[1,3] <- 1; Q[2,3] <- 1; Q[3,3] <- 1
  Q[1,1] <- A[1,1] - 1; Q[2,2] <- A[2,2] - 1
  C <- matrix(c(0,0,1), nrow=3)
  p <- solve(t(Q)) %*% C
  return(p)
}

Nhat <- function(k, delta, H, alpha) {
  N = sum((1 - 1 * (delta==k) - alpha) * H)
  return(N)
}

DRule <- function(lambda, H0, H1, H2) {
  m <- length(H0)
  delta <- rep(0,m)
  R1 <- H0 - lambda[1] * H1
  R2 <- H0 - lambda[2] * H2
  delta[R1 <= 0 & R1 <= R2] = 1
  delta[R2 <= 0 & R2 <= R1] = -1
  return(delta)
}
