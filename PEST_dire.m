%This function is for estimating the non-null proportion of guassian
%mixture distribution
%x is the dataset
% u and sigma is the known or estimated parameters of null distribution
function [p0,p1,p2]=PEST_dire(x,u,sigma)
z=(x-u)/sigma;
xi= 0:0.01:1;
gamma=0.5;
tmax=sqrt(2*gamma*log(length(x)));
tt=linspace(0,tmax,100);
epshat=zeros(1,length(tt));
for k1=1:length(tt)
    t=tt(k1);
    f=t*xi;
    f=exp(f.^2/2);
    w=(1-abs(xi));
    co=0.*xi;
for k2=1:length(xi)
    co(k2)=mean(cos(t*xi(k2).*z));
end
    epshat(k1)=1-sum(w.*f.*co)/sum(w);
end
p0=1-max(epshat);
[~,ind]=sort(abs(x));
x_sort=x(ind);
m=length(x);
x_oc=x_sort(p0*m+1:m);
p1=sum(x_oc>0)/m;p2=1-p0-p1;%sum(x_oc<0)/m;
end