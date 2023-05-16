function g=hmu(p,n,lam,X,mu)
  sigma=1/sqrt(n);
  Time=linspace(-1/lam,1/lam,100);
  result=0;
  for t=1:(length(Time)-1)
    result=result+mean(cos(Time(t)*(mu-X)))*exp(0.5*sigma^2*Time(t)^2)/lam/50;
  end
  result=result-2*(1-p)*sin(mu/lam)/mu;
  result=result/(2*pi*p);
  g=max(result,0);
end
