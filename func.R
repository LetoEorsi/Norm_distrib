Gauss <-function(sigx,sigy,rho,N){

x1 <- runif(N,0,1)
x2 <- runif(N,0,1)

y1 <- sqrt(-2*log(x1))*cos(2*pi*x2)
y2 <- sqrt(-2*log(x1))*sin(2*pi*x2)

z1 = sigx*y1;
z2 = rho*sigy*y1+sigy*sqrt(1-rho^2)*y2
return(data.frame(xvar = z1,yvar = z2));
}

GaussU <-function(sigx,sigy,N){
  
  x1 <- runif(N,0,1)
  x2 <- runif(N,0,1)
  
  y1 <- sqrt(-2*log(x1))*cos(2*pi*x2)
  y2 <- sqrt(-2*log(x1))*sin(2*pi*x2)
  
  rho <- runif(N,-1,1)
  
  z1 = sigx*y1;
  z2 = rho*sigy*y1+sigy*sqrt(1-rho^2)*y2
  return(data.frame(xvar = z1,yvar = z2));
}

GaussR <-function(sigx,sigy,N){
  
  x1 <- runif(N,0,1)
  x2 <- runif(N,0,1)
  
  y1 <- sqrt(-2*log(x1))*cos(2*pi*x2)
  y2 <- sqrt(-2*log(x1))*sin(2*pi*x2)
  
  rho <- floor(runif(N,0,2))*2-1
  
  z1 = sigx*y1;
  z2 = rho*sigy*y1+sigy*sqrt(1-rho^2)*y2
  return(data.frame(xvar = z1,yvar = z2));
}

GaussB <-function(sigx,sigy,nu,mu,N){
  
  x1 <- runif(N,0,1)
  x2 <- runif(N,0,1)
  
  y1 <- sqrt(-2*log(x1))*cos(2*pi*x2)
  y2 <- sqrt(-2*log(x1))*sin(2*pi*x2)
  
  rho = 2*BetaN(nu,mu,N)-1
 
  
  z1 = sigx*y1;
  z2 = rho*sigy*y1+sigy*sqrt(1-rho^2)*y2
  return(data.frame(xvar = z1,yvar = z2));
}

SumInf <- function(dat){
  z1 <- dat$xvar
  z2 <- dat$yvar
  m1=mean(z1)
  m2=mean(z2)
  d1=var(z1)
  d2=var(z2)
  C=var(z1,z2)
  sin2 = 2*C/sqrt((d1-d2)^2+4*C^2)
  cos2 = (d1-d2)/sqrt((d1-d2)^2+4*C^2)
  tg=sin2/(1+cos2)
  
  costg=1/sqrt(1+tg^2)
  sintg=tg/sqrt(1+tg^2)
  
  r=sqrt((n/2)*(d1+d2-sqrt((d1-d2)^2+4*C^2)))
    
  return(data.frame(Mx=m1,My=m2,Dx=d1,Dy=d2,Cov=C,ry=c(C/d1,C/d2),rx=c(m2-C/d1*m1,m1-C/d2*m2),cy=tg,cx=m2-tg*m1,evp=c(costg,sintg),R=r))
}

Beta <- function(nu,mu){
  
  z1=runif(1,0,1)
  z2=runif(1,0,1)
  k=1
  
  while((k<100)&&((z1^(1/nu)+z2^(1/mu))>1 )){
    z1=runif(1,0,1)
    z2=runif(1,0,1)
    k=k+1
    if(k==100){
      return(-2);
    }
  }
  
  
  return((z1^(1/nu))/(z1^(1/nu)+z2^(1/mu)));
  
}

BetaN <- function(nu,mu,N){
  xi=rep(0,N)
  for(i in 1:N){
    xi[i]=-2
    while(xi[i]==-2){
      xi[i]=Beta(nu,mu);
    }
  }
  return(xi)
}





