N=1000

u=rep(0,N)
sx=1
sy=1
r=0.5

if(sx == sy){
  tg=c(1,-1)
}else{
tg2=(2*sx*sy*r)/(sx^2-sy^2)
if(r>0){
  tg=c(-1/tg2+sqrt(1+1/tg2^2),-1/tg2-sqrt(1+1/tg2^2))  
}else{
  tg=c(-1/tg2-sqrt(1+1/tg2^2),-1/tg2+sqrt(1+1/tg2^2))
}
}
costg=1/sqrt(1+tg^2)
sintg=tg/sqrt(1+tg^2)

d <- Gauss(sx,sy,r,N+1)
z=rep(0,N)
for(n in 2:(N+1))
  {
  #dat <- SumInf(Gauss(sx,sy,0,1000))
  dat <- SumInf(head(d,n))
  z[n-1]=dat$R;
  #z[n-1]=sqrt( (((dat)$evp)[1] - costg[1])^2+(((dat)$evp)[2]-sintg[1])^2)
  #u[n-1]=z[n-1]*sqrt((n)/(log(log(max(exp(1),n)))))
                     
}
plot(z,type="l")
#plot(u,type = "s")
#s <- data.frame( k = c(rep("0.25",N+1),rep("0.5",N+1),rep("0.7",N+1),rep("1",N+1)),angle = c(z2,z1,z,z3))
#s <- data.frame(angle = z)
#ggplot(s,aes(x=angle))+
#geom_density()


