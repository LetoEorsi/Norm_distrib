dat <- Gauss(1,1,0.5,10000)
S <- SumInf(dat)

library(ggplot2)
ggplot(dat, aes(x=xvar, y=yvar)) +
  geom_point(shape=1)+   
  geom_abline(intercept = S$cx, slope = S$cy,color="red")+
  geom_abline(intercept = (S$rx)[1], slope = (S$ry)[1] ,color="blue")+
  geom_abline(intercept = -(S$rx)[2]/(S$ry)[2], slope = 1/(S$ry)[2] ,color="magenta")
