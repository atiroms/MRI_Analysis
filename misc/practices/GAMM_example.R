library(mgcv)
## simple examples using gamm as alternative to gam
set.seed(0) 
n <- 400
sig <- 2
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
e <- rnorm(n, 0, sig)
y <- f + e
b <- gamm(y~s(x0)+s(x1)+s(x2)+s(x3))
plot(b$gam,pages=1)
summary(b$lme) # details of underlying lme fit
summary(b$gam) # gam style summary of fitted model
anova(b$gam) 

b <- gamm(y~te(x0,x1)+s(x2)+s(x3)) 
op <- par(mfrow=c(2,2))
plot(b$gam)
par(op)

## Add a factor to the linear predictor, to be modelled as random
fac <- rep(1:4,n/4)
f <- f + fac*3
fac<-as.factor(fac)

g<-exp(f/5)
y<-rpois(rep(1,n),g)
b2<-gamm(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson,random=list(fac=~1))
plot(b2$gam,pages=1)

## now an example with autocorrelated errors....
x <- 0:(n-1)/(n-1)
f <- 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10-1.396
e <- rnorm(n,0,sig)
for (i in 2:n) e[i] <- 0.6*e[i-1] + e[i]
y <- f + e
op <- par(mfrow=c(2,2))
b <- gamm(y~s(x,k=20),correlation=corAR1())
plot(b$gam);lines(x,f-mean(f),col=2)
b <- gamm(y~s(x,k=20))
plot(b$gam);lines(x,f-mean(f),col=2)
b <- gam(y~s(x,k=20))
plot(b);lines(x,f-mean(f),col=2)

## more complicated autocorrelation example - AR errors
## only within groups defined by `fac'
e <- rnorm(n,0,sig)
for (i in 2:n) e[i] <- 0.6*e[i-1]*(fac[i-1]==fac[i]) + e[i]
y <- f + e
b <- gamm(y~s(x,k=20),correlation=corAR1(form=~1|fac))
plot(b$gam);lines(x,f-mean(f),col=2)
par(op) 

## more complex situation with nested random effects and within
## group correlation 

set.seed(0)
n.g <- 10
n<-n.g*10*4
sig <- 2
## simulate smooth part
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f <- 2 * sin(pi * x0)
f <- f + exp(2 * x1) - 3.75887
f <- f+0.2*x2^11*(10*(1-x2))^6+10*(10*x2)^3*(1-x2)^10-1.396
## simulate nested random effects....
fa <- as.factor(rep(1:10,rep(4*n.g,10)))
ra <- rep(rnorm(10),rep(4*n.g,10))
fb <- as.factor(rep(rep(1:4,rep(n.g,4)),10))
rb <- rep(rnorm(4),rep(n.g,4))
for (i in 1:9) rb <- c(rb,rep(rnorm(4),rep(n.g,4)))
## simulate auto-correlated errors within groups
e<-array(0,0)
for (i in 1:40) {
  eg <- rnorm(n.g, 0, sig)
  for (j in 2:n.g) eg[j] <- eg[j-1]*0.6+ eg[j]
  e<-c(e,eg)
}
y <- f + ra + rb + e
dat<-data.frame(y=y,x0=x0,x1=x1,x2=x2,x3=x3,fa=fa,fb=fb)
## fit model .... 
b <- gamm(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
            s(x3,bs="cr"),data=dat,random=list(fa=~1,fb=~1),
          correlation=corAR1())
plot(b$gam,pages=1)

## and a "spatial" example
library(nlme);set.seed(1)
test1<-function(x,z,sx=0.3,sz=0.4)
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
                 0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n<-200
old.par<-par(mfrow=c(2,2))
x<-runif(n);z<-runif(n);
xs<-seq(0,1,length=30);zs<-seq(0,1,length=30)
pr<-data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth <- matrix(test1(pr$x,pr$z),30,30)
contour(xs,zs,truth)  # true function
f <- test1(x,z)  # true expectation of response
## Now simulate correlated errors...
cstr <- corGaus(.1,form = ~x+z)  
cstr <- Initialize(cstr,data.frame(x=x,z=z))
V <- corMatrix(cstr) # correlation matrix for data
Cv <- chol(V)
e <- t(Cv) %*% rnorm(n)*0.05 # correlated errors
## next add correlated simulated errors to expected values
y <- f + e ## ... to produce response
b<- gamm(y~s(x,z,k=50),correlation=corGaus(.1,form=~x+z))
plot(b$gam) # gamm fit accounting for correlation
# overfits when correlation ignored.....  
b1 <- gamm(y~s(x,z,k=50));plot(b1$gam) 
b2 <- gam(y~s(x,z,k=50));plot(b2)
par(old.par)