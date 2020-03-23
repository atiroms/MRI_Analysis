# confer below
# https://stackoverflow.com/questions/38516139/gam-r-variance-explained-by-variable


df<-data.frame(a=seq(100))
df$b<-sin(df$a)
df$c<-runif(100)
df$y<-3*df$a+2*df$b+df$c

library(mgcv)

mod<-gam(y~a+b,data=df,method="REML")
mod_a<-gam(y~a,data=df,method="REML")
mod_b<-gam(y~b,data=df,method="REML")
mod_0<-gam(y~1,data=df)

ve_a<-(deviance(mod_b)-deviance(mod))/deviance(mod_0)
ve_b<-(deviance(mod_a)-deviance(mod))/deviance(mod_0)
