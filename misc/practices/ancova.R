# mod6 is the right way

library(multcomp)

data("mtcars")
head(mtcars)

mtcars$am<-as.factor(mtcars$am)
mtcars$vs<-as.factor(mtcars$vs)

mod1<-lm(mpg~am*vs+hp,data=mtcars)
summary(mod1)

mod2<-lm(mpg~am:vs+hp,data=mtcars)
summary(mod2)

mod3<-lm(mpg~am:vs+hp-1,data=mtcars)
summary(mod3)

mod4<-lm(mpg~am:vs+hp+1,data=mtcars)
summary(mod4)

mod5<-aov(mpg~am:vs+hp-1,data=mtcars)
summary(mod5)

contrast1<-matrix(0L,nrow=3,ncol=5)
contrast1[1,2]<-contrast1[2,3]<-contrast1[3,4]<-1
contrast1[1,3]<-contrast1[2,4]<-contrast1[3,5]<--1
contrast1

test1<-glht(mod3, linfct = contrast1)
summary(test1)

mod6<-aov(mpg~am:vs+hp,data=mtcars)
summary(mod6)

mod7<-lm(mpg~am:vs+hp,data=mtcars)
summary(mod7)

sum5<-summary(mod5)[[1]]
sum6<-summary(mod6)[[1]]
