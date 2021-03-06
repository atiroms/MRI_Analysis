# mod10, posthoc10_2 are right

library(multcomp)

data("mtcars")
head(mtcars)
summary(mtcars)

mtcars$am<-as.factor(mtcars$am)
mtcars$vs<-as.factor(mtcars$vs)
mtcars$gear<-as.factor(mtcars$gear)
mtcars$am_vs<-as.factor(paste(as.character(mtcars$am),as.character(mtcars$vs),sep="_"))

summary(mtcars)

mod6<-aov(mpg~am:vs+hp,data=mtcars)
summary(mod6)
posthoc6<-TukeyHSD(mod6,which="am:vs")[[1]]
posthoc6_2<-glht(mod10, linfct = mcp("am:vs" = 'Tukey'))

mod10<-aov(mpg~am_vs+hp,data=mtcars)
summary(mod10)
posthoc10<-TukeyHSD(mod10,which="am_vs")[[1]] # TukeyHSD() does NOT support ANCOVA
posthoc10_2<-glht(mod10, linfct = mcp(am_vs = 'Tukey')) # Right way
posthoc10_2<-summary(posthoc10_2)$test

####
mean_0_0<-mean(mtcars[mtcars$am==0 & mtcars$vs==0,"mpg"])
mean_1_1<-mean(mtcars[mtcars$am==1 & mtcars$vs==1,"mpg"])
diff_mean<-mean_1_1-mean_0_0        # this mean was equal to mod10 mean

####

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

mod7<-lm(mpg~am:vs+hp,data=mtcars)
summary(mod7)

sum5<-summary(mod5)[[1]]
sum6<-summary(mod6)[[1]]

mod8<-aov(mpg~gear+hp,data=mtcars)
sum8<-summary(mod8)[[1]]
posthoc8<-TukeyHSD(mod8,which='gear')[[1]]

mod9<-aov(mpg~gear,data=mtcars)
posthoc9<-TukeyHSD(mod9)
