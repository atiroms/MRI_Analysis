library(mgcv)
library(dplyr)
library(ggplot2)
library(itsadug)
library(ggrepel)
#library(voxel)
library(purrr)


#df_str<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
#df_clinical<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

df_str<-read.csv("D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
df_clinical<-read.csv("D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

#df_str<-read.csv("C:/Users/NICT_WS/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
#df_clinical<-read.csv("C:/Users/NICT_WS/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

df_str<-df_str[which(df_str['measure']=='volume'),]
df_str_w1<-df_str[which(df_str['wave']==1),]
df_str_w2<-df_str[which(df_str['wave']==2),]
df_str_w1<-left_join(df_str_w1,df_clinical[,c('ID_pnTTC','Sex','W1_Age_at_MRI','W1_Tanner_Max')],by='ID_pnTTC')
df_str_w2<-left_join(df_str_w2,df_clinical[,c('ID_pnTTC','Sex','W2_Age_at_MRI','W2_Tanner_Max')],by='ID_pnTTC')
colnames(df_str_w1)[which(colnames(df_str_w1)=='Sex')]<-'sex'
colnames(df_str_w1)[which(colnames(df_str_w1)=='W1_Age_at_MRI')]<-'age'
colnames(df_str_w1)[which(colnames(df_str_w1)=='W1_Tanner_Max')]<-'tanner'
colnames(df_str_w2)[which(colnames(df_str_w2)=='Sex')]<-'sex'
colnames(df_str_w2)[which(colnames(df_str_w2)=='W2_Age_at_MRI')]<-'age'
colnames(df_str_w2)[which(colnames(df_str_w2)=='W2_Tanner_Max')]<-'tanner'

df_str<-rbind(df_str_w1,df_str_w2)

df_str<-df_str[which(df_str['roi']=='dk_01001'),]
df_str<-df_str[which(!is.na(df_str['tanner'])),]
df_str$ID_pnTTC<-as.factor(df_str$ID_pnTTC)
df_str$sex<-as.factor(df_str$sex)

for (i in seq(nrow(df_str))){
  if (df_str[i,'sex']==1){
    df_str[i,'tanner_male']<-df_str[i,'tanner']
    df_str[i,'male']<-1
    df_str[i,'female']<-0
  }else if (df_str[i,'sex']==2){
    df_str[i,'tanner_female']<-df_str[i,'tanner']
    df_str[i,'male']<-0
    df_str[i,'female']<-1
  }
}

#model<-gam(value ~ s(age) + s(tanner,k=5) + s(ID_pnTTC,bs='re'),data=df_str)
#model<-gam(value ~ s(age) + s(tanner) + s(ID_pnTTC,bs='re'),data=df_str)
#model<-gam(value ~ age + tanner + s(ID_pnTTC,bs='re'),data=df_str)
#model2<-gam(value ~ age*tanner + s(ID_pnTTC,bs='re'),data=df_str)

#model3<-gam(value ~ age + tanner_male + tanner_female + age:tanner_male + age:tanner_female + s(ID_pnTTC,bs='re'),data=df_str,na.action =na.pass)
#model3<-gam(value ~ age + tanner_male + tanner_female  + s(ID_pnTTC,bs='re'),data=df_str,na.action=na.exclude)
#model3<-gam(value ~ s(age) + s(tanner,k=5,by=sex) + s(ID_pnTTC,bs='re'),data=df_str)
#summary(model3)

#model4<-gam(value ~ s(age) + s(tanner,k=5):male + s(tanner,k=5):female + s(ID_pnTTC,bs='re'),data=df_str)
#summary(model4)

#model5<-gam(value ~ s(age) + tanner:male + tanner:female + s(ID_pnTTC,bs='re'),data=df_str)
model6<-gam(value ~ age + tanner:sex + s(ID_pnTTC,bs='re'),data=df_str)
model7<-gam(value ~ age + tanner:sex + age:tanner:sex + s(ID_pnTTC,bs='re'),data=df_str)

model8<-gam(value ~ age + tanner:male + tanner:female + s(ID_pnTTC,bs='re'),data=df_str)
model9<-gam(value ~ age + tanner:male + age:tanner:male + tanner:female + age:tanner:female + s(ID_pnTTC,bs='re'),data=df_str)

#model4<-gam(value ~ age + tanner + s(ID_pnTTC,bs='re'),data=df_str)
#model5<-gam(value ~ age + s(tanner,k=2) + s(ID_pnTTC,bs='re'),data=df_str)

anova.gam(model6,model7,test="F")
anova(model6,model7,test="F")
AIC(model,model2)
step(model2)

#str_formula<-"value ~ s(age) + s(tanner_max,k=5) + s(ID_pnTTC,bs='re')"
#formula<-as.formula(str_formula)
#model<-gam(formula,data=df_str)

#formula_lm<-as.formula(str_formula_lm)
summary(model)
#summary.gam(model)$s.table
gam.check(model)


## plotting results
temp.data<-model$model
#class(temp.data)
smooth.cov<-'tanner_max'
plot.df <- data.frame(x = seq(min(temp.data[smooth.cov]),
                             max(temp.data[smooth.cov]),
                             length.out=200))
names(plot.df) <- smooth.cov
for (i in names(temp.data)[-1]) {
  if (i != smooth.cov) {
    if (any(class(temp.data[i][,1])[1] == c("numeric", "integer","boolean"))) {
      plot.df[, dim(plot.df)[2] + 1] <- mean(temp.data[i][,1])
      names(plot.df)[dim(plot.df)[2]] <- i
    }
    else if (any(class(temp.data[i][,1])[1] == c("character", "factor","ordered"))) {
      plot.df[, dim(plot.df)[2] + 1] <- temp.data[i][1,1]
      names(plot.df)[dim(plot.df)[2]] <- i
    }
  }
}

for (i in 1:dim(plot.df)[2]) {
  if (class(plot.df[,i])[1] == "ordered" |  class(plot.df[,i])[1] == "factor") {
    warning("There are one or more factors in the model fit, please consider plotting by group since plot might be unprecise")
  }
}
plot.df = cbind(plot.df, as.data.frame(predict.gam(model, plot.df, se.fit = TRUE)))

plot <- (ggplot(data=plot.df, aes(x=plot.df[,1]))
         + geom_line(aes(y=fit), size=1)
         + geom_ribbon(data=plot.df, aes(ymax = fit+1.96*se.fit, ymin = fit-1.96*se.fit, linetype=NA), alpha = .2)
         + geom_point(data = temp.data, aes(x=as_vector(temp.data[smooth.cov]), y=temp.data[,1]))
         + ggtitle("GAMM model")
         + ylab("Structural measure")
         + xlab("Tanner stage")
         + theme_light())



## another way of plotting
df_str<-df_str
df_str$resid<-resid(model)
df_str$predict<-predict(model)
plot <- ggplot(df_str, aes(x = tanner_max, y = value)) + 
  geom_point(shape = 19, stroke = 1) +
  xlab("tanner_max") + ylab("value") +
  theme_bw()

plot<-plot+geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = df_str, size = 1, se = FALSE)


#### longitudinal clinical data analysis ####
df_clinical_w1<-df_clinical[,c('ID_pnTTC','Sex','W1_Age_at_MRI','W1_Tanner_Max','W1_Tanner_Full')]
colnames(df_clinical_w1)[which(colnames(df_clinical_w1)=='Sex')]<-'sex'
colnames(df_clinical_w1)[which(colnames(df_clinical_w1)=='W1_Age_at_MRI')]<-'age'
colnames(df_clinical_w1)[which(colnames(df_clinical_w1)=='W1_Tanner_Max')]<-'tanner_max'
colnames(df_clinical_w1)[which(colnames(df_clinical_w1)=='W1_Tanner_Full')]<-'tanner_full'
df_clinical_w1$wave<-1
df_clinical_w2<-df_clinical[,c('ID_pnTTC','Sex','W2_Age_at_MRI','W2_Tanner_Max','W2_Tanner_Full')]
colnames(df_clinical_w2)[which(colnames(df_clinical_w2)=='Sex')]<-'sex'
colnames(df_clinical_w2)[which(colnames(df_clinical_w2)=='W2_Age_at_MRI')]<-'age'
colnames(df_clinical_w2)[which(colnames(df_clinical_w2)=='W2_Tanner_Max')]<-'tanner_max'
colnames(df_clinical_w2)[which(colnames(df_clinical_w2)=='W2_Tanner_Full')]<-'tanner_full'
df_clinical_w2$wave<-2
df_clinical_rbind<-rbind(df_clinical_w1,df_clinical_w2)
df_clinical_rbind$ID_pnTTC<-as.factor(df_clinical_rbind$ID_pnTTC)
df_clinical_rbind$wave<-as.factor(df_clinical_rbind$wave)
df_clinical_rbind$sex<-as.factor(df_clinical_rbind$sex)


ggplot(df_clinical_rbind) +
  #aes(x=age,y=tanner_max,color=wave, label=ID_pnTTC) +
  aes(x=age,y=tanner_full,color=sex,shape=wave, label=ID_pnTTC) +
  scale_colour_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral")) +
  scale_shape_manual(name=NULL,labels=c("1st wave","2nd wave"),values=c(3,4)) +
  geom_point(size=4) +
  geom_text_repel(size=2) +
  geom_path(aes(group=ID_pnTTC),color="black",size=0.5,alpha=0.5) +
  ggtitle("Tanner stage vs Age") +
  xlab("Age (day)") +
  ylab("Tanner Stage") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),legend.justification=c(0,1), legend.position=c(0.05,0.95),legend.direction="horizontal",panel.grid.minor=element_blank())
  


mod_clinical<-gam(tanner_max ~ s(age) + s(ID_pnTTC,bs='re'),data=df_clinical_rbind)
summary(mod_clinical)
plot(mod_clinical,select=2)
plot_smooth(mod_clinical,view='age')
  