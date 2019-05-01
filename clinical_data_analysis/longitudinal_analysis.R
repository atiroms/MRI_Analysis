#**************************************************
# Description =====================================
#**************************************************

# R script to analyze longitudinal clinical data.


#**************************************************
# Parameters ======================================
#**************************************************
df_clinical<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(dplyr)
library(ggplot2)
library(itsadug)
library(ggrepel)


#**************************************************
# Dataframe Setup =================================
#**************************************************
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


#**************************************************
# Longitudinal change Graph =======================
#**************************************************
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


#**************************************************
# GAMM ============================================
#**************************************************
mod_clinical<-gam(tanner_max ~ s(age) + s(ID_pnTTC,bs='re'),data=df_clinical_rbind)
summary(mod_clinical)
plot(mod_clinical,select=2)
plot_smooth(mod_clinical,view='age')