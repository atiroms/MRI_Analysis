library(mgcv)
library(dplyr)
library(ggplot2)
library(itsadug)


df_str<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
df_clinical<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

#df_str<-read.csv("D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
#df_clinical<-read.csv("D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

df_str<-df_str[which(df_str['measure']=='volume'),]
df_str_w1<-df_str[which(df_str['wave']==1),]
df_str_w2<-df_str[which(df_str['wave']==2),]
df_str_w1<-left_join(df_str_w1,df_clinical[,c('ID_pnTTC','W1_Age_at_MRI','W1_Tanner_Max')],by='ID_pnTTC')
df_str_w2<-left_join(df_str_w2,df_clinical[,c('ID_pnTTC','W2_Age_at_MRI','W2_Tanner_Max')],by='ID_pnTTC')
colnames(df_str_w1)[which(colnames(df_str_w1)=='W1_Age_at_MRI')]<-'age'
colnames(df_str_w1)[which(colnames(df_str_w1)=='W1_Tanner_Max')]<-'tanner_max'
colnames(df_str_w2)[which(colnames(df_str_w2)=='W2_Age_at_MRI')]<-'age'
colnames(df_str_w2)[which(colnames(df_str_w2)=='W2_Tanner_Max')]<-'tanner_max'

df_str<-rbind(df_str_w1,df_str_w2)

df_str<-df_str[which(df_str['roi']=='dk_01001'),]
df_str<-df_str[which(!is.na(df_str['tanner_max'])),]
df_str$ID_pnTTC<-as.factor(df_str$ID_pnTTC)

model<-gam(value ~ s(age) + s(tanner_max,k=5) + s(ID_pnTTC,bs='re'),data=df_str)

summary(model)


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


ggplot(df_clinical_rbind) +
  #aes(x=age,y=tanner_max,color=wave) +
  aes(x=age,y=tanner_full,color=wave) +
  geom_point() +
  geom_path(aes(group=ID_pnTTC,color=NULL)) +
  ggtitle("Tanner stage vs Age")


mod_clinical<-gam(tanner_max ~ s(age) + s(ID_pnTTC,bs='re'),data=df_clinical_rbind)
summary(mod_clinical)
plot(mod_clinical,select=2)
plot_smooth(mod_clinical,view='age')
  