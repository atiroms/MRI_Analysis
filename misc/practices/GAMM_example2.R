library(mgcv)
library(dplyr)


df_str<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract/output/fs_measure.csv")
df_clinical<-read.csv("C:/Users/atiro/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB.csv")

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
