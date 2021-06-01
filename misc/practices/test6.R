path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full<-NULL

dir_in<-"421_fc_aroma"
dir_out<-"427_fc_gamm_aroma_test17"

library(data.table)
library(mgcv)
source(file.path(getwd(),"util/function.R"))

paths_<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)

####

memory.limit(1000000)

atlas<-"aal116"

print(paste("Preparing FC data: ",atlas,sep=""))
df_fc<-as.data.frame(fread(file.path(paths_$input,"output",
                                     paste("atl-",atlas,"_fc.csv",sep=""))))
df_fc<-df_fc[df_fc$ses!="2-1",]

####

list_covar<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                 "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                 "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex"))
list_tanner_<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                  "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                  "label"="Tanner stage (gonadal)"),
                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                  "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                  "label"="Tanner stage (adrenal)"))
idx_tanner<-names(list_tanner_)[1]
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]

list_wave_<-c(1,2)
subset_subj_<- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                            list("key"="W1_rsfMRIexist","condition"="==1"),
                            list("key"="W1_Censor","condition"="<126")),
                   "2"=list(list("key"="W2_T1QC","condition"="==1"),
                            list("key"="W2_rsfMRIexist","condition"="==1"),
                            list("key"="W2_Censor","condition"="<126")))


####

# Prepare clinical data
data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T,
                                   prefix=paste("var-",idx_tanner,sep=""),print_terminal=F)
df_clin<-data_clin$df_clin

# Join FC and clinical data
df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                -which(colnames(df_fc)=="p"))]
df_clin$wave<-as.character(df_clin$wave)
df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex')){
  if (key %in% colnames(df_join)){
    df_join[,key]<-as.factor(df_join[,key])
  }
}
df_join$value<-as.numeric.factor(df_join$value)

####

df_join_subset<-df_join

list_node<-sort(unique(c(df_join$from,df_join$to)))
list_node<-list_node[1:10]
df_join_subset<-df_join[(df_join$from %in% list_node) & (df_join$to %in% list_node),]

df_join_subset$edge<-paste(df_join_subset$from,df_join_subset$to,sep='_')
df_join_subset$edge<-factor(df_join_subset$edge,ordered=F)
df_join_subset$from<-df_join_subset$to<-NULL
df_join_subset$value<-as.numeric(df_join_subset$value)

####

#mod1<-gam(value ~ age + s(tanner,k=3) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
#summary(mod1)
#plot(mod1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

mod2<-gam(value ~ age + te(edge,tanner,bs=c('re','tp')) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod2)
plot(mod2, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
vis.gam(mod2, view=c("edge","tanner"), color="cm", theta=0,phi=0)
vis.gam(mod2, view=c("edge","tanner"), color="cm", theta=45)

####
# BAM implementation of above
mod3<-bam(value ~ age + te(edge,tanner,bs=c('re','tp')) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod3)
plot(mod3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

mod4<-gam(value ~ age + s(tanner,id=1,by=edge,k=3) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod4)
vis.gam(mod4, view=c("edge","tanner"), color="cm", theta=45)

####

mod5<-gam(value ~ age + s(tanner,id=1,by=edge,k=5) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod5)
vis.gam(mod5, view=c("edge","tanner"), color="cm", theta=45)

####

mod6<-gam(value ~ age + s(tanner,bs='cr',id=1,by=edge,k=5) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod6)
vis.gam(mod6, view=c("edge","tanner"), color="cm", theta=45)

####

mod7<-gam(value ~ age + s(tanner,bs='cr',id=1,by=edge,k=5) + s(ID_pnTTC,bs='re'),data=df_join_subset,method="REML")
summary(mod7)
vis.gam(mod7, view=c("edge","tanner"), color="cm", theta=45)