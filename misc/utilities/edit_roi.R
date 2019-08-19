library(dplyr)


dir_exp<-"C:/Users/atiro/Dropbox/temp/Shen_to_ICBM/roi_dictionary/source"

df_out<-data.frame(matrix(nrow=268,ncol=0))
df_out$number<-1:268
df_out$atlas<-"shen268"
df_out$id<-paste(df_out$atlas,sprintf("%05d",df_out$number),sep="_")

df_group<-read.csv(file.path(dir_exp,"shen268CommunityAffiliation.csv"))
df_roiconvert<-read.csv(file.path(dir_exp,"shen268CommunityNames.csv"))

df_group<-left_join(df_group,df_roiconvert,by="id_group")

df_out<-cbind(df_out,df_group)
df_out$label<-paste(as.character(df_out$number),df_out$name_group_long,sep=" ")
df_out$group_1<-df_out$group_2<-df_out$group_3<-df_out$name_group
df_out$label_short<-NA
df_out2<-df_out
df_out<-df_out[,c("id","atlas","number","label_short","label","group_1","group_2","group_3")]

write.csv(df_out,file.path(dir_exp,"roi_shen.csv"),row.names=F)

df_out2$label_xcp<-paste(df_out2$name_group_xcp,as.character(df_out2$number),sep="_")
df_out2<-df_out2[,c("number","label_xcp")]

write.csv(df_out2,file.path(dir_exp,"label_xcp_shen.csv"),row.names=F)

####

working_dir <- "G:/MRI/Statistics/Structural_FS"
setwd(working_dir)

a<-read.csv("ROI_CONN.csv")
a$short<-as.character(a$short)


e<-data.frame(matrix(nrow=nrow(a),ncol=1))
for (i in 1:nrow(a)){

  start <- regexpr("\\(", a[i,2])
  stop <- regexpr("\\)", a[i,2])
  if (start>=0){
    b<-a[i,2]
    c<-substr(b,start+1,stop-1)
    a[i,3]<-c
    d<-substr(b,1,start-2)
    a[i,2]<-d
  }
  
  a[i,1]<-sprintf("%05d",as.numeric(a[i,1]))
  e[i,1]<-sprintf("CNA_%05d",as.numeric(a[i,1]))
}
a<-cbind(e,a)
colnames(a)[1]<-"ID_long"





a<-read.csv("ROI_FreeSurfer.csv")

e<-data.frame(matrix(nrow=nrow(a),ncol=1))
for (i in 1:nrow(a)){
  
  a[i,1]<-sprintf("%05d",as.numeric(a[i,1]))
  e[i,1]<-sprintf("FS_%05d",as.numeric(a[i,1]))
}
a<-cbind(e,a)
colnames(a)[1]<-"ID_long"

write.csv(a,file.path(working_dir,"ROI_FreeSurfer2.csv"),row.names=F)



a<-read.csv("ROI_CONN_N.csv")
a<-a[-1]

e<-data.frame(matrix(nrow=nrow(a),ncol=1))
for (i in 1:nrow(a)){
  
  a[i,1]<-sprintf("%05d",as.numeric(a[i,1]))
  e[i,1]<-sprintf("CNN_%05d",as.numeric(a[i,1]))
}
a<-cbind(e,a)
colnames(a)[1]<-"ID_long"

write.csv(a,file.path(working_dir,"ROI_CONN_N2.csv"),row.names=F)




a<-read.csv("ROI_Power.csv")

e<-data.frame(matrix(nrow=nrow(a),ncol=1))
for (i in 1:nrow(a)){
  
  a[i,1]<-sprintf("%05d",as.numeric(a[i,1]))
  e[i,1]<-sprintf("PWR_%05d",as.numeric(a[i,1]))
}
a<-cbind(e,a)
colnames(a)[1]<-"ID_long"

write.csv(a,file.path(working_dir,"ROI_Power2.csv"),row.names=F)




a<-read.csv("ROI_FreeSurfer2.csv")

e<-sprintf("fs.cluster%03d",a$ID)

a<-cbind(e,a)

write.csv(a,file.path(working_dir,"ROI_FreeSurfer3.csv"),row.names=F)

####
working_dir <- "D:/MRI/Statistics/CommonData"
setwd(working_dir)
a<-read.csv("CSUB_Clinical_Data.csv")

write.csv(a,file.path(working_dir,"CSUB_Clinical_Data.csv"),row.names=F)

####
working_dir <- "D:/MRI/Statistics/CommonData"
setwd(working_dir)
a<-read.csv("ROI_All.csv")

write.csv(a,file.path(working_dir,"ROI_All.csv"),row.names=F)

####
