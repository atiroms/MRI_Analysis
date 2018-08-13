#

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
