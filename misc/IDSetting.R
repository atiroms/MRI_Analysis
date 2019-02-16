##

working_dir <- "D:/atiroms/Dropbox/MRI/Statistics/CommonData"
input_filename<-"rsfMRI.csv"
output_filename<-"rsfMRI2.csv"

input<-read.csv(file.path(working_dir,input_filename))

output<-data.frame(matrix(ncol=ncol(input),nrow=351))
colnames(output)<-colnames(input)
output$ID_pnTTC<-1:nrow(output)

for (i in 1:nrow(input)){
  output[input[i,"ID_pnTTC"],]<-input[i,]
}

write.csv(output, file.path(working_dir,output_filename),row.names=T)

##

working_dir <- "D:/atiroms/Dropbox/MRI/Statistics/CommonData"
input_filename<-"CSUB2_02.csv"
output_filename<-"CSUB2_03.csv"

input<-read.csv(file.path(working_dir,input_filename))

input$X<-as.character(input$X)
input$pn.TTC.01<-as.character(input$pn.TTC.01)
input$T1<-as.character(input$T1)
input$rsfMRI<-as.character(input$rsfMRI)
input$MRS<-as.character(input$MRS)
input$QC_T1.RAW_1<-as.numeric(as.character(input$QC_T1.RAW_1))
input$QC_ANOMALY_1<-as.numeric(as.character(input$QC_ANOMALY_1))
input$QC備考<-as.character(input$QC備考)

input$ID<-substr(input$X,6,10)
input$ID<-as.numeric(input$ID)

output<-data.frame(matrix(ncol=15,nrow=351))
colnames(output)<-colnames(input)
for (i in 1:nrow(input)){
  output[input[i,"ID"],]<-input[i,]
}
write.csv(output, file.path(working_dir,output_filename),row.names=T)


#

working_dir <- "G:/MRI/Info"
input_filename<-"ID5.csv"
output_filename<-"ID6.csv"


input_filename2<-"CONN_voxel_QC.csv"
setwd(working_dir)

input<-read.csv(file.path(working_dir,input_filename))
input2<-read.csv(file.path(working_dir,input_filename2))

a<-input


counter<-0
for (i in 1:dim(input)[1]){
  if (!is.na(input[i,4])){
    if (input[i,4]==1){
      counter<-counter+1
      input[i,3]<-counter
    }
  }
}



for (i in 1:dim(a)[1]){
  if (is.na(a[i,4])){
    a[i,5]<-NA
  }else{
    a[i,5]<-1
  }
}

for (i in 1:dim(a)[1]){
  if (is.na(a[i,3])){
    a[i,3]<-0
  }else{
    a[i,3]<-1
  }
}

for (i in 1:dim(a)[1]){
  if (is.na(a[i,6])){
    a[i,7]<-NA
  }else if(a[i,6]==1){
    a[i,7]<-a[i,6]
  }
}

for (i in 1:dim(a)[1]){
  if (!is.na(a[i,4])){
    if (a[i,4]==1){
      a[i,7]
    }
  }
}



output<-data.frame(matrix(ncol=2,nrow=dim(input)[1]))
for (i in input2$ID_T1QC_rsfMRIexist){
  input_index<-which(input2$ID_T1QC_rsfMRIexist==i)
  input_value20<-input2$QC20[input_index]
  input_value10<-input2$QC10[input_index]
  target_index<-which(input$ID_T1QC_rsfMRIexist==i)
  output[target_index,1]<-input_value20
  output[target_index,2]<-input_value10
}
colnames(output)<-c("QC20","QC10")


a<-cbind(input,output)
colnames(a)[5]<-"T1_exist"
b<-data.frame(a$ID_pnTTC,a$ID_TTC,a$T1_exist,a$T1QC,a$ID_T1QC)

#####

working_dir <- "C:/Users/atiro/Dropbox/MRI/Statistics/CommonData"
input_file<-"CSUB.csv"
output_file<-"CSUB.csv"

output<-read.csv(file.path(working_dir,input_file))

output$ID_W1_T1QC[which(output$W1_T1QC==1)]<-seq(length(which(output$W1_T1QC==1)))
output$ID_W1_T1QC[which(output$W1_T1QC==0 | is.na(output$W1_T1QC))]<-NA

output$ID_W1_T1QC_rsfMRIexist[which(output$W1_T1QC_rsfMRIexist==1)]<-seq(length(which(output$W1_T1QC_rsfMRIexist==1)))
output$ID_W1_T1QC_rsfMRIexist[which(output$W1_T1QC_rsfMRIexist==0 | is.na(output$W1_T1QC_rsfMRIexist))]<-NA

output$ID_W1_T1QC_rsfMRIexist_CONNvoxelQC20[which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC20==1)]<-seq(length(which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC20==1)))
output$ID_W1_T1QC_rsfMRIexist_CONNvoxelQC20[which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC20==0 | is.na(output$W1_T1QC_rsfMRIexist_CONNvoxelQC20))]<-NA

output$ID_W1_T1QC_rsfMRIexist_CONNvoxelQC10[which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC10==1)]<-seq(length(which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC10==1)))
output$ID_W1_T1QC_rsfMRIexist_CONNvoxelQC10[which(output$W1_T1QC_rsfMRIexist_CONNvoxelQC10==0 | is.na(output$W1_T1QC_rsfMRIexist_CONNvoxelQC10))]<-NA



library(plyr)
output$Sex<-mapvalues(output$Sex, from=c(1,2), to=c("male", "female"))
write.csv(output, file.path(working_dir,output_file),row.names=F)

