library(R.matlab)
library(ggplot2)
#library(abind)

working_dir<-"D:/atiroms/MRI/iTTC_rsfMRI_C/CONN/Timeseries/V2V_04"
setwd(working_dir)
input_filename<-"PCA.Timeseries.mat"
matlabdata<-readMat(input_filename)
n_participant<-length(matlabdata$data)



extractnesteddata<-function(input,i){
  output<-input[[i]]
  for (j in 1:2){
    output<-output[[1]]
  }
  timepoint<-1:dim(output[[1]])[1]
  output<-as.data.frame(output)
  output<-cbind(participant=i,timepoint=timepoint,output)
  return(output)
}


timeseries<-extractnesteddata(matlabdata$data, 1)
#  output<-array(output, dim=c(1, nrow(output), ncol(output)))
for (i in 2:length(matlabdata$data)){
  idata<-extractnesteddata(matlabdata$data, i)
  timeseries<-rbind(timeseries,idata)
#  output<-abind(output, array(idata, dim=c(1,nrow(idata),ncol(idata))), along=1)
}


ggplot(timeseries, aes(X1, X2,group=participant,color=participant)) + geom_path()
ggplot(timeseries, aes(X2, X3,group=participant,color=participant)) + geom_path()
