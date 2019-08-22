
#path<-"C:/Users/atiro/Dropbox/temp/prep_graphvar"
path<-"C:/Users/NICT_WS/Dropbox/temp/prep_graphvar"

library(R.matlab)

#df<-data.frame(col1=c(1,2),col2=c(3,4))

mat_out<-matrix(1:9,nrow=3,ncol=3)

writeMat(file.path(path,"out.mat"),matrix=mat_out)
