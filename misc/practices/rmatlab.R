
path<-"C:/Users/atiro/Dropbox/temp/prep_graphvar"

library(R.matlab)

df<-data.frame(col1=c(1,2),col2=c(3,4))

writeMat(file.path(path,"out.mat"),dataframe=df)
