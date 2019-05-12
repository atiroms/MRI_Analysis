
path_exp<-"P:/MRI/pnTTC/Preproc"
list_src<-list(list("dir"="51_c1_1_ts","ses"=1),
               list("dir"="51_c1_2_ts","ses"=1),
               list("dir"="52_c2_1_ts","ses"=2),
               list("dir"="52_c2_2_ts","ses"=2))
dir_dst<-"53_ts_acompcor"

list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
ses<-1


for (atlas in list_atlas){
  print(paste("Calculating atlas: ",atlas,sep=""))
  df_out<-data.frame()
  for (src in list_src){
    dir_src<-src$dir
    ses<-src$ses
    path_file_in<-file.path(path_exp,dir_src,"output",paste("ts_",atlas,".csv",sep=""))
    df_out_add<-read.csv(path_file_in)
    df_out_add<-cbind(ses=ses,df_out_add)
    df_out<-rbind(df_out,df_out_add)
  }
  file_out<-paste("atl-",atlas,"_ts.csv",sep="")
  path_file_out<-file.path(path_exp,dir_dst,"output",file_out)
  write.csv(df_out,path_file_out,row.names=F)
}
print("Finish")