path_exp<-"C:/Users/atiro/Dropbox/MRI/pnTTC/Info/TTC_pnTTC_ID_data"

ID_integrate<-read.csv(file.path(path_exp,"ID_integrate.csv"))

ID_dict<-read.csv(file.path(path_exp,"ID_dict.csv"))

for (i in seq(nrow(ID_dict))){
  ID_integrate[which(paste(ID_integrate$ID_string,"C-02",sep="")==ID_dict[i,"ID_string"]),"ID_TTC_2"]<-ID_dict[i,"ID_TTC_2"]  
}

for (i in seq(nrow(ID_integrate))){
  if(!is.na(ID_integrate$ID_TTC_1[i])){
    ID_integrate$ID_TTC_integrate[i]<-ID_integrate$ID_TTC_1[i]
  }else if (!is.na(ID_integrate$ID_TTC_2[i])){
    ID_integrate$ID_TTC_integrate[i]<-ID_integrate$ID_TTC_2[i]
  }else{
    ID_integrate$ID_TTC_integrate[i]<-NA
  }
}

write.csv(ID_integrate,file.path(path_exp,"ID_integrate2.csv"),row.names=F)


ID_W2_exist<-read.csv(file.path(path_exp,"ID_W2_exist.csv"))

for (i in seq(nrow(ID_W2_exist))){
  ID_W2_exist$ID_TTC[i]<-ID_integrate$ID_TTC_integrate[which(paste(ID_integrate$ID_string,"C-02",sep="")==ID_W2_exist$ID_string[i])]
}

write.csv(ID_W2_exist,file.path(path_exp,"ID_W2_exist2.csv"),row.names=F)

nrow(ID_W2_exist)

length(unique(ID_W2_exist$ID_TTC))
unique(ID_W2_exist$ID_TTC)
