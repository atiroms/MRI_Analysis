path_exp<-"/home/atiroms/Dropbox/MRI/pnTTC/Info/QC_New"

df_in<-read.csv(file.path(path_exp,'w2_source.csv'))

df_in$ID_pnTTC<-gsub("CSUB-","",df_in$ID_str)
df_in$ID_pnTTC<-gsub("C-02","",df_in$ID_pnTTC)
df_in$ID_pnTTC<-as.integer(df_in$ID_pnTTC)

df_out<-data.frame(matrix(nrow=max(df_in$ID_pnTTC),ncol=ncol(df_in)))
colnames(df_out)<-colnames(df_in)
df_out$ID_pnTTC<-seq(max(df_in$ID_pnTTC))

for (id_sub in df_in$ID_pnTTC){
  df_out[which(df_out$ID_pnTTC==id_sub),]<-df_in[which(df_in$ID_pnTTC==id_sub),]
}

df_out<-cbind(df_out[,"ID_pnTTC"],df_out[,which(colnames(df_out)!="ID_pnTTC")])
colnames(df_out)[1]<-"ID_pnTTC"

write.csv(df_out,file.path(path_exp,"w2_spaced.csv"),row.names=F)
