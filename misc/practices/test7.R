# test standardization of FC

df_fc<-data_fc$df_fc

df_mean_sd<-data.frame()
list_wave<-sort(unique(df_fc$ses))
for (wave in list_wave){
  list_id_subj<-sort(unique(df_fc[df_fc$ses==wave,"ID_pnTTC"]))
  for (id_subj in list_id_subj){
    df_fc_subset<-df_fc[df_fc$ses==wave & df_fc$ID_pnTTC==id_subj,]
    mean_fc<-mean(df_fc_subset$z_r)
    sd_fc<-sd(df_fc_subset$z_r)
    df_mean_sd<-rbind(df_mean_sd,
                      data.frame(ses=wave,ID_pnTTC=id_subj,mean=mean_fc,sd=sd_fc))
  }
}


df_fc_std<-data.frame()
list_wave<-sort(unique(df_fc$ses))
for (wave in list_wave){
  list_id_subj<-sort(unique(df_fc[df_fc$ses==wave,"ID_pnTTC"]))
  for (id_subj in list_id_subj){
    df_fc_subset<-df_fc[df_fc$ses==wave & df_fc$ID_pnTTC==id_subj,]
    mean_fc<-mean(df_fc_subset$z_r)
    sd_fc<-sd(df_fc_subset$z_r)
    df_fc_subset$z_r<-(df_fc_subset$z_r-mean_fc)/sd_fc
    df_fc_std<-rbind(df_fc_std,df_fc_subset[,c("ses","ID_pnTTC","from","to","z_r")])
  }
}

