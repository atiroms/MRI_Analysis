gamm_fc_list_tanner<-gamm_fc_list_tanner[1]
list_atlas<-list_atlas[1]

####

paths_=paths
subset_subj_=gamm_fc_subset_subj
list_wave_=list_wave
list_atlas_=list_atlas
key_group_='group_3'
list_covar_tanner_=gamm_fc_list_covar_tanner
list_tanner_=gamm_fc_list_tanner
list_mod_tanner_=gamm_fc_list_mod_tanner
list_plot_tanner_=gamm_fc_list_plot_tanner
list_covar_hormone_=gamm_fc_list_covar_hormone
list_hormone_=gamm_fc_list_hormone
list_mod_hormone_=gamm_fc_list_mod_hormone
list_plot_hormone_=gamm_fc_list_plot_hormone

####

print("Starting gamm_fc_multi().")
nullobj<-func_createdirs(paths_,str_proc="gamm_fc_multi()",copy_log=T)
dict_roi <- func_dict_roi(paths_)

####

atlas<-list_atlas[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
df_fc<-as.data.frame(fread(file.path(paths_$input,"output",
                                     paste("atl-",atlas,"_fc.csv",sep=""))))

# Prepare dataframe of ROIs
list_roi<-sort(unique(c(as.character(df_fc$from),as.character(df_fc$to))))
df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"

# prepare dataframe of group-wise FC averages
list_group<-unique(as.character(df_roi$group))
df_fc_temp<-df_fc
df_fc_temp$z_r[which(is.nan(df_fc_temp$z_r))]<-0
df_fc_temp<-inner_join(df_fc_temp,df_roi[,c("id","group")],by=c("from"="id"))
colnames(df_fc_temp)[colnames(df_fc_temp)=="group"]<-"from_group"
df_fc_temp<-inner_join(df_fc_temp,df_roi[,c("id","group")],by=c("to"="id"))
colnames(df_fc_temp)[colnames(df_fc_temp)=="group"]<-"to_group"
df_subj<-NULL
list_subj<-sort(unique(df_fc$ID_pnTTC))
for (id_subj in list_subj){
  list_ses<-sort(unique(df_fc[df_fc$ID_pnTTC==id_subj,"ses"]))
  list_ses<-list_ses[list_ses!="2-1"]
  df_subj<-rbind(df_subj,data.frame(ID_pnTTC=id_subj,ses=list_ses))
}

df_fc_grp<-NULL
for (idx_subj_ses in seq(dim(df_subj)[1])){
  #print(paste(df_subj[idx_subj_ses,"ID_pnTTC"],df_subj[idx_subj_ses,"ses"]))
  df_fc_subset1<-df_fc_temp[df_fc_temp$ID_pnTTC==df_subj[idx_subj_ses,"ID_pnTTC"]
                            & df_fc_temp$ses==df_subj[idx_subj_ses,"ses"],]
  for (idx_grp1 in seq(length(list_group))){
    for (idx_grp2 in seq(idx_grp1,length(list_group))){
      df_fc_subset2<-rbind(df_fc_subset1[df_fc_subset1$from_group==list_group[idx_grp1]
                                         & df_fc_subset1$to_group==list_group[idx_grp2],],
                          df_fc_subset1[df_fc_subset1$from_group==list_group[idx_grp2]
                                        & df_fc_subset1$to_group==list_group[idx_grp1],])
      df_fc_grp<-rbind(df_fc_grp,
                       cbind(ID_pnTTC=df_subj[idx_subj_ses,"ID_pnTTC"],ses=df_subj[idx_subj_ses,"ses"],
                             from=list_group[idx_grp1],to=list_group[idx_grp2],
                             z_r=mean(df_fc_subset2$z_r)))
    }
  }
}

####

idx_tanner<-names(list_tanner_)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]

####

list_mod<-list_mod_tanner_
list_plot<-list_plot_tanner_
idx_var<-idx_tanner
list_type_p_=list_type_p
thr_p_=thr_p

####

data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T,
                                   prefix=paste("var-",idx_var,sep=""),print_terminal=F)
df_clin<-data_clin$df_clin

# Join FC and clinical data
df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                -which(colnames(df_fc)=="p"))]
df_clin$wave<-as.character(df_clin$wave)
df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex')){
  if (key %in% colnames(df_join)){
    df_join[,key]<-as.factor(df_join[,key])
  }
}

# Calculate ROI-wise GAMM of FC
data_gamm<-iterate_gamm(df_join,df_roi,list_mod,calc_parallel=F,calc_identical=F)
df_plot<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod,
                    list_plot,calc_seed_level=F)

# Join group-wise FC and clinical data
colnames(df_fc_grp)[colnames(df_fc_grp)=="z_r"]<-"value"
colnames(df_fc_grp)[colnames(df_fc_grp)=="ses"]<-"wave"
df_join_grp<-inner_join(df_fc_grp,df_clin,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex')){
  if (key %in% colnames(df_join_grp)){
    df_join_grp[,key]<-as.factor(df_join_grp[,key])
  }
}

# Calculate Group-wise GAMM of FC
list_group<-unique(as.character(df_roi$group))
df_grp<-data.frame(id=list_group,label=str_to_title(gsub("_"," ",as.character(list_group))))
data_gamm_grp<-iterate_gamm(df_join_grp,df_grp,list_mod,calc_parallel=F,calc_identical=T)
df_plot_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_grp,list_mod,
                        list_plot,calc_seed_level=F)

# Save results
write.csv(data_gamm$df_out_gamm,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
write.csv(data_gamm$df_out_aic,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
write.csv(df_plot,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_out_gamm,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_out_aic,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
write.csv(df_plot_grp,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot_grp.csv",sep="")),row.names = F)

# Graphical output of ROI- and group-wise GAMM of FC
plot_gam_fc(paths_,df_gam=df_plot,df_gam_grp_sign=df_plot_grp,df_gam_grp_abs=NULL,atlas,
            list_mod,list_plot,list_type_p=list_type_p_,thr_p=thr_p_,waves=NULL,idx_var)
#plot_gam_fc(df_plot,df_roi,analysis="roi",atlas,list_mod,list_plot,
#            list_type_p_,thr_p,paths_,suffix_=paste("var-",idx_var,sep=""))