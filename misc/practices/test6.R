paths_=paths
list_waves_=ca_fc_list_waves
subset_subj_=ca_fc_subset_subj
list_sex_=ca_fc_list_sex
list_atlas_=list_atlas
list_covar_tanner_=ca_fc_list_covar_tanner
list_tanner_=ca_fc_list_tanner
list_covar_hormone_=ca_fc_list_covar_hormone
list_hormone_=ca_fc_list_hormone
list_dim_ca_=list_dim_ca
ratio_vis_=ratio_vis


print("Starting ca_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs_multi()",copy_log=T)
# Increase memory limit
memory.limit(1000000)
df_cor<-NULL

atlas<-"aal116"

# Group-wise analysis of factor-MRI matrix
print(paste("Calculating group-wise contribution of factors, atlas:",atlas,sep=" "))
dict_roi<-func_dict_roi(paths_)
dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]
list_group<-unique(dict_roi$group_3)

waves<-names(list_waves_)[1]
wave_mri<-list_waves_[[waves]]$wave_mri


df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output",
                                          paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var.csv",sep=""))))
df_ica_mri<-as.data.frame(fread(file.path(paths_$output,"output",
                                          paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep=""))))




#df_pca_mri_grp<-df_ica_mri_grp<-NULL

####
df_ca_mri<-df_pca_mri
dim<-40
####

df_ca_mri<-inner_join(df_ca_mri,dict_roi[,c("id","group_3")],by=c("from"="id"))
colnames(df_ca_mri)[colnames(df_ca_mri)=="group_3"]<-"from_group"
df_ca_mri<-inner_join(df_ca_mri,dict_roi[,c("id","group_3")],by=c("to"="id"))
colnames(df_ca_mri)[colnames(df_ca_mri)=="group_3"]<-"to_group"

df_ca_mri_grp<-NULL
for (label_sex in names(list_sex_)){
  for (idx_grp1 in seq(length(list_group))){
    for (idx_grp2 in seq(idx_grp1,length(list_group))){
      df_ca_mri_subset<-rbind(df_ca_mri[df_ca_mri$sex==label_sex
                                        & df_ca_mri$from_group==list_group[idx_grp1]
                                        & df_ca_mri$to_group==list_group[idx_grp2],],
                              df_ca_mri[df_ca_mri$sex==label_sex
                                        & df_ca_mri$from_group==list_group[idx_grp2]
                                        & df_ca_mri$to_group==list_group[idx_grp1],])
      df_ca_mri_subset<-df_ca_mri_subset[,sprintf("comp_%03d",seq(dim))]
      df_ca_mri_grp<-rbind(df_ca_mri_grp,
                           cbind(sex=label_sex,abs=F,from=list_group[idx_grp1],to=list_group[idx_grp2],
                                 t(colMeans(df_ca_mri_subset))))
      df_ca_mri_subset_abs<-abs(df_ca_mri_subset)
      df_ca_mri_grp<-rbind(df_ca_mri_grp,
                           cbind(sex=label_sex,abs=T,from=list_group[idx_grp1],to=list_group[idx_grp2],
                                 t(colMeans(df_ca_mri_subset_abs))))
      
    }
  }
}