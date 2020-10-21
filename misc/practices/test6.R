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

print(paste("Calculating group-wise contribution to factors, atlas:",atlas,sep=" "))
dict_roi<-func_dict_roi(paths_)
dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]
list_group<-unique(dict_roi$group_3)
for (waves in names(list_waves_)){
  wave_mri<-list_waves_[[waves]]$wave_mri
  path_pca_mri_grp<-file.path(paths_$output,"output",
                              paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var_grp.csv",sep=""))
  path_ica_mri_grp<-file.path(paths_$output,"output",
                              paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var_grp.csv",sep=""))
  if (!(file.exists(path_pca_mri_grp) & file.exists(path_ica_mri_grp))){
    df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output",
                                              paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var.csv",sep=""))))
    df_ica_mri<-as.data.frame(fread(file.path(paths_$output,"output",
                                              paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep=""))))
    df_pca_mri_grp<-df_ica_mri_grp<-data.frame()
    # PCA
    dim<-max(list_dim_ca_)
    df_pca_mri_grp<-group_factor(df_pca_mri,dim,dict_roi,list_group,list_sex_)
    write.csv(df_pca_mri_grp,path_pca_mri_grp,row.names=F)
    # ICA
    for (dim in list_dim_ca_){
      df_ica_mri_grp<-rbind.fill(df_ica_mri_grp,
                                 group_factor(df_pca_mri,dim,dict_roi,list_group,list_sex_))
    }
    write.csv(df_ica_mri_grp,path_ica_mri_grp,row.names=F)
  }
}