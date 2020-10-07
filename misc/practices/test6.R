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



# Calculate factor attibution-clinical relation
for (waves in names(list_waves_)){
  wave_clin<-list_waves_[[waves]]$wave_clin
  wave_mri<-list_waves_[[waves]]$wave_mri
  
  print(paste("Clinical wave: ",wave_clin,", MRI wave: ",wave_mri,", loading PCA/ICA results.",sep=""))
  
  df_pca_subj<-read.csv(file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_subj.csv",sep="")))
  df_ica_subj<-read.csv(file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep="")))
  
  # Calculate factor-clinical correlation
  subset_subj_temp<-list(subset_subj_[[as.character(wave_mri)]])
  names(subset_subj_temp)<-wave_clin
  #1 Tanner stage
  for (idx_tanner in names(list_tanner_)){
    print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
    list_covar<-list_covar_tanner_
    list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
    n_covar<-length(list_covar)
    suffix<-paste("ses-",waves,"_var-",idx_tanner,sep="")
    
    data_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                       rem_na_clin=T,prefix=suffix,print_terminal=F)
    df_clin<-data_clin$df_clin
    df_clin$wave<-NULL
    
    # Calculate correlation between component attribution and clinical covariate
    data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                            n_covar=n_covar,list_sex=list_sex_,method="pca",
                            wave_mri=wave_mri,wave_clin=wave_clin,
                            idx_var=idx_tanner,label_var=list_tanner_[[idx_tanner]]$label)
    df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_tanner,
                               method="pca",data_cor$df_cor_flat))
    data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                            n_covar=n_covar,list_sex=list_sex_,method="ica",
                            wave_mri=wave_mri,wave_clin=wave_clin,
                            idx_var=idx_tanner,label_var=list_tanner_[[idx_tanner]]$label)
    df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_tanner,
                               method="ica",data_cor$df_cor_flat))
  } # End for Tanner stages
  
  #2 Hormones
  for (idx_hormone in names(list_hormone_)){
    print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
    list_covar<-list_covar_hormone_
    list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
    n_covar<-length(list_covar)
    suffix<-paste("ses-",waves,"_var-",idx_hormone,sep="")
    
    data_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                       rem_na_clin=T,prefix=suffix,print_terminal=F)
    df_clin<-data_clin$df_clin
    df_clin$wave<-NULL
    
    # Calculate correlation between component attribution and clinical covariate
    data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                            n_covar=n_covar,list_sex=list_sex_,method="pca",
                            wave_mri=wave_mri,wave_clin=wave_clin,
                            idx_var=idx_hormone,label_var=list_hormone_[[idx_hormone]]$label)
    df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_hormone,
                               method="pca",data_cor$df_cor_flat))
    data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                            n_covar=n_covar,list_sex=list_sex_,method="ica",
                            wave_mri=wave_mri,wave_clin=wave_clin,
                            idx_var=idx_hormone,label_var=list_hormone_[[idx_hormone]]$label)
    df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_hormone,
                               method="ica",data_cor$df_cor_flat))
  } # End for hormones
} # End for waves