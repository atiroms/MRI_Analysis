paths_=paths
list_wave_mri_=ca_fc_list_wave_mri
list_wave_clin_=ca_fc_list_wave_clin
subset_subj_=ca_fc_subset_subj
list_sex_=ca_fc_list_sex
list_atlas_=list_atlas
list_covar_tanner_=ca_fc_list_covar_tanner
list_tanner_=ca_fc_list_tanner
list_covar_hormone_=ca_fc_list_covar_hormone
list_hormone_=ca_fc_list_hormone
list_dim_ca_=list_dim_ca


####

print("Starting ca_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs_multi()",copy_log=T)
# Increase memory limit for later ICA calculation
memory.limit(1000000)
df_cor<-NULL

####

atlas<-list_atlas_[1]

####

# Load and examine FC data
print(paste("Loading FC of atlas: ",atlas,sep=""))
df_conn<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep=""))))

# Create graph edge dataframe and node list
df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
n_edge<-dim(df_edge)[1]
list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
n_node<-length(list_node)

####

label_wave_mri<-names(list_wave_mri_)[1]

####

wave_mri<-list_wave_mri_[[label_wave_mri]]
print(paste("MRI wave: ",wave_mri,", loading PCA/ICA results.",sep=""))
df_pca_mri<-read.csv(file.path(paths_$output,"output","temp",
                               paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var.csv",sep="")))
df_pca_mri_grp<-read.csv(file.path(paths_$output,"output","temp",
                                   paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var_grp.csv",sep="")))

####

label_sex<-names(list_sex_)[1]

####

dim_ca<-max(list_dim_ca_)
df_pca_mri_subset<-df_pca_mri[df_pca_mri$sex==label_sex & df_pca_mri$dim==dim_ca,]
df_pca_mri_grp_subset<-df_pca_mri_grp[df_pca_mri_grp$sex==label_sex & df_pca_mri_grp$dim==dim_ca,]
df_pca_mri_subset$sex<-df_pca_mri_subset$dim<-df_pca_mri_grp_subset$sex<-df_pca_mri_grp_subset$dim<-NULL
# Visualize factor-FC matrix in heatmap plot
plot_ca_fc_heatmap(paths_=paths_,df_pca_mri_subset,df_pca_mri_grp_subset,atlas=atlas,dim_ca=dim_ca,
                   method="pca",label_sex=label_sex,ses=wave_mri)

####

label_wave_mri<-names(list_wave_mri_)[1]

####

wave_mri<-list_wave_mri_[[label_wave_mri]]
df_pca_subj<-read.csv(file.path(paths_$output,"output","temp",
                                paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_subj.csv",sep="")))
df_ica_subj<-read.csv(file.path(paths_$output,"output","temp",
                                paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep="")))

####

label_wave_clin<-names(list_wave_clin_)[1]

####

wave_clin<-list_wave_clin_[[label_wave_clin]]
print(paste("Clinical wave: ",wave_clin,", MRI wave: ",wave_mri,", calculating factor-clinical correlation.",sep=""))
# Pickup subsetting condition of MRI wave, and rename it to clinical wave
subset_subj_temp<-list(subset_subj_[[as.character(wave_mri)]])
names(subset_subj_temp)<-wave_clin

####

#1 Tanner stage
for (idx_tanner in names(list_tanner_)){
  #print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
  list_covar<-list_covar_tanner_
  list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
  n_covar<-length(list_covar)
  prefix<-paste("ses-c",as.character(wave_clin),"m",as.character(wave_mri),"_var-",idx_tanner,sep="")
  data_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                     rem_na_clin=T,prefix=prefix,print_terminal=F)
  df_clin<-data_clin$df_clin
  df_clin$wave<-NULL
  
  # Calculate correlation between component attribution and clinical covariate
  data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                          n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="pca",
                          wave_mri=wave_mri,wave_clin=wave_clin,
                          idx_var=idx_tanner,label_var=list_tanner_[[idx_tanner]]$label)
  df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
                             variable=idx_tanner,method="pca",data_cor$df_cor_flat))
  data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                          n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="ica",
                          wave_mri=wave_mri,wave_clin=wave_clin,
                          idx_var=idx_tanner,label_var=list_tanner_[[idx_tanner]]$label)
  df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
                             variable=idx_tanner,method="ica",data_cor$df_cor_flat))
} # End for Tanner stages

#2 Hormones
for (idx_hormone in names(list_hormone_)){
  #print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
  list_covar<-list_covar_hormone_
  list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
  n_covar<-length(list_covar)
  prefix<-paste("ses-c",as.character(wave_clin),"m",as.character(wave_mri),"_var-",idx_hormone,sep="")
  data_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                     rem_na_clin=T,prefix=prefix,print_terminal=F)
  df_clin<-data_clin$df_clin
  df_clin$wave<-NULL
  
  # Calculate correlation between component attribution and clinical covariate
  data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                          n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="pca",
                          wave_mri=wave_mri,wave_clin=wave_clin,
                          idx_var=idx_hormone,label_var=list_hormone_[[idx_hormone]]$label)
  df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
                             variable=idx_hormone,method="pca",data_cor$df_cor_flat))
  data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                          n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="ica",
                          wave_mri=wave_mri,wave_clin=wave_clin,
                          idx_var=idx_hormone,label_var=list_hormone_[[idx_hormone]]$label)
  df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
                             variable=idx_hormone,method="ica",data_cor$df_cor_flat))
} # End for hormones


####

#ca_fc_list_wave_mri<-ca_fc_list_wave_mri[1]
#ca_fc_list_wave_clin<-ca_fc_list_wave_clin[1]
#ca_fc_list_tanner<-ca_fc_list_tanner[1]
#ca_fc_list_hormone<-ca_fc_list_hormone[1]
list_atlas<-list_atlas[1]
