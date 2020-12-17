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


####

#ca_fc_list_wave_mri<-ca_fc_list_wave_mri[1]
#ca_fc_list_wave_clin<-ca_fc_list_wave_clin[1]
ca_fc_list_tanner<-ca_fc_list_tanner[1]
ca_fc_list_hormone<-ca_fc_list_hormone[1]
list_atlas<-list_atlas[1]
