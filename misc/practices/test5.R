paths_=paths
list_waves_=ca_fc_list_waves
subset_subj_=ca_fc_subset_subj
list_atlas_=list_atlas
list_covar_tanner_=ca_fc_list_covar_tanner
list_tanner_=ca_fc_list_tanner
list_covar_hormone_=ca_fc_list_covar_hormone
list_hormone_=ca_fc_list_hormone
list_dim_ca_=list_dim_ca
plot_result=F



print("Starting ca_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs_multi()",copy_log=T)
df_cor<-df_ca_subj_bind<-df_ca_var_bind<-df_ca_vaf_bind<-NULL



atlas<-list_atlas_[1]



# Load and examine FC data
print(paste("Loading FC of atlas: ",atlas,sep=""))
df_conn<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep=""))))

# Create graph edge dataframe and node list
df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
n_edge<-dim(df_edge)[1]
list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
n_node<-length(list_node)

wave_mri_done<-NULL



waves<-names(list_waves_)[[1]]



wave_clin<-list_waves_[[waves]]$wave_clin
wave_mri<-list_waves_[[waves]]$wave_mri
print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))
if (wave_mri=="2-1"){
  label_wave_mri<-"m21"
}else{
  label_wave_mri<-paste("m",wave_mri,sep="")
}


wave_mri_done<-c(wave_mri_done,wave_mri)

# Prepare subject subsetting condition (MRI QC criteria) according to specified waves
subset_subj_temp<-list(subset_subj_[[as.character(wave_mri)]])
names(subset_subj_temp)<-wave_mri
data_clin<-func_clinical_data_long(paths_,wave_mri,subset_subj_temp,list_covar=NULL,
                                   rem_na_clin=F,prefix=paste("wave-",label_wave_mri,sep=""),
                                   print_terminal=F)
df_clin<-data_clin$df_clin
colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"

# Create list of subjects who meet subsetting condition and whose MRI data exist
df_conn_ses<-df_conn[df_conn$ses==wave_mri,]
list_subj_mri<-unique(df_conn_ses$ID_pnTTC)
list_subj_qc<-unique(df_clin[df_clin$ses==wave_mri,]$ID_pnTTC)
list_subj_calc<-intersect(list_subj_mri,list_subj_qc)
n_subj_del<-length(list_subj_qc)-length(list_subj_calc)
if (n_subj_del>0){
  print(paste("MRI data absent in",as.character(n_subj_del),"subjects.",sep=" "))
}

# Cbind FC data (Fisher-z transform of FC) as input for PCA function
df_conn_calc<-data.frame(matrix(nrow=n_edge,ncol=0))
df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
colnames(df_clin_exist)<-colnames(df_clin)
for (id_subj in list_subj_calc){
  df_conn_subj<-df_conn_ses[which(df_conn_ses$ID_pnTTC==id_subj),]
  df_conn_calc<-cbind(df_conn_calc,df_conn_subj[["z_r"]])
  df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ID_pnTTC==id_subj,])
}
colnames(df_conn_calc)<-as.character(seq(ncol(df_conn_calc)))
rownames(df_conn_calc)<-NULL
# Transpose connection dataframe (rows >> data for each subject/session, columns >> data for each edge)
df_conn_calc<-as.data.frame(t(df_conn_calc))
df_conn_ses<-NULL
gc()

# Calculate PCA of FC
dim_ca<-max(list_dim_ca_)



df_src=df_conn_calc
df_var=df_edge
df_indiv=df_clin_exist
dim_ca=dim_ca
calc_corr=F



if (sum(is.na(df_src))>0){
  # Estimate number of dimensions
  if (is.null(dim_ca)){
    ncp_estim<-estim_ncpPCA(df_src,ncp.max=ncol(df_src))$ncp
    ncp_calc<-ncp_estim
  }else{
    ncp_estim<-estim_ncpPCA(df_src,ncp.max=dim_ca)$ncp
    if (ncp_estim==dim_ca){
      print(paste("PCA data dimension may be greater than: ",as.character(ncp_estim),sep=""))
    }
    ncp_calc<-dim_ca
  }
  # Impute data
  df_src<-imputePCA(df_src,ncp=ncp_calc)$completeObs
}else{
  ncp_calc<-dim_ca
}

print(paste("Calculating PCA, dimension: ",as.character(ncp_calc),sep=""))
# PCA calculation
data_pca<-PCA(df_src,scale.unit = TRUE, ncp = ncp_calc, graph = FALSE)

# Component-imaging variable matrix
# Row: MRI variable, Column: component(factor)
df_comp_mri<-data.frame(data_pca$var$coord)
if(!is.null(df_var)){
  df_comp_mri<-cbind(df_var,df_comp_mri)
  colnames(df_comp_mri)<-c(colnames(df_var),sprintf("comp_%03d",1:ncp_calc))
}else{
  colnames(df_comp_mri)<-sprintf("comp_%03d",1:ncp_calc)
}
rownames(df_comp_mri)<-NULL





df_comp_mri<-rownames_to_column(df_comp_mri,"id")
df_comp_mri$id<-as.numeric(df_comp_mri$id)
for (id_comp in 1:ncp_calc){
  df_comp_mri$abs<-abs(df_comp_mri[[sprintf("comp_%03d",id_comp)]])
  df_comp_mri<-df_comp_mri[order(df_comp_mri$abs,decreasing=T),]
  df_comp_mri[[sprintf("rank_%03d",id_comp)]]<-1:nrow(df_comp_mri)
}
df_comp_mri<-df_comp_mri[order(df_comp_mri$id),]
df_comp_mri<-df_comp_mri[c(colnames(df_var),
                            sprintf("comp_%03d",1:ncp_calc),
                            sprintf("rank_%03d",1:ncp_calc))]
