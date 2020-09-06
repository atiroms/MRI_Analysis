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

if (wave_mri %in% wave_mri_done){
  print(paste("Clinical: ", wave_clin,", MRI: ",wave_mri,", loading PCA/ICA results.",sep=""))
  df_pca_subj<-read.csv(file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_subj.csv",sep="")))
  df_ica_subj<-read.csv(file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep="")))
}else{
  wave_mri_done<-c(wave_mri_done,wave_mri)
  df_pca_mri<-df_pca_subj<-df_pca_vaf<-df_ica_mri<-df_ica_subj<-df_ica_vaf<-NULL
  for (label_sex in names(list_sex_)){
    
    # Prepare subject subsetting condition (MRI QC criteria and sex) according to specified mri wave
    subset_subj_temp<-list(c(subset_subj_[[as.character(wave_mri)]],
                             list(list("key"="Sex","condition"=list_sex_[[label_sex]]))))
    names(subset_subj_temp)<-wave_mri
    data_clin<-func_clinical_data_long(paths_,wave_mri,subset_subj_temp,list_covar=NULL,
                                       rem_na_clin=F,
                                       prefix=paste("ses-m",wave_mri,"_sex-",label_sex,sep=""),
                                       print_terminal=F)
    df_clin<-data_clin$df_clin
    colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
    print(paste("Clinical: ", wave_clin,", MRI: ",wave_mri,
                ", Sex: ",label_sex," PCA/ICA.",sep=""))
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    df_conn_ses<-df_conn[df_conn$ses==wave_mri,]
    list_subj_mri<-unique(df_conn_ses$ID_pnTTC)
    list_subj_qc<-unique(df_clin[df_clin$ses==wave_mri,]$ID_pnTTC)
    list_subj_calc<-intersect(list_subj_mri,list_subj_qc)
    
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
    data_pca<-func_pca(df_src=df_conn_calc,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca,calc_corr=F)
    df_pca_mri<-rbind(df_pca_mri,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_mri))
    df_pca_subj<-rbind(df_pca_subj,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_subj))
    df_pca_vaf<-rbind(df_pca_vaf,cbind(sex=label_sex,dim=dim_ca,data_pca$df_vaf))
    df_ca_var_bind<-rbind(df_ca_var_bind,cbind(atlas=atlas,method="pca",ses=wave_mri,df_pca_mri))
    df_ca_subj_bind<-rbind(df_ca_subj_bind,cbind(atlas=atlas,method="pca",ses=wave_mri,df_pca_subj))
    df_ca_vaf_bind<-rbind(df_ca_vaf_bind,cbind(atlas=atlas,method="pca",ses=wave_mri,df_pca_vaf))
    
    data_pca<-NULL
    gc()
    
    # Calculate ICA of FC
    for (dim_ca in list_dim_ca_){
      data_ica<-func_ica(df_src=df_conn_calc,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca,calc_corr=F)
      df_ica_mri<-rbind.fill(df_ica_mri,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_mri))
      df_ica_subj<-rbind.fill(df_ica_subj,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_subj))
      df_ica_vaf<-rbind.fill(df_ica_vaf,cbind(sex=label_sex,dim=dim_ca,data_ica$df_vaf))
    }
    df_ca_var_bind<-rbind(df_ca_var_bind,cbind(atlas=atlas,method="ica",ses=wave_mri,df_ica_mri))
    df_ca_subj_bind<-rbind(df_ca_subj_bind,cbind(atlas=atlas,method="ica",ses=wave_mri,df_ica_subj))
    df_ca_vaf_bind<-rbind.fill(df_ca_vaf_bind,cbind(atlas=atlas,method="ica",ses=wave_mri,df_ica_vaf))
    
    data_ica<-NULL
    gc()
  } # end of loop over sex
  write.csv(df_pca_mri,file.path(paths_$output,"output",
                                 paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var.csv",sep="")),row.names=F)
  write.csv(df_pca_subj,file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_subj.csv",sep="")),row.names=F)
  write.csv(df_pca_vaf,file.path(paths_$output,"output",
                                 paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_vaf.csv",sep="")),row.names=F)
  write.csv(df_ica_mri,file.path(paths_$output,"output",
                                 paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep="")),row.names=F)
  write.csv(df_ica_subj,file.path(paths_$output,"output",
                                  paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep="")),row.names=F)
  write.csv(df_ica_vaf,file.path(paths_$output,"output",
                                 paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_vaf.csv",sep="")),row.names=F)
} # end if wave_mri is not in wave_mri_done



subset_subj_temp<-list(subset_subj_[[as.character(wave_mri)]])
names(subset_subj_temp)<-wave_clin
idx_tanner<-names(list_tanner_)[[1]]



print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
n_covar<-length(list_covar)
suffix<-paste("ses-",waves,"_var-",idx_tanner,sep="")

data_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                   rem_na_clin=T,prefix=suffix,print_terminal=F)
df_clin<-data_clin$df_clin
df_clin$wave<-NULL



df_comp_subj=df_pca_subj
list_sex<-list_sex_


list_dim<-sort(unique(df_comp_subj$dim))
df_cor_rbind<-df_cor_flat_rbind<-NULL

dim<-list_dim[1]
label_sex<-names(list_sex)[[1]]



df_comp_subj_subset<-df_comp_subj[df_comp_subj$dim==dim & df_comp_subj$sex==label_sex,
                                  c("ID_pnTTC",sprintf("comp_%03d",1:dim))]
df_join<-inner_join(df_clin,df_comp_subj_subset,by="ID_pnTTC")
df_join$ID_pnTTC<-NULL
data_cor<-func_cor(df_join)
df_cor<-data_cor$cor
df_cor<-df_cor[(n_covar+1):nrow(df_cor),1:n_covar]
df_cor<-rownames_to_column(df_cor,"comp")
df_cor$comp<-sapply(sapply(df_cor$comp,substr,start=6,stop=8),as.integer)
df_cor<-cbind(dim=dim,df_cor)
df_cor_flat<-data_cor$cor_flat
df_cor_flat<-df_cor_flat[df_cor_flat$from %in% colnames(df_clin) & df_cor_flat$to %in% colnames(df_comp_subj_dim),]
df_cor_flat<-df_cor_flat[,c("from","to","r","p")]
colnames(df_cor_flat)<-c("covar","comp","r","p")
df_cor_flat$comp<-sapply(sapply(df_cor_flat$comp,substr,start=6,stop=8),as.integer)
df_cor_flat<-cbind(dim=dim,df_cor_flat)
df_cor_rbind<-rbind(df_cor_rbind,cbind(sex=label_sex,df_cor))
df_cor_flat_rbind<-rbind(df_cor_flat_rbind,cbind(sex=label_sex,df_cor_flat))