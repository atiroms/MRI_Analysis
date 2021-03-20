#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"analyze/connection.R"))


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full<-NULL
#path_exp_full<-"/media/atiroms/HDD_05/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in<-"431_fc_aroma_gsr"
#dir_out<-"433_fc_gam_aroma_gsr"

list_dim_ca<-c(10,20,40)
#list_dim_ca<-10
#ratio_vis<-0.01

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100x7","schaefer200x7","schaefer400x7",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
#list_atlas<-c("aal116","gordon333","ho112","power264",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
#list_atlas<-c("aal116","ho112")
#list_atlas<-c("aal116","power264","shen268")
list_atlas<-"aal116"

#list_type_p=c("p","p_bh","seed_p_bh")
list_type_p=c("p_bh") # Benjamini-Hochberg method of FDR
thr_p <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(data.table)
library(stringr)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)



# OBSOLETE ======


##**************************************************
## Component analyses of FC in cross section =======
##**************************************************
# OBSOLETE
#ca_fc_cs<-function(paths_=paths,list_atlas_=list_atlas,wave_clin_=wave_clin,wave_mri_=wave_mri,
#                   list_covar_=list_covar,subset_subj_=subset_subj,list_dim_ca_=list_dim_ca,
#                   plot_result=F,suffix_=""){
#  
#  print("Starting ca_fc().")
#  #memory.limit(200000)
#  nullobj<-func_createdirs(paths_,str_proc="ca_fc()",copy_log=T)
#  
#  # Load and subset clinical data according to specified subsetting condition and covariate availability
#  print('Loading clinical data.')
#  data_clin<-func_clinical_data_long(paths_,wave_clin_,subset_subj_,list_covar_,
#                                     rem_na_clin=F,suffix=suffix_)
#  df_clin<-data_clin$df_clin
#  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
#  
#  for (atlas in list_atlas_){
#    # Load and examine FC data
#    print(paste("Loading FC of atlas: ",atlas,sep=""))
#    
#    # Create graph edge dataframe and node list
#    #df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
#    df_conn<-fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
#    df_conn<-df_conn[df_conn$ses==wave_mri_,]
#    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
#    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
#    n_edge<-dim(df_edge)[1]
#    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
#    n_node<-length(list_node)
#    
#    # Create list of subjects who meet subsetting condition and whose MRI data exist
#    list_subj_mri<-unique(df_conn[df_conn$ses==wave_mri_,]$ID_pnTTC)
#    list_subj_clin<-unique(df_clin[df_clin$ses==wave_clin_,]$ID_pnTTC)
#    list_subj_calc<-intersect(list_subj_mri,list_subj_clin)
#    print(paste("MRI data absent in",
#                as.character(length(list_subj_clin)-length(list_subj_calc)),
#                "subjects.",sep=" "))
#    
#    # Cbind FC data (Fisher-z transform of FC) as input for PCA function
#    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
#    df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
#    colnames(df_clin_exist)<-colnames(df_clin)
#    for (id_subj in list_subj_calc){
#      df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
#      df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
#      df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ID_pnTTC==id_subj,])
#    }
#    colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
#    rownames(df_conn_cbind)<-NULL
#    # Transpose connection dataframe (rows >> data for each subject/session, columns >> data for each edge)
#    df_conn<-as.data.frame(t(df_conn_cbind))
#    df_conn_cbind<-NULL
#    gc()
#    
#    # Calculate PCA of FC
#    print("Calculating PCA of FC.")
#    dim_ca<-max(list_dim_ca_)
#    data_pca<-func_pca(df_src=df_conn,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca)
#    write.csv(data_pca$df_comp_mri,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_variable.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_comp_subj,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_subject.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_vaf,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_variance.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_comp_clin_flat,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_correlation.csv",sep="")),row.names=F)
#    
#    # Plot PCA results
#    if (plot_result){
#      print("Plotting PCA of FC.")
#      #list_plot_pca<-plot_ca(df_src=data_pca$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_pca$dim)
#      #for (i_dim in names(list_plot_pca)){
#      #  for (name_covar in names(list_plot_pca[[i_dim]])){
#      #    plot<-list_plot_pca[[i_dim]][[name_covar]]
#      #    plot<-(plot
#      #           + ggtitle("PCA of FC"))
#      #    ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_pca.eps",sep=""),plot=plot,device=cairo_ps,
#      #           path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
#      #  }
#      #}
#    }
#    data_pca<-NULL
#    list_plot_pca<-NULL
#    gc()
#    
#    # Calculate ICA of FC
#    print("Calculating ICA of FC.")
#    df_comp_mri_rbind<-data.frame()
#    df_comp_subj_rbind<-data.frame()
#    df_variance_rbind<-data.frame()
#    df_comp_clin_flat_rbind<-data.frame()
#    for (dim_ca in list_dim_ca_){
#      data_ica<-func_ica(df_src=df_conn,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca)
#      df_comp_mri_rbind<-rbind.fill(df_comp_mri_rbind,cbind(dim=dim_ca,data_ica$df_comp_mri))
#      df_comp_subj_rbind<-rbind.fill(df_comp_subj_rbind,cbind(dim=dim_ca,data_ica$df_comp_subj))
#      df_variance_rbind<-rbind.fill(df_variance_rbind,cbind(dim=dim_ca,data_ica$df_variance))
#      df_comp_clin_flat_rbind<-rbind.fill(df_comp_clin_flat_rbind,cbind(dim=dim_ca,data_ica$df_comp_clin_flat))
#    }
#    write.csv(df_comp_mri_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_variable.csv",sep="")),row.names=F)
#    write.csv(df_comp_subj_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_subject.csv",sep="")),row.names=F)
#    write.csv(df_variance_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_variance.csv",sep="")),row.names=F)
#    write.csv(df_comp_clin_flat_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_correlation.csv",sep="")),row.names=F)
#    
#    # Plot ICA results
#    if (plot_result){
#      print("Plotting ICA of FC.")
#      list_plot_ica<-plot_ca(df_src=data_ica$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_ica$dim)
#      for (i_dim in names(list_plot_ica)){
#        for (name_covar in names(list_plot_ica[[i_dim]])){
#          plot<-list_plot_ica[[i_dim]][[name_covar]]
#          plot<-(plot
#                 + ggtitle("ICA of FC"))
#          ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_ica.eps",sep=""),plot=plot,device=cairo_ps,
#                 path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
#        }
#      }
#    }
#    data_ica<-NULL
#    list_plot_ica<-NULL
#    gc()
#  }
#  print("Finished ca_fc().")
#}


##**************************************************
## Iterate gam_fc_cs() over clinical variables =====
## and waves =======================================
##**************************************************
## succeeded in gam_fc_cs()
#gam_fc_cs_multi<-function(paths_=paths,list_atlas_=list_atlas,
#                          list_waves_=gam_fc_list_waves,subset_subj_=gam_fc_subset_subj,
#                          list_covar_tanner_=gam_fc_list_covar_tanner,list_tanner_=gam_fc_list_tanner,
#                          list_mod_tanner_=gam_fc_list_mod_tanner,list_plot_tanner_=gam_fc_list_plot_tanner,
#                          list_covar_hormone_=gam_fc_list_covar_hormone,list_hormone_=gam_fc_list_hormone,
#                          list_mod_hormone_=gam_fc_list_mod_hormone,list_plot_hormone_=gam_fc_list_plot_hormone,
#                          list_type_p_=list_type_p,thr_p_=thr_p){
#  print("Starting gam_fc_cs_multi()")
#  nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)
#  memory.limit(1000000)
#  
#  for (label_waves in names(list_waves_)){
#    wave_clin<-list_waves_[[label_waves]]$wave_clin
#    wave_mri<-list_waves_[[label_waves]]$wave_mri
#    waves<-list_waves_[label_waves]
#    #print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))
#    
#    # Prepare subject subsetting condition (MRI QC criteria) according to specified waves
#    subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
#    subset_subj_temp<-list(subset_subj_temp)
#    names(subset_subj_temp)<-wave_clin
#
#    #1 Tanner stage
#    for (idx_tanner in names(list_tanner_)){
#      print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,", Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
#      list_covar<-list_covar_tanner_
#      list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
#      #suffix<-paste("ses-",label_waves,"_var-",idx_tanner,sep="")
#      
#      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
#                         list_atlas_=list_atlas_,
#                         list_mod_=list_mod_tanner_,list_plot_=list_plot_tanner_,
#                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
#                         waves_=waves,idx_var_=idx_tanner)
#      nullobj<-NULL
#      gc()
#    } # finished looping over Tanner stages
#    
#    
#    #2 Hormones
#    for (idx_hormone in names(list_hormone_)){
#      print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,", Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
#      list_covar<-list_covar_hormone_
#      list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
#      #suffix<-paste("ses-",label_waves,"_var-",idx_hormone,sep="")
#      
#      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
#                         list_atlas_=list_atlas,
#                         list_mod_=list_mod_hormone_,list_plot_=list_plot_hormone_,
#                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
#                         waves_=waves,idx_var_=idx_hormone)
#      nullobj<-NULL
#      gc()
#    } # finished looping over Hormones
#    
#  } # finished looping over waves
#  
#  print("Combining results.")
#  df_gam<-df_gam_grp_sign<-df_gam_grp_abs<-df_gam_aic<-df_gam_grp_sign_aic<-df_gam_grp_abs_aic<-data.frame()
#  for (label_waves in names(list_waves_)){
#    for (idx_var in c(names(list_tanner_),names(list_hormone_))){
#      for (atlas in list_atlas_){
#        prefix<-paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var,sep="")
#        df_prefix<-data.frame("atlas"=atlas,"ses"=label_waves,"idx_var"=idx_var)
#        df_gam<-rbind(df_gam,cbind(df_prefix,
#                                   as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam.csv',sep=''))))))
#        df_gam_aic<-rbind(df_gam_aic,cbind(df_prefix,
#                                   as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam_aic.csv',sep=''))))))
#        df_gam_grp_sign<-rbind(df_gam_grp_sign,cbind(df_prefix,
#                                                     as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam_grp_sign.csv',sep=''))))))
#        df_gam_grp_abs<-rbind(df_gam_grp_abs,cbind(df_prefix,
#                                                   as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam_grp_abs.csv',sep=''))))))
#        df_gam_grp_sign_aic<-rbind(df_gam_grp_sign_aic,cbind(df_prefix,
#                                                     as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam_grp_sign_aic.csv',sep=''))))))
#        df_gam_grp_abs_aic<-rbind(df_gam_grp_abs_aic,cbind(df_prefix,
#                                                   as.data.frame(fread(file.path(paths_$output,'output','temp',paste(prefix,'_gam_grp_abs_aic.csv',sep=''))))))
#        
#      }
#    }
#  }
#  write.csv(df_gam,file.path(paths_$output,"output","result","gam.csv"),row.names = F)
#  write.csv(df_gam_aic,file.path(paths_$output,"output","result","gam_aic.csv"),row.names = F)
#  write.csv(df_gam_grp_sign,file.path(paths_$output,"output","result","gam_grp_sign.csv"),row.names = F)
#  write.csv(df_gam_grp_abs,file.path(paths_$output,"output","result","gam_grp_abs.csv"),row.names = F)
#  write.csv(df_gam_grp_sign_aic,file.path(paths_$output,"output","result","gam_grp_sign_aic.csv"),row.names = F)
#  write.csv(df_gam_grp_abs_aic,file.path(paths_$output,"output","result","gam_grp_abs_aic.csv"),row.names = F)
#  print("Finished gam_fc_cs_multi()")
#}
#
#
##**************************************************
## Additive/Linear model of FC in cross-section ====
##**************************************************
#gam_fc_cs<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
#                    list_atlas_=list_atlas,
#                    list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
#                    list_type_p_=list_type_p,thr_p_=thr_p,waves_,idx_var_
#                    ){
#  #print("Starting gam_fc_cs().")
#  #nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
#  dict_roi <- func_dict_roi(paths_)
#  
#  # Load and subset clinical data according to specified subsetting condition and covariate availability
#  #print('Loading clinical data.')
#  label_waves<-names(waves_)
#  wave_clin<-waves_[[1]]$wave_clin
#  wave_mri<-waves_[[1]]$wave_mri
#  data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin,subset_subj_,
#                                     list_covar=list_covar_,rem_na_clin=T,
#                                     prefix=paste("ses-",label_waves,"_var-",idx_var_,sep=""),
#                                     print_terminal=F)
#  df_clin<-data_clin$df_clin
#  
#  for (atlas in list_atlas_){
#    # check if the last file in the loop is already generated
#    path_file_check<-file.path(paths_$output,"output","temp",
#                               paste("atl-",atlas,"_ses-",label_waves,
#                                     "_var-",idx_var_,"_gam_grp_abs_src.csv",sep=""))
#    if (!file.exists(path_file_check)){
#      # Load ROI-wise FC data
#      df_fc<-as.data.frame(fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep=''))))
#      df_join<-join_fc_clin_cs(df_fc,df_clin,wave_clin,wave_mri)
#      
#      # Calculate and save ROI-wise GAMM of FC
#      print(paste('Calculating GAM, atlas: ',atlas,sep=''))
#      list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
#      df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
#      df_roi$id<-as.character(df_roi$id)
#      colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
#      data_gamm<-iterate_gamm(df_join,df_roi,list_mod_,calc_parallel=T,calc_identical=F)
#      df_gam<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod_,list_plot_,calc_seed_level=F)
#      write.csv(df_gam,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam.csv",sep=""))
#                ,row.names = F)
#      write.csv(data_gamm$df_out_aic,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_aic.csv",sep="")),
#                row.names = F)
#      write.csv(df_join,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_src.csv",sep="")),
#                row.names=F)
#      
#      # Calculate group-wise average of FC
#      df_fc<-df_fc[df_fc$ses==wave_mri,]
#      df_fc<-left_join(df_fc,df_roi,by=c("from"="id"))
#      colnames(df_fc)[colnames(df_fc)=="group"]<-"from_grp"
#      df_fc<-left_join(df_fc,df_roi,by=c("to"="id"))
#      colnames(df_fc)[colnames(df_fc)=="group"]<-"to_grp"
#      list_id_subj<-sort(unique(df_fc$ID_pnTTC))
#      list_group<-unique(df_roi$group)
#      df_fc_grp<-data.frame()
#      for (id_subj in list_id_subj){
#        df_fc_subset1<-df_fc[df_fc$ID_pnTTC==id_subj,]
#        for (idx_grp1 in seq(length(list_group))){
#          for (idx_grp2 in seq(idx_grp1,length(list_group))){
#            label_grp1<-as.character(list_group[idx_grp1])
#            label_grp2<-as.character(list_group[idx_grp2])
#            df_fc_subset2<-df_fc_subset1[df_fc_subset1$from_grp==label_grp1 & df_fc_subset1$to_grp==label_grp2,]
#            df_fc_grp<-rbind(df_fc_grp,data.frame("ID_pnTTC"=id_subj,"from"=label_grp1,"to"=label_grp2,
#                                                  "mean_z_r"=mean(df_fc_subset2$z_r),
#                                                  "mean_abs_z_r"=mean(abs(df_fc_subset2$z_r))))
#          }
#        }
#      }
#      
#      # Join group-wise FC and clinical data
#      df_fc_grp$mean_z_r[which(is.nan(df_fc_grp$mean_z_r))]<-0
#      df_fc_grp$mean_abs_z_r[which(is.nan(df_fc_grp$mean_abs_z_r))]<-0
#      df_clin_wave<-df_clin[df_clin$wave==wave_clin,]
#      df_clin_wave$wave<-NULL
#      df_join<-inner_join(df_fc_grp,df_clin_wave,by='ID_pnTTC')
#      for (key in c('ID_pnTTC','sex')){
#        if (key %in% colnames(df_join)){
#          df_join[,key]<-as.factor(df_join[,key])
#        }
#      }
#      
#      # Calculate and save group-wise GAMM of FC
#      df_grp<-data.frame("id"=list_group,"label"=str_to_title(gsub("_"," ",list_group)))
#      df_join_sign<-df_join_abs<-df_join
#      colnames(df_join_sign)[colnames(df_join_sign)=="mean_z_r"]<-"value"
#      df_join_sign$mean_abs_z_r<-NULL
#      data_gamm_grp_sign<-iterate_gamm(df_join_sign,df_grp,list_mod_,calc_parallel=T,calc_identical=T)
#      df_gam_grp_sign<-add_mltcmp(data_gamm_grp_sign$df_out_gamm,df_grp,list_mod_,list_plot_,calc_seed_level=F)
#      write.csv(df_gam_grp_sign,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_sign.csv",sep="")),
#                row.names = F)
#      write.csv(data_gamm_grp_sign$df_out_aic,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_sign_aic.csv",sep="")),
#                row.names = F)
#      write.csv(df_join_sign,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_sign_src.csv",sep="")),
#                row.names = F)
#      colnames(df_join_abs)[colnames(df_join_abs)=="mean_abs_z_r"]<-"value"
#      df_join_abs$mean_z_r_<-NULL
#      data_gamm_grp_abs<-iterate_gamm(df_join_abs,df_grp,list_mod_,calc_parallel=T,calc_identical=T)
#      df_gam_grp_abs<-add_mltcmp(data_gamm_grp_abs$df_out_gamm,df_grp,list_mod_,list_plot_,calc_seed_level=F)
#      write.csv(df_gam_grp_abs,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_abs.csv",sep="")),
#                row.names = F)
#      write.csv(data_gamm_grp_abs$df_out_aic,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_abs_aic.csv",sep="")),
#                row.names = F)
#      write.csv(df_join_abs,
#                file.path(paths_$output,"output","temp",
#                          paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_abs_src.csv",sep="")),
#                row.names = F)
#      
#      # Graphical output of ROI-wise and group-wise GAMM of FC
#      plot_gam_fc(paths_,df_gam,df_gam_grp_sign,df_gam_grp_abs,atlas,
#                  list_mod=list_mod_,list_plot=list_plot_,
#                  list_type_p=list_type_p_,thr_p=thr_p_,waves=waves_,idx_var=idx_var_)
#    }
#  } # End for atlas
#}
#