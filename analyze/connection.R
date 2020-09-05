#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full<-NULL


#dir_in<-dir_out<-"401_fc_acompcor"
dir_in<-"421_fc_aroma"
dir_out<-"422_fp_aroma"
list_atlas<-c("aal116","glasser360","gordon333","power264",
              "schaefer100x7","schaefer200x7","schaefer400x7",
              "schaefer100x17","schaefer200x17","schaefer400x17",
              "shen268")
path_exp_full<-"/media/atiroms/HDD_04/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in<-"421_fc_aroma"
#dir_out<-"427_fc_gamm_aroma"
#list_atlas<-"aal116"
#list_atlas<-"power264"
#list_atlas<-c("aal116","gordon333","power264","schaefer400x7","shen268")
#path_exp_full<-"/media/veracrypt1/MRI_img/pnTTC/puberty/stats/func_XCP"

#path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_CONN"
#dir_in<-"56.2_fc"
#dir_out<-"56.3_gamm_fc"
#list_atlas<-c("cnn","hoa","power264")
#list_atlas<-"cnn"
#list_atlas<-"hoa"

list_wave <- c(1,2)

list_type_p=c("p","p_bh","seed_p_bh")
#list_type_p="p"
thr_p <- 0.05

list_cost<-seq(0.15,0.40,0.01)
absolute<-T
threshold<-NA

list_dim_ca<-c(5,10,20,40)


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(GGally)
library(igraph)
#library(qgraph)
library(ggrepel)
library(colorRamps)
library(tidyverse)
library(parallel)
library(mgcv)
library(car)
library(plyr)
library(dplyr)
library(data.table)
library(pbapply)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/gta_function.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# GLM/GAM of FCs ==================================
#**************************************************

gamm_fc_core<-function(paths_,df_fc,atlas,df_roi,list_wave_,subset_subj_,
                       list_covar,list_mod,list_plot,idx_var,
                       list_type_p_=list_type_p,thr_p_=thr_p
                       ){
  # Prepare clinical data
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T,
                                     prefix=paste("var-",idx_var,sep=""),print_terminal=F)
  df_clin<-data_clin$df_clin
  
  # Join fc and clinical data
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
  data_gamm<-iterate_gamm(df_join,df_roi,list_mod)
  
  # Calculate multiple comparison-corrected p values
  df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod,
                           list_plot,calc_seed_level=T)
  
  # Save results
  write.csv(data_gamm$df_out_gamm,
            file.path(paths_$output,"output",
                      paste("atl-",atlas,"_var-",idx_var,"_fc_gamm.csv",sep="")),row.names = F)
  write.csv(data_gamm$df_out_aic,
            file.path(paths_$output,"output",
                      paste("atl-",atlas,"_var-",idx_var,"_fc_gamm_aic.csv",sep="")),row.names = F)
  write.csv(df_plot_gamm,
            file.path(paths_$output,"output",
                      paste("atl-",atlas,"_var-",idx_var,"_fc_gamm_plot.csv",sep="")),row.names = F)
  
  # Graphical output of ROI-wise GAMM of FC
  plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
              list_type_p_,thr_p,paths_,suffix_=paste("var-",idx_var,sep=""))
}

gamm_fc_multi<-function(paths_=paths,subset_subj_=gamm_fc_subset_subj,list_wave_=list_wave,
                        list_atlas_=list_atlas,key_group_='group_3',
                        list_covar_tanner_=gamm_fc_list_covar_tanner,list_tanner_=gamm_fc_list_tanner,
                        list_mod_tanner_=gamm_fc_list_mod_tanner,list_plot_tanner_=gamm_fc_list_plot_tanner,
                        list_covar_hormone_=gamm_fc_list_covar_hormone,list_hormone_=gamm_fc_list_hormone,
                        list_mod_hormone_=gamm_fc_list_mod_hormone,list_plot_hormone_=gamm_fc_list_plot_hormone){
  
  print("Starting gamm_fc_multi().")
  nullobj<-func_createdirs(paths_,str_proc="gamm_fc_multi()",copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    print(paste("Loading FC of atlas: ",atlas,sep=""))
    df_fc<-as.data.frame(fread(file.path(paths_$input,"output",
                                         paste("atl-",atlas,"_fc.csv",sep=""))))
    
    # Prepare dataframe of ROIs
    list_roi<-sort(unique(c(as.character(df_fc$from),as.character(df_fc$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_tanner in names(list_tanner_)){
      print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
      list_covar<-list_covar_tanner_
      list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
      gamm_fc_core(paths_,df_fc,atlas,df_roi,list_wave_,subset_subj_,
                   list_covar,list_mod_tanner_,list_plot_tanner_,idx_tanner)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(list_hormone_)){
      print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
      list_covar<-list_covar_hormone_
      list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
      gamm_fc_core(paths_,df_fc,atlas,df_roi,list_wave_,subset_subj_,
                   list_covar,list_mod_hormone_,list_plot_hormone_,idx_hormone)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  df_gamm<-df_aic<-df_gamm_plot<-data.frame()
  for (atlas in list_atlas_){
    #1 Tanner stage
    for (idx_tanner in names(list_tanner_)){
      df_gamm_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                 paste("atl-",atlas,"_var-",idx_tanner,
                                                       "_fc_gamm.csv",sep=""))))
      df_gamm<-rbind(df_gamm,cbind(atlas=atlas,variable=idx_tanner,df_gamm_add))
      df_aic_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                paste("atl-",atlas,"_var-",idx_tanner,
                                                      "_fc_gamm_aic.csv",sep=""))))
      df_aic<-rbind(df_aic,cbind(atlas=atlas,variable=idx_tanner,df_aic_add))
      df_gamm_plot_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                      paste("atl-",atlas,"_var-",idx_tanner,
                                                            "_fc_gamm_plot.csv",sep=""))))
      df_gamm_plot<-rbind(df_gamm_plot,cbind(atlas=atlas,variable=idx_tanner,df_gamm_plot_add))
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(list_hormone_)){
      df_gamm_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                 paste("atl-",atlas,"_var-",idx_hormone,
                                                       "_fc_gamm.csv",sep=""))))
      df_gamm<-rbind(df_gamm,cbind(atlas=atlas,variable=idx_hormone,df_gamm_add))
      df_aic_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                paste("atl-",atlas,"_var-",idx_hormone,
                                                      "_fc_gamm_aic.csv",sep=""))))
      df_aic<-rbind(df_aic,cbind(atlas=atlas,variable=idx_hormone,df_aic_add))
      df_gamm_plot_add<-as.data.frame(fread(file.path(paths_$output,"output",
                                                      paste("atl-",atlas,"_var-",idx_hormone,
                                                            "_fc_gamm_plot.csv",sep=""))))
      df_gamm_plot<-rbind(df_gamm_plot,cbind(atlas=atlas,variable=idx_hormone,df_gamm_plot_add))
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  write.csv(df_gamm,file.path(paths_$output,"output","fc_gamm.csv"),row.names = F)
  write.csv(df_aic,file.path(paths_$output,"output","fc_gamm_aic.csv"),row.names = F)
  write.csv(df_gamm_plot,file.path(paths_$output,"output","fc_gamm_plot.csv"),row.names = F)
}

join_fc_clin<-function(df_fc,df_clin){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                  -which(colnames(df_fc)=="p"))]
  df_clin$wave<-as.character(df_clin$wave)
  
  # Join clinical and FC data frames
  print('Joining clinical and FC data.')
  df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  return(df_join)
}

gamm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                  list_wave_=list_wave,list_atlas_=list_atlas,
                  #list_measure_=list_measure,list_str_group_=list_str_group,
                  list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
                  list_type_p_=list_type_p,thr_p_=thr_p
                  ){
  print("Starting gamm_fc().")
  nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
  
  for (atlas in list_atlas_){
    
    #****************************
    # ROI-wise FC GAMM calculation
    #****************************
    # Load ROI-wise FC data
    print(paste('Loading FC data, atlas:',atlas,sep=' '))
    df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_join<-join_fc_clin(df_fc,df_clin)
    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save ROI-wise GAMM of FC
    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
    write.csv(data_gamm$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
                             list_mod,list_plot,calc_seed_level=T)
    write.csv(df_plot_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_plt.csv",sep="")),row.names = F)
    
    # Graphical output of ROI-wise GAMM of FC
    plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)
    
    #****************************
    # Group-wise FC GAMM calculation
    #****************************
    # Load group-wise FC data
    print(paste('Loading group FC data, atlas:',atlas,sep=' '))
    df_fc_grp<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc_grp.csv',sep='')))
    df_join_grp<-join_fc_clin(df_fc_grp,df_clin)
    write.csv(df_join_grp,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save group-wise GAMM of FC
    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
    list_roi_grp<-sort(unique(c(as.character(df_join_grp$from),as.character(df_join_grp$to))))
    #df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    df_roi_grp<-data.frame(id=list_roi_grp,label=capitalize(list_roi_grp),group="group")
    data_gamm_grp<-iterate_gamm(df_join_grp,df_roi_grp,list_mod_)
    write.csv(data_gamm_grp$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm.csv",sep="")),row.names = F)
    write.csv(data_gamm_grp$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_roi_grp,analysis="grp",atlas,
                                 list_mod,list_plot,calc_seed_level=T)
    write.csv(df_plot_gamm_grp,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_plt.csv",sep="")),row.names = F)
    
    # Graphical output of group-wise GAMM of FC
    plot_gam_fc(df_plot_gamm_grp,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)
    
    #****************************
    # Multi-scale FC GAMM calculation
    #****************************
    # Subset ROI-wise GAMM result to include only within-group connections
    df_gamm_ms<-NULL
    for (group in list_roi_grp){
      list_roi_within_grp<-as.character(df_roi[df_roi$group==group,"id"])
      df_gamm_ms_add<-data_gamm$df_out_gamm[which(is.element(as.character(data_gamm$df_out_gamm[,"from"]),list_roi_within_grp)
                                            & is.element(as.character(data_gamm$df_out_gamm[,"to"]),list_roi_within_grp)),]
      df_gamm_ms_add<-cbind(group=group,df_gamm_ms_add)
      df_gamm_ms<-rbind(df_gamm_ms,df_gamm_ms_add)
    }
    
    # Combine within-group ROI-wise GAMM results and between-group GAMM results
    df_gamm_ms<-rbind(df_gamm_ms,cbind(group="group",data_gamm_grp$df_out_gamm))
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm_ms<-add_mltcmp(df_gamm_ms,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
                                calc_seed_level=F)
    write.csv(df_plot_gamm_ms,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-ms_gamm_plt.csv",sep="")),row.names = F)
    
    # Split data into ROI-wise and group-wise GAMM results, graphical output
    for (group in list_roi_grp){
      df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group==group,-1]
      df_roi_split<-df_roi[df_roi$group==group,]
      label_analysis<-paste("ms_grp-",group,sep="")
      plot_gamm_fc(df_plot_gamm_ms_split,df_roi_split,analysis=label_analysis,atlas,list_mod,list_plot,
                   list_type_p_,thr_p,paths_)
    }
    df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group=="group",-1]
    plot_gam_fc(df_plot_gamm_ms_split,df_roi_grp,analysis="ms_grp-group",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)

  }
  print('Finished gamm_fc().')
}


#**************************************************
# Calculate longitudinal FC difference ============
#**************************************************

diff_fc<-function(paths_=paths,list_atlas_=list_atlas,key_roigroup="group_3"){
  print("Starting diff_fc().")
  #nullobj<-func_createdirs(paths_,str_proc="diff_fc()")
  
  for (atlas in list_atlas_){
    path_file_fc<-file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep=""))
    if (!file.exists(path_file_fc)){
      print(paste("Atlas:",atlas,"FC data does not exist.",sep=" "))
    }else{
    
      # Load connection data
      dt_conn<-fread(path_file_fc)
      df_conn<-as.data.frame(dt_conn)
      dt_conn<-NULL
      gc()
      
      if (sum(df_conn$ses=="2-1")>0){
        print(paste("Atlas:",atlas,"FC difference already calculated.",sep=" "))
      }else{
        
        # Create dataframe of existing graph edges
        df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
        df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
        df_edge$from<-as.character(df_edge$from)
        df_edge$to<-as.character(df_edge$to)
        
        # Examine existing subject IDs and sessions in connection data
        df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
        colnames(df_ses_subj)<-c("ses","ID_pnTTC")
        for (ses in c(1,2)){
          df_ses_subj<-rbind(df_ses_subj,
                             data.frame(ses=ses,ID_pnTTC=sort(unique(df_conn[df_conn$ses==ses,"ID_pnTTC"]))))
        }
        list_subj_exist_long<-intersect(df_ses_subj[df_ses_subj$ses==1,"ID_pnTTC"],
                                        df_ses_subj[df_ses_subj$ses==2,"ID_pnTTC"])
        print(paste("Atlas: ",atlas,", ",
                    as.character(length(list_subj_exist_long))," subjects with longitudinal data.",
                    sep=""))
        
        # Calculate longitudinal fc difference
        df_out<-data.frame()
        for (id_subj in list_subj_exist_long){
          df_ses1<-df_conn[(df_conn$ses==1 & df_conn$ID_pnTTC==id_subj),c("from","to","r","z_r")]
          colnames(df_ses1)[colnames(df_ses1)=="r"]<-"ses1_r"
          colnames(df_ses1)[colnames(df_ses1)=="z_r"]<-"ses1_z_r"
          df_ses2<-df_conn[(df_conn$ses==2 & df_conn$ID_pnTTC==id_subj),c("from","to","r","z_r")]
          colnames(df_ses2)[colnames(df_ses2)=="r"]<-"ses2_r"
          colnames(df_ses2)[colnames(df_ses2)=="z_r"]<-"ses2_z_r"
          df_diff<-left_join(df_ses1,df_ses2,by=c("from","to"))
          df_diff$r<-df_diff$ses2_r-df_diff$ses1_r
          df_diff$z_r<-df_diff$ses2_z_r-df_diff$ses1_z_r
          df_diff<-data.frame(ses="2-1",ID_pnTTC=id_subj,df_diff[,c("from","to","r","z_r")])
          df_out<-rbind(df_out,df_diff)
        }
        df_out<-transform(df_out,ses=as.character(ses))
        df_out<-rbind.fill(df_conn,df_out)
        write.csv(df_out,file.path(paths_$output,"output",paste("atl-",atlas,"_fc.csv",sep="")),row.names=F)
      }
    }
  }
  print('Finished diff_fc()')
}


#**************************************************
# Component analyses of FC ========================
#**************************************************
# OBSOLETE
ca_fc<-function(paths_=paths,list_atlas_=list_atlas,list_wave_=list_wave,
                list_covar_=list_covar,subset_subj_=subset_subj,list_dim_ca_=list_dim_ca,
                plot_result=F){
  print("Starting ca_fc().")
  memory.limit(200000)
  nullobj<-func_createdirs(paths_,str_proc="ca_fc()",copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  for (atlas in list_atlas_){
    # Load and examine FC data
    print(paste("Loding FC of atlas: ",atlas,sep=""))
    
    # Create graph edge dataframe and node list
    #df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_conn<-fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
    n_edge<-dim(df_edge)[1]
    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
    n_node<-length(list_node)
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      id_subj_exist<-unique(df_conn_ses$ID_pnTTC)
      id_subj_subset<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
      id_subj_exist<-intersect(id_subj_exist,id_subj_subset)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist)
    }
    
    # Cbind FC data (Fisher-z transform of FC) as input for PCA function
    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
    df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
    colnames(df_clin_exist)<-colnames(df_clin)
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
        df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ses==ses & df_clin$ID_pnTTC==id_subj,])
      }
    }
    colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
    rownames(df_conn_cbind)<-NULL
    # Transpose connection dataframe (rows >> data for each subject/session, columns >> data for each edge)
    df_conn<-as.data.frame(t(df_conn_cbind))
    df_conn_cbind<-NULL
    gc()
    
    # Calculate PCA of FC
    print("Calculating PCA of FC.")
    dim_ca<-max(list_dim_ca_)
    data_pca<-func_pca(df_src=df_conn,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca)
    write.csv(data_pca$df_comp_mri,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variable.csv",sep="")),row.names=F)
    write.csv(data_pca$df_comp_subj,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_subject.csv",sep="")),row.names=F)
    write.csv(data_pca$df_variance,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variance.csv",sep="")),row.names=F)
    write.csv(data_pca$df_comp_clin_flat,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_correlation.csv",sep="")),row.names=F)
    
    # Plot PCA results
    if (plot_result){
      print("Plotting PCA of FC.")
      list_plot_pca<-plot_ca(df_src=data_pca$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_pca$dim)
      for (i_dim in names(list_plot_pca)){
        for (name_covar in names(list_plot_pca[[i_dim]])){
          plot<-list_plot_pca[[i_dim]][[name_covar]]
          plot<-(plot
                 + ggtitle("PCA of FC"))
          ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_pca.eps",sep=""),plot=plot,device=cairo_ps,
                 path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
        }
      }
    }
    data_pca<-NULL
    list_plot_pca<-NULL
    gc()
    
    # Calculate ICA of FC
    print("Calculating ICA of FC.")
    df_comp_mri_rbind<-data.frame()
    df_comp_subj_rbind<-data.frame()
    df_variance_rbind<-data.frame()
    df_comp_clin_flat_rbind<-data.frame()
    for (dim_ca in list_dim_ca_){
      data_ica<-func_ica(df_src=df_conn,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca)
      df_comp_mri_rbind<-rbind.fill(df_comp_mri_rbind,cbind(dim=dim_ca,data_ica$df_comp_mri))
      df_comp_subj_rbind<-rbind.fill(df_comp_subj_rbind,cbind(dim=dim_ca,data_ica$df_comp_subj))
      df_variance_rbind<-rbind.fill(df_variance_rbind,cbind(dim=dim_ca,data_ica$df_variance))
      df_comp_clin_flat_rbind<-rbind.fill(df_comp_clin_flat_rbind,cbind(dim=dim_ca,data_ica$df_comp_clin_flat))
    }
    write.csv(df_comp_mri_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_variable.csv",sep="")),row.names=F)
    write.csv(df_comp_subj_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_subject.csv",sep="")),row.names=F)
    write.csv(df_variance_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_variance.csv",sep="")),row.names=F)
    write.csv(df_comp_clin_flat_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_correlation.csv",sep="")),row.names=F)
    
    # Plot ICA results
    if (plot_result){
      print("Plotting ICA of FC.")
      list_plot_ica<-plot_ca(df_src=data_ica$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_ica$dim)
      for (i_dim in names(list_plot_ica)){
        for (name_covar in names(list_plot_ica[[i_dim]])){
          plot<-list_plot_ica[[i_dim]][[name_covar]]
          plot<-(plot
                 + ggtitle("ICA of FC"))
          ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_ica.eps",sep=""),plot=plot,device=cairo_ps,
                 path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
        }
      }
    }
    data_ica<-NULL
    list_plot_ica<-NULL
    gc()
  }
  print("Finished ca_fc().")
}


#**************************************************
# Fingerprinting ==================================
#**************************************************

# Core function for parallelization of fp_fc()
fp_fc_core<-function(data_zr){
  measure<-"fc"
  atlas<-data_zr$atlas
  paths<-data_zr$paths
  group_1<-data_zr$group[[1]]
  group_2<-data_zr$group[[2]]
  df_zr<-data_zr$df_zr
  df_ses_subj<-data_zr$df_ses_subj
  n_edge<-dim(df_zr)[1]
  
  # Calculate correlation matrix
  data_fingerprint<-func_cor(input=df_zr)
  df_fp_subnet<-data_fingerprint$cor_flat
  
  # Rename correlation matrix to sessions and subjects
  df_fp_subnet$from_ses<-df_fp_subnet$from_ID_pnTTC<-df_fp_subnet$to_ses<-df_fp_subnet$to_ID_pnTTC<-NA
  for (i in seq(dim(df_fp_subnet)[1])){
    from_id<-df_fp_subnet[[i,"from"]]
    to_id<-df_fp_subnet[[i,"to"]]
    df_fp_subnet[[i,"from_ses"]]<-df_ses_subj[[from_id,"ses"]]
    df_fp_subnet[[i,"from_ID_pnTTC"]]<-df_ses_subj[[from_id,"ID_pnTTC"]]
    df_fp_subnet[[i,"to_ses"]]<-df_ses_subj[[to_id,"ses"]]
    df_fp_subnet[[i,"to_ID_pnTTC"]]<-df_ses_subj[[to_id,"ID_pnTTC"]]
  }
  df_fp_subnet$measure<-measure
  df_fp_subnet$group_1<-group_1
  df_fp_subnet$group_2<-group_2
  df_fp_subnet<-df_fp_subnet[c("measure","group_1","group_2","from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r","z_r")]
  
  # rbind to output dataframe
  #df_fp<-rbind(df_fp,df_fp_subnet)
  
  # Prepare dataframe for fingerprint correlation plot
  df_fp_plot<-data_fingerprint$cor
  list_name_subj_ses<-paste(sprintf("%05d",df_ses_subj$ID_pnTTC),as.character(df_ses_subj$ses),sep="_")
  colnames(df_fp_plot)<-rownames(df_fp_plot)<-list_name_subj_ses
  
  # Heatmap plot of fp correlation matrix
  plot_fp_heatmap<-plot_cor_heatmap(input=df_fp_plot)
  suppressMessages(plot_fp_heatmap<-(plot_fp_heatmap
                                     + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                     + ggtitle(paste("FP Cor,",atlas,measure,group_1,group_2,sep=" "))
                                     + theme(plot.title = element_text(hjust = 0.5),
                                             axis.title=element_blank())))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp1-",group_1,"_grp2-",group_2,"_fp.eps",sep=""),plot=plot_fp_heatmap,device=cairo_ps,
         path=file.path(paths$output,"output","plot"),dpi=300,height=10,width=10,limitsize=F)
  
  return(df_fp_subnet)
}

# Main function for fingerprint computing

fp_fc<-function(paths_=paths,list_wave_=list_wave,list_atlas_=list_atlas,key_roigroup="group_3"){
  print("Starting fp_fc().")
  nullobj<-func_createdirs(paths_,str_proc="fp_fc()")
  dict_roi<-func_dict_roi(paths_)
  dict_roi<-data.frame(id=as.character(dict_roi$id),group=as.character(dict_roi[,key_roigroup]),stringsAsFactors = F)
  
  for (atlas in list_atlas_){
    if (file.exists(file.path(paths_$output,"output",paste("atl-",atlas,"_fp.csv",sep="")))){
      print(paste("Atlas: ",atlas,", fingerprint already calculated.",sep=""))
    }else{
      # Load connection data
      df_conn<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep=""))))
      df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
      df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
      df_edge$from<-as.character(df_edge$from)
      df_edge$to<-as.character(df_edge$to)
      
      # Examine existing subject IDs and sessions in connection data
      df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
      colnames(df_ses_subj)<-c("ses","ID_pnTTC")
      #list_ses_exist <- sort(unique(df_conn$ses))
      for (ses in list_wave_){
        df_ses_subj<-rbind(df_ses_subj,
                           data.frame(ses=ses,ID_pnTTC=sort(unique(df_conn[df_conn$ses==ses,"ID_pnTTC"]))))
      }
      
      # Add node subgroup column to df_edge
      df_edge<-left_join(df_edge,dict_roi,by=c("from"="id"))
      colnames(df_edge)[colnames(df_edge)=="group"]<-"from_group"
      df_edge<-left_join(df_edge,dict_roi,by=c("to"="id"))
      colnames(df_edge)[colnames(df_edge)=="group"]<-"to_group"
      
      # List groups of existing nodes
      list_group<-sort(unique(c(df_edge[,"from_group"],df_edge[,"to_group"])))
      if (!("whole" %in% list_group)){
        list_group<-c("whole",list_group)
      }
      n_group<-length(list_group)
      print(paste("Atlas: ",atlas, ", ", as.character(n_group)," groups:",sep=""))
      print(list_group)
      
      # Split and combine z_r data for each subgroup of networks for parallel computing
      list_data_zr<-list()
      
      for (idx_group_1 in seq(n_group)){
        for (idx_group_2 in seq(idx_group_1,n_group)){
          group_1<-list_group[idx_group_1]
          group_2<-list_group[idx_group_2]
          # 1. whole <-> whole
          # 2. (whole-group_2) <-> group_2 and group_2 <-> (whole-group_2)
          # 3. group_1 <-> group_2 and group_2 <-> group_1 (including group_1=group_2)
          if (group_1=="whole"){
            if (group_2=="whole"){
              # 1. whole <-> whole
              df_edge_group<-df_edge
            }else{
              # 2. (whole-group_2) <-> group_2 and group_2 <-> (whole-group_2)
              df_edge_group<-rbind(df_edge[df_edge$from_group!=group_2 & df_edge$to_group==group_2,],
                                   df_edge[df_edge$from_group==group_2 & df_edge$to_group!=group_2,])
            }
          }else{
            # 3. group_1 <-> group_2 and group_2 <-> group_1 (including group_1=group_2)
            if (group_1==group_2){
              df_edge_group<-df_edge[df_edge$from_group==group_1 & df_edge$to_group==group_1,]
            }else{
              df_edge_group<-rbind(df_edge[df_edge$from_group==group_1 & df_edge$to_group==group_2,],
                                   df_edge[df_edge$from_group==group_2 & df_edge$to_group==group_1,])
            }
          }
          
          df_edge_group$idx<-seq(nrow(df_edge_group))
          n_edge_group<-dim(df_edge_group)[1]
          
          if (n_edge_group<5){
            print(paste("Atlas: ",atlas,", Group: ",group_1," and ",group_2, ", Edges: ",as.character(n_edge_group)," < 5, fp calculation skipped.",sep=""))
          }else{
            # Create combined dataframe of Z-transformed correlation coefficients
            # according to pre-calculated edge and subject data
            print(paste("Atlas: ",atlas,", Group: ",group_1," and ",group_2,", preparing data.",sep=""))
            df_conn_cbind<-data.frame(matrix(nrow=n_edge_group,ncol=0))
            
            for (idx_subj in seq(nrow(df_ses_subj))){
              #print(paste("Ses: ",df_ses_subj[idx_subj,"ses"],", Subj: ",df_ses_subj[idx_subj,"ID_pnTTC"],sep=""))
              df_conn_subset<-df_conn[df_conn$ID_pnTTC==df_ses_subj[idx_subj,"ID_pnTTC"]
                                      & df_conn$ses==df_ses_subj[idx_subj,"ses"],c("from","to","z_r")]
              df_conn_subset$from<-as.character(df_conn_subset$from)
              df_conn_subset$to<-as.character(df_conn_subset$to)
              df_conn_subset<-inner_join(df_conn_subset,df_edge_group,by=c("from","to"))
              df_conn_subset<-df_conn_subset[order(df_conn_subset$idx),"z_r"]
              df_conn_cbind<-cbind(df_conn_cbind,df_conn_subset)
            }
            colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
            rownames(df_conn_cbind)<-NULL
            
            list_data_zr<-c(list_data_zr,list(list("group"=c(group_1,group_2),"atlas"=atlas,"paths"=paths_,
                                                   "df_zr"=df_conn_cbind,"df_ses_subj"=df_ses_subj)))
          }
        }
      }
      
      # Parallel fingerprint correlation computing over groups of subnetworks
      print(paste("Atlas: ",atlas,", calculating FC fingerprint correlation in parallel.",sep=""))
      n_cluster<-min(floor(detectCores()*3/4),length(list_data_zr))
      clust<-makeCluster(n_cluster)
      clusterExport(clust,
                    varlist=c("func_cor",
                              "plot_cor_heatmap","rcorr","func_fisherz","rownames_to_column","gather",
                              "ggplot","aes","geom_tile","scale_fill_gradientn",
                              "matlab.like2","scale_y_discrete","scale_x_discrete",
                              "theme_light","theme","element_text","element_blank",
                              "ggtitle","ggsave"),
                    envir=environment())
      list_df_fp<-pblapply(list_data_zr,fp_fc_core,cl=clust)
      stopCluster(clust)
      
      # Output dataframe
      df_fp<-NULL
      for (df_fp_subnet in list_df_fp){
        if (!is.null(df_fp_subnet)){
          df_fp<-rbind(df_fp,df_fp_subnet)
        }
      }
      
      # Save fingerprint correlation
      write.csv(df_fp,file.path(paths_$output,"output",paste("atl-",atlas,"_fp.csv",sep="")),row.names=F)
    }
  }
  print("Finished fp_fc().")
}


#**************************************************
# GTA functionalities =============================
#**************************************************

edges2igraph<-function(df_conn,df_edge,list_node,dict_roi){
  edges<-data.frame(matrix(ncol=3,nrow=dim(df_edge)[1]))
  edges[,1:2]<-df_edge[,c("from","to")]
  for (i in 1:dim(df_edge)[1]){
    edges[i,3]<-as.numeric(df_conn[intersect(which(df_conn$from==df_edge[i,"from"]),
                                             which(df_conn$to==df_edge[i,"to"])),"r"])
  }
  colnames(edges)<-c("from","to","weight")
  nodes<-data.frame(id=list_node)
  for (i in seq(length(list_node))){
    nodes[i,"label"]<-as.character(dict_roi[dict_roi$id==as.character(nodes[i,"id"]),"label"])
  }
  output <- graph.data.frame(d = edges, vertices = nodes, directed = F)
  return(output)
}


#**************************************************
# Binary GTA ======================================
#**************************************************

# Subset edges according to desired cost
subset_edge<-function(input_igraph, input_cost,n_node,n_edge){
  n_edges4cost<-as.integer(n_node*(n_node-1)/2*input_cost)
  edges2delete<-head(order(E(input_igraph)$weight),(n_edge-n_edges4cost))
  output<-delete.edges(input_igraph,edges2delete)
  return(output)
}

# Calculate binary graph metrics
gta_bin_metrics<-function(input_igraph){
  node<-names(V(input_igraph))
  metrics<-data.frame(matrix(nrow=0,ncol=3))
  colnames(metrics)<-c("node","metric","value")
  ## graph-level metrics
  # characteristic path length
  metrics<-rbind(metrics,cbind(node="graph",metric="characteristic path length",
                               value=average.path.length(input_igraph)))
  # global efficiency
  eff<-1/(shortest.paths(input_igraph))
  eff[!is.finite(eff)]<-0
  metrics<-rbind(metrics,cbind(node="graph",metric="global efficiency",
                               value=mean(eff,na.rm=TRUE)))
  # global clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="global clustering coefficient",
                               value=transitivity(input_igraph)))
  # average clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="average clustering coefficient",
                               value=transitivity(input_igraph,type="average")))
  # local efficiency
  # modularity
  # small-worldness
  suppressWarnings(metrics<-rbind(metrics,cbind(node="graph",metric="small-world index",
                                                value=smallworldIndex(input_igraph)$index)))
  
  ## node-level metrics
  # degree centrality
  metrics<-rbind(metrics,cbind(node=node,metric="degree centrality",
                               value=centr_degree(input_igraph)$res))
  # betweenness centrality
  metrics<-rbind(metrics,cbind(node=node,metric="betweenness centrality",
                               value=centr_betw(input_igraph)$res))
  # eigenvector centrality
  metrics<-rbind(metrics,cbind(node=node,metric="eigenvector centrality",
                               value=eigen_centrality(input_igraph)$vector))
  
  rownames(metrics)<-NULL
  return(metrics)
}


gta_bin<-function(paths_=paths,
                  list_atlas_=list_atlas,
                  list_cost_=list_cost){
  print("Starting to calculate binary GTA.")
  #data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  
  for (atlas in list_atlas_){
    print(paste("Calculate atlas: ",atlas,sep=""))
    file_conn<-paste("atl-",atlas,"_fc.csv",sep="")
    df_conn<-read.csv(file.path(paths_$input,"output",file_conn))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    n_edge<-dim(df_edge)[1]
    list_node<-unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to))))
    list_node<-list_node[order(list_node)]
    n_node<-length(list_node)
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_conn_ses$ID_pnTTC))
    }
    #df_dst<-data.frame()
    list_file_tmp<-NULL
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
      #for (id_subj in list_id_subj_exist[[ses]][c(1,2)]){
        print(paste("Calculating Wave: ",as.character(ses), ", Subject: ",as.character(id_subj),sep=""))
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        igraph_subj<-edges2igraph(df_conn=df_conn_subj,df_edge=df_edge,list_node=list_node,dict_roi=dict_roi)
        
        df_metric_subj<-data.frame()
        for (cost in list_cost_){
          igraph_subj_subset<-subset_edge(igraph_subj,cost,n_node,n_edge)
          E(igraph_subj_subset)$weight<-1
          metrics<-gta_bin_metrics(igraph_subj_subset)
          metrics<-cbind(cost=cost,metrics)
          df_metric_subj<-rbind(df_metric_subj,metrics)
        }
        df_metric_subj$value<-as.numeric.factor(df_metric_subj$value)
        df_metric<-df_metric_subj[which(df_metric_subj$cost==list_cost_[1]),c("node","metric")]
        average<-data.frame()
        for (i in seq(dim(df_metric)[1])){
          average<-rbind(average,
                         cbind(cost="average",node=as.character(df_metric[i,"node"]),
                               metric=as.character(df_metric[i,"metric"]),
                               value=mean(df_metric_subj[intersect(which(df_metric_subj$node==df_metric[i,"node"]),
                                                                   which(df_metric_subj$metric==df_metric[i,"metric"])),"value"])))
        }
        df_metric_subj<-rbind(df_metric_subj,average)
        rownames(df_metric_subj)<-NULL
        df_metric_subj<-cbind(ses=ses,ID_pnTTC=id_subj,df_metric_subj)
        file_metric_tmp<-paste("TMP_atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d",id_subj),"_gta_bin.csv",sep="")
        path_file_metric_tmp<-file.path(paths_$output,"output",file_metric_tmp)
        write.csv(df_metric_subj,path_file_metric_tmp,row.names=F)
        list_file_tmp<-c(list_file_tmp,path_file_metric_tmp)
        #df_dst<-rbind(df_dst,df_metric_subj)
      }
    }
    df_dst<-data.frame()
    for (path_file_metric_tmp in list_file_tmp){
      df_tmp<-read.csv(path_file_metric_tmp)
      df_dst<-rbind(df_dst,df_tmp)
      file.remove(path_file_metric_tmp)
      print(paste("Finished binding: ",path_file_metric_tmp,sep=""))
    }
    file_dst<-paste("atl-",atlas,"_gta_bin.csv",sep="")
    write.csv(df_dst,file.path(paths_$output,"output",file_dst),row.names=F)
  }
  print("Finished gta_bin().")
}


#**************************************************
# Weighted GTA ====================================
#**************************************************

# add metric to output list in weighted GTA
AddMetric<-function(input){
  output<-data.frame(matrix(nrow=0,ncol=3))
  if (!is.null(input$graph)){
    output_add<-cbind(node="graph",metric=input$name[[1]],value=input$graph)
    output<-rbind(output,output_add)
  }
  if (!is.null(input$node)){
    output_add<-cbind(node=names(input$node),
                      metric=input$name[[1]],value=input$node)
    output<-rbind(output,output_add)
  }
  colnames(output)<-c("node","metric","value")
  return(output)
}

WeightedMetric<-function(input_igraph){
  metrics<-data.frame(matrix(nrow=0,ncol=3))
  distance<-WeightedDistance(input_igraph)$distance
  
  metrics<-rbind(metrics,AddMetric(WeightedCharPath(input_distance=distance)))
  #metrics<-rbind(metrics,AddMetric(WeightedEccentricity(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedRadius(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedDiameter(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedGlobalEfficiency(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedClustCoef(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedTransitivity(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedLocalEfficiency(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedModularity(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedStrength(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedClosenessCentrality(input_distance = distance)))
  #metrics<-rbind(metrics,AddMetric(WeightedBetweennessCentrality(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedEigenvectorCentrality(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedNeighborDegree(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedAssortativityCoef(input = input_igraph)))
  
  colnames(metrics)<-c("node","metric","value")
  rownames(metrics)<-NULL
  return(metrics)
}


gta_weight<-function(absolute=T,
                     threshold=NA,
                     paths_=paths,
                     list_atlas_=list_atlas,
                     list_wave_=list_wave,
                     subset_subj_=subset_subj){
  print("Starting gta_weight().")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  df_clin<-func_clinical_data_long(paths_,list_wave_)
  data_subset_clin<-func_subset_clin(df_clin,
                                     list_wave_,subset_subj_,
                                     list_covar=NULL,
                                     rem_na_clin=T)
  df_clin_subset<-data_subset_clin$df_clin

  for (atlas in list_atlas_){
    print(paste("Calculate atlas: ",atlas,sep=""))
    file_conn<-paste("atl-",atlas,"_fc.csv",sep="")
    df_conn<-read.csv(file.path(paths_$input,"output",file_conn))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    n_edge<-dim(df_edge)[1]
    list_node<-unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to))))
    list_node<-list_node[order(list_node)]
    n_node<-length(list_node)
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      id_subj_exist<-unique(df_conn_ses$ID_pnTTC)
      id_subj_subset<-df_clin_subset[df_clin_subset$wave==ses,"ID_pnTTC"]
      id_subj_exist<-intersect(id_subj_exist,id_subj_subset)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist)
    }
    #df_dst<-data.frame()
    list_file_tmp<-NULL
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        #for (id_subj in list_id_subj_exist[[ses]][c(1,2)]){
        print(paste("Calculating Wave: ",as.character(ses), ", Subject: ",as.character(id_subj),sep=""))
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        igraph_subj<-edges2igraph(df_conn=df_conn_subj,df_edge=df_edge,list_node=list_node,dict_roi=dict_roi)
        if (absolute){
          E(igraph_subj)$weight<-abs(E(igraph_subj)$weight)
        }
        E(igraph_subj)$weight[is.na(E(igraph_subj)$weight)]<-0
        df_metric_subj<-WeightedMetric(igraph_subj)
        df_metric_subj$value<-as.numeric.factor(df_metric_subj$value)
        rownames(df_metric_subj)<-NULL
        df_metric_subj<-cbind(ses=ses,ID_pnTTC=id_subj,df_metric_subj)
        file_metric_tmp<-paste("TMP_atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d",id_subj),"_gta_weight.csv",sep="")
        path_file_metric_tmp<-file.path(paths_$output,"output",file_metric_tmp)
        write.csv(df_metric_subj,path_file_metric_tmp,row.names=F)
        list_file_tmp<-c(list_file_tmp,path_file_metric_tmp)
        #df_dst<-rbind(df_dst,df_metric_subj)
      }
    }
    df_dst<-data.frame()
    for (path_file_metric_tmp in list_file_tmp){
      df_tmp<-read.csv(path_file_metric_tmp)
      df_dst<-rbind(df_dst,df_tmp)
      file.remove(path_file_metric_tmp)
      print(paste("Finished binding: ",path_file_metric_tmp,sep=""))
    }
    file_dst<-paste("atl-",atlas,"_gta_weight.csv",sep="")
    write.csv(df_dst,file.path(paths_$output,"output",file_dst),row.names=F)
  }
  print("Finished gta_weight().")
}


#**************************************************
# FC-FC correlation ===============================
#**************************************************
# for comparison of preprocessing methods

fc_corr<-function(paths_=paths,subset_subj_=subset_subj){
  print("Starting to calculate FC-FC correlation.")
  df_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=F)
  figs<-list()
  for (id_subj in df_clinical$list_id_subj){
    print(paste("Starting to calculate",as.character(id_subj),sep=" "))
    list_path_file_input<-NULL
    id_study_first<-0
    for (id_study in seq(length(paths_$dir_in))){
      file_input<-paste("fc_",sprintf("%05d", id_subj),"_rp.csv",sep="")
      path_file_input<-file.path(paths_$input[id_study],"output",file_input)
      if (file.exists(path_file_input)){
        list_path_file_input<-c(list_path_file_input,path_file_input)
        if(id_study_first==0){
          id_study_first<-id_study
        }
      }else{
        list_path_file_input<-c(list_path_file_input,NA)
      }
    }
    if(id_study_first==0){
      print("No input available.")
    }
    
    for (id_study in seq(length(paths_$dir_in))){
      path_file_input<-list_path_file_input[id_study]
      if(!is.na(path_file_input)){
        df_fc<-read.csv(path_file_input)
        if(id_study==id_study_first){
          df_fc_allstudy<-data.frame(matrix(ncol=length(paths_$dir_in)+2,nrow=nrow(df_fc)))
          colnames(df_fc_allstudy)<-c("from","to",paths_$dir_in)
        }
        df_fc_allstudy[,c("from","to",paths_$dir_in[id_study])]<-df_fc[,c("from","to","r")]
      }
    }
    
    fig<-ggpairs(df_fc_allstudy[,c(-1,-2)],
                 upper=list(continuous=custom_corr_heatmap),
                 #lower=list(continuous=wrap("points",alpha=0.01,size=0.001,stroke = 0, shape = ".")),
                 lower=list(continuous=custom_smooth),
                 diag=list(continuous=custom_densityDiag),
                 title=paste(sprintf("%05d",id_subj),"fc_corr",sep="_"))
    ggsave(paste(sprintf("%05d",id_subj),"fc_corr.eps",sep="_"),plot=fig,device=cairo_ps,
           path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    figs<-c(figs,list(fig))
    print(paste("Finished calculating subject",as.character(id_subj),sep=" "))
  }
  names(figs)<-as.character(df_clinical$list_id_subj)
  print("Finished calculating FC-FC correlation.")
  return(figs)
}

