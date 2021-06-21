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
#path_exp_full<-"/media/atiroms/SSD_03/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in<-"411_fc_acompcor_gsr"
#dir_out<-"414_fc_gamm_acompcor_gsr_test2" 

dir_in<-"421_fc_aroma"
dir_out<-"423.3_fc_gam_diff_aroma_test12" 
#dir_out<-"422.1_fp_aroma_test2"
#dir_out<-"425.1_fc_ca_aroma_test2"
#dir_out<-"424_fc_gamm_aroma_test31" 


#**************************************************
# Libraries =======================================
#**************************************************
library(easypackages)
#libraries(ggplot2,GGally,igraph,ggrepel,colorRamps,tidyverse,parallel,mgcv,car,plyr,dplyr,data.table,pbapply,stringr)
libraries("ggplot2","colorRamps","tidyverse","parallel","mgcv","dplyr","data.table","pbapply","stringr","lmerTest")


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/gta_function.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# GLM/GAM of FC longitudinal difference ===========
#**************************************************

gam_fc_diff_core<-function(paths,data_fc,atlas,param,list_covar,list_mod,list_term,idx_var,
                           calc_parallel=T,test_mod=F
){
  # Prepare clinical data and demean
  data_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,prefix=paste("var-",idx_var,sep=""),print_terminal=F)
  list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
  df_clin<-data_clin$df_clin
  df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]                         # select longitudinal data
  df_clin<-func_omit_decreasing(df_clin,var_check=param$omit_decreasing)        # omit decreasing data in var_check
  df_clin<-dplyr::rename(df_clin,c("ses"="wave"))
  df_clin<-func_clinical_data_diffmean(df_clin,list_id_subj,list_covar)         # calculate longitudinal difference and mean
  df_clin<-data.frame(wave="2-1",func_std_clin(df_clin,separate_sex=T)$df_clin) # standardize
  fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_clin.csv",sep="")),row.names=F)
  
  label_wave<-"2-1"
  # Calculate model
  data_gamm<-func_calc_gamm(paths,param,data_fc,df_clin,atlas,list_covar,list_mod,list_term,idx_var,label_wave,calc_parallel,test_mod)
  
  if (param$param_nbs$tfnbs){
    if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
    #if (calc_parallel){clust<-makeCluster(floor(detectCores()*1/2))}else{clust<-makeCluster(1)}
    list_sex<-param$list_sex
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                                  "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC",
                                  "param","func_bfs","%nin%","left_join","ggsave"),
                  envir=environment())
    
    # Calculate threshold-free network based statistics
    data_tfnbs<-func_iterate_tfnbs2(paths,clust,df_gamm,df_anova,df_deltah_in=NULL,var_exp_perm_in=NULL,data_fc,plot_result=T,return_nbs=T,
                                    atlas,param,list_mod,list_term,idx_var,label_wave)
    # Permutation test
    func_tfnbs_permutation(paths,clust,data_fc,df_clin,data_tfnbs_in=data_tfnbs,calc_parallel,plot_result=T,
                           atlas,param,list_mod,list_term,idx_var,label_wave)
    stopCluster(clust)
    
  }else{
    # Threshold and plot graph edges
    data_plot<-func_threshold_gamm(paths,param,data_gamm,data_fc,atlas,list_covar,list_mod,list_term,idx_var,label_wave)
    # Detect sub-network by breadth-first approach
    data_bfs<-func_detect_subnet(paths,param,data_plot,data_gamm,data_fc,atlas,list_covar,list_mod,list_term,idx_var,label_wave,plot_result=F)
    # Permutation test
    data_nbs<-func_nbs_permutation(paths,param,df_clin,data_bfs,data_fc, atlas,list_covar,list_mod,list_term,idx_var,label_wave,calc_parallel,plot_result=T)
  }
}


gam_fc_diff<-function(param=param_gam_fc_diff){
  print("Starting gam_fc_diff().")
  nullobj<-func_createdirs(paths,str_proc="gam_fc_diff()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in param$list_atlas){
    print(paste("Preparing FC data: ",atlas,sep=""))
    data_fc<-prep_data_fc2(paths,atlas,param$key_group,list_wave="2-1",include_grp=T,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    fwrite(data_fc$df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_src_fc.csv",sep="")),row.names=F)
    fwrite(data_fc$df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_src_fc_grp.csv",sep="")),row.names=F)
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_var in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_var]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_var]]
      list_covar$tanner$dtype<-param$dtype_tanner
      gam_fc_diff_core(paths,data_fc,atlas,param,list_covar,
                       list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_var)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_var in names(param$list_hormone)){
      print(paste("Atlas: ",atlas,", Hormone type: ",param$list_hormone[[idx_var]][["label"]],sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_var]]
      gam_fc_diff_core(paths,data_fc,atlas,param,list_covar,
                       list_mod=param$list_mod_hormone,list_term=param$list_term_hormone,idx_var=idx_var)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  list_var<-c(param$list_tanner,param$list_hormone)
  func_combine_result(paths,param$list_atlas,list_var,"2-1",list(list("measure"="")),c("gamm","plot","gamm_anova","gamm_aic","gamm_grp","plot_grp","gamm_anova_grp","gamm_aic_grp","bfs_edge","bfs_node","bfs_size","bfs_pred","perm_max","perm_thr","perm_fwep","perm_sign","tfnbs"))
  
  print("Finished gam_fc_diff().")
}


#**************************************************
# Fingerprinting ==================================
#**************************************************

# Core function for parallelization of fp_fc()
fp_fc_core<-function(data_zr){
  atlas<-data_zr$atlas
  paths<-data_zr$paths
  group_1<-data_zr$group[[1]]
  group_2<-data_zr$group[[2]]
  df_zr<-data_zr$df_zr
  df_ses_subj<-data_zr$df_ses_subj
  n_edge<-dim(df_zr)[1]
  
  # Calculate correlation matrix
  data_fp<-func_cor(input=df_zr)
  df_fp<-data_fp$cor_flat
  
  # Rename correlation matrix to sessions and subjects
  df_ses_subj$id<-as.character(seq(nrow(df_ses_subj)))
  df_fp<-inner_join(df_fp,df_ses_subj,by=c("from"="id"))
  df_fp<-rename(df_fp,"from_ses"="ses","from_ID_pnTTC"="ID_pnTTC")
  df_fp<-inner_join(df_fp,df_ses_subj,by=c("to"="id"))
  df_fp<-rename(df_fp,"to_ses"="ses","to_ID_pnTTC"="ID_pnTTC")
  df_fp$group_1<-group_1
  df_fp$group_2<-group_2
  df_fp<-df_fp[c("group_1","group_2","from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r","z_r")]
  
  for (id_row in seq(nrow(df_ses_subj))){
    ses<-df_ses_subj[id_row,"ses"]
    id_subj<-df_ses_subj[id_row,"ID_pnTTC"]
    df_fp_subset<-rbind(df_fp[df_fp$from_ses==ses & df_fp$from_ID_pnTTC==id_subj,],
                        df_fp[df_fp$to_ses==ses & df_fp$to_ID_pnTTC==id_subj,])
    mean_zr<-mean(df_fp_subset$z_r,na.rm=T)
    df_ses_subj[id_row,"mean_z_r"]<-mean_zr
  }
  
  # Prepare dataframe for fingerprint correlation plot
  df_fp_plot<-data_fp$cor
  list_name_subj_ses<-paste(sprintf("%05d",df_ses_subj$ID_pnTTC),as.character(df_ses_subj$ses),sep="_")
  colnames(df_fp_plot)<-rownames(df_fp_plot)<-list_name_subj_ses
  
  # Heatmap plot of fp correlation matrix
  plot_fp_heatmap<-plot_cor_heatmap(input=df_fp_plot)
  suppressMessages(plot_fp_heatmap<-(plot_fp_heatmap
                                     + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                     + ggtitle(paste("FP Cor,",atlas,group_1,group_2,sep=" "))
                                     + theme(plot.title = element_text(hjust = 0.5),
                                             axis.title=element_blank())))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_grp1-",group_1,"_grp2-",group_2,"_fp.eps",sep=""),plot=plot_fp_heatmap,device=cairo_ps,
         path=file.path(paths$output,"output","plot"),dpi=300,height=10,width=10,limitsize=F)
  
  df_mean_fp<-data.frame("group_1"=group_1,"group_2"=group_2,df_ses_subj[,c("ses","ID_pnTTC","mean_z_r")])
  
  return(list("df_fp"=df_fp,"df_mean_fp"=df_mean_fp))
}

# Main function for fingerprint computing

#fp_fc<-function(paths_=paths,list_wave_=list_wave,list_atlas_=list_atlas,key_roigroup="group_3"){
fp_fc<-function(paths_=paths,list_atlas_=list_atlas,param=param_fp_fc){
  print("Starting fp_fc().")
  nullobj<-func_createdirs(paths_,str_proc="fp_fc()",list_param=param)
  dict_roi<-func_dict_roi(paths_)
  dict_roi<-data.frame(id=as.character(dict_roi$id),group=as.character(dict_roi[,param$key_group]),stringsAsFactors = F)
  
  df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar=NULL,rem_na_clin=T,
                                   prefix="",print_terminal=F)$df_clin
  df_clin<-rename(df_clin,"ses"="wave")
  df_clin$qc<-1
  
  for (atlas in list_atlas_){
    if (file.exists(file.path(paths_$output,"output",paste("atl-",atlas,"_fp.csv",sep="")))){
      print(paste("Atlas: ",atlas,", fingerprint already calculated.",sep=""))
    }else{
      # Load connection data
      data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=param$list_wave,include_grp=param$calc_group,
                             abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
      df_conn<-data_fc$df_fc; df_edge<-data_fc$df_edge
      #df_conn[is.na(df_conn$z_r),"z_r"]<-0
      
      # Examine existing subject IDs and sessions in connection data
      df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
      colnames(df_ses_subj)<-c("ses","ID_pnTTC")
      #list_ses_exist <- sort(unique(df_conn$ses))
      for (ses in param$list_wave){
        df_ses_subj<-rbind(df_ses_subj,
                           data.frame(ses=ses,ID_pnTTC=sort(unique(df_conn[df_conn$ses==ses,"ID_pnTTC"]))))
      }
      
      # Add node subgroup column to df_edge
      df_edge<-left_join(df_edge,dict_roi,by=c("from"="id"))
      colnames(df_edge)[colnames(df_edge)=="group"]<-"from_group"
      df_edge<-left_join(df_edge,dict_roi,by=c("to"="id"))
      colnames(df_edge)[colnames(df_edge)=="group"]<-"to_group"
      
      if (param$calc_group){
        # List groups of existing nodes
        list_group<-sort(unique(c(df_edge[,"from_group"],df_edge[,"to_group"])))
        if (!("whole" %in% list_group)){
          list_group<-c("whole",list_group)
        }
        n_group<-length(list_group)
        print(paste("Atlas: ",atlas, ", ", as.character(n_group)," groups:",sep=""))
        print(list_group)
      }else{
        list_group<-"whole"
        n_group<-1
      }
      
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
                              "inner_join","rename",
                              "ggplot","aes","geom_tile","scale_fill_gradientn",
                              "matlab.like2","scale_y_discrete","scale_x_discrete",
                              "theme_light","theme_linedraw","theme","element_text","element_blank",
                              "ggtitle","ggsave"),
                    envir=environment())
      list_data_fp<-pblapply(list_data_zr,fp_fc_core,cl=clust)
      stopCluster(clust)
      
      # Output dataframe
      df_fp<-df_mean_fp<-NULL
      for (data_fp in list_data_fp){
        if (!is.null(data_fp)){
          df_fp<-rbind(df_fp,data_fp$df_fp)
          df_mean_fp<-rbind(df_mean_fp,data_fp$df_mean_fp)
        }
      }
      
      df_mean_fp<-left_join(df_mean_fp,df_clin,by=c("ses","ID_pnTTC"))
      df_mean_fp[is.na(df_mean_fp$qc),"qc"]<-0
      
      # Save fingerprint correlation
      fwrite(df_fp,file.path(paths_$output,"output","result",paste("atl-",atlas,"_fp.csv",sep="")),row.names=F)
      fwrite(df_mean_fp,file.path(paths_$output,"output","result",paste("atl-",atlas,"_mean_fp.csv",sep="")),row.names=F)
    }
  }
  print("Finished fp_fc().")
}



#**************************************************
# PCA of FC =======================================
#**************************************************
ca_fc<-function(paths_=paths,list_atlas_=list_atlas,param=param_ca_fc){
  print("Starting ca_fc()")
  nullobj<-func_createdirs(paths_,str_proc="ca_fc()",copy_log=T,list_param=param)
  # Increase memory limit for later ICA calculation
  memory.limit(1000000)
  for (atlas in list_atlas_){
    # Load and examine FC data
    data_fc<-prep_data_fc2(paths,atlas,param$key_group,list_wave=param$list_wave_mri,include_grp=F,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    df_fc<-data_fc$df_fc; df_edge<-data_fc$df_edge
    
    # Calculate PCA/ICA factors (if not yet)
    for (wave_mri in param$list_wave_mri){
      path_pca_subj<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep=""))
      if (!(file.exists(path_pca_subj))){
        df_pca_mri<-df_pca_subj<-df_pca_vaf<-NULL
        for (label_sex in names(param$list_sex)){
          print(paste("Calculating Atlas: ",atlas,", MRI wave: ",wave_mri,", Sex: ",label_sex,sep=""))
          
          # Prepare subject subsetting condition (MRI QC criteria and sex) according to specified MRI wave
          # Existence of clinical variables are not considered here
          subset_subj_temp<-list(c(param$subset_subj[[wave_mri]],list(list("key"="Sex","condition"=param$list_sex[[label_sex]]))))
          names(subset_subj_temp)<-wave_mri
          df_clin<-func_clinical_data_long(paths_,wave_mri,subset_subj_temp,list_covar=NULL,rem_na_clin=F,
                                           prefix=paste("wav-m",wave_mri,"_sex-",label_sex,sep=""),print_terminal=F)$df_clin
          df_clin<-dplyr::rename(df_clin,"ses"="wave")
          
          df_fc_ses<-df_fc[df_fc$ses==wave_mri,]
          df_fc_calc<-df_clin_exist<-NULL
          for (id_subj in df_clin$ID_pnTTC){
            df_fc_subj<-df_fc_ses[df_fc_ses$ID_pnTTC==id_subj,"z_r"]
            df_fc_calc<-rbind(df_fc_calc,df_fc_subj)
            df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ID_pnTTC==id_subj,])
          }
          rownames(df_fc_calc)<-as.character(seq(nrow(df_fc_calc))); colnames(df_fc_calc)<-as.character(seq(ncol(df_fc_calc)))
          colnames(df_clin_exist)<-colnames(df_clin)
          df_clin_exist$ses<-NULL
          gc()
          
          # Calculate PCA of FC
          data_pca<-func_pca(df_src=df_fc_calc,df_var=data_fc$df_edge,df_indiv=df_clin_exist,dim_ca=param$dim_ca,calc_corr=F)
          
          # Save results
          df_pca_mri<-rbind(df_pca_mri,cbind(sex=label_sex,dim=param$dim_ca,data_pca$df_comp_mri))
          df_pca_subj<-rbind(df_pca_subj,cbind(sex=label_sex,dim=param$dim_ca,data_pca$df_comp_subj))
          df_pca_vaf<-rbind(df_pca_vaf,cbind(sex=label_sex,dim=param$dim_ca,data_pca$df_vaf))
          data_pca<-NULL
          gc()
          
        } # end for over sex
        fwrite(df_pca_mri,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep="")),row.names=F)
        fwrite(df_pca_subj,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep="")),row.names=F)
        fwrite(df_pca_vaf,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_vaf.csv",sep="")),row.names=F)
      } # end if not PCA/ICA factors are already calculated
    } # end for waves
    # end calculating PCA/ICA factors
    
    # Group-wise average of factor-MRI matrix
    print(paste("Calculating group-wise contribution to factors:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      path_pca_mri_grp<-file.path(paths_$output,"output","temp", paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var_grp.csv",sep=""))
      if (!(file.exists(path_pca_mri_grp))){
        df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
        df_pca_mri_grp<-NULL
        # PCA
        df_pca_mri_grp<-group_factor(df_pca_mri,param$dim_ca,data_fc$df_roi,data_fc$df_grp,param$list_sex)
        fwrite(df_pca_mri_grp,path_pca_mri_grp,row.names=F)
      }
    }
    # end calculating group-wise average
    
    # Generate visual output of MRI factors
    print(paste("Generating heatmap plot of PCA/ICA factors:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
      df_pca_mri_grp<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var_grp.csv",sep=""))))
      # Visual output of PCA/ICA factors
      for (label_sex in names(param$list_sex)){
        # PCA
        dim_ca<-param$dim_ca
        df_pca_mri_subset<-df_pca_mri[df_pca_mri$sex==label_sex & df_pca_mri$dim==param$dim_ca,]
        df_pca_mri_grp_subset<-df_pca_mri_grp[df_pca_mri_grp$sex==label_sex & df_pca_mri_grp$dim==param$dim_ca,]
        df_pca_mri_subset$sex<-df_pca_mri_subset$dim<-df_pca_mri_grp_subset$sex<-df_pca_mri_grp_subset$dim<-NULL
        # Visualize factor-FC matrix in heatmap plot
        plot_ca_fc_heatmap(paths_=paths_,data_fc,df_pca_mri_subset,df_pca_mri_grp_subset,atlas=atlas,dim_ca=param$dim_ca,
                           method="pca",label_sex=label_sex,ses=wave_mri)
      }
    } # End for MRI wave
    # End of generating factor heatmap plot
    
    # Calculate node strength in each factor
    print(paste("Calculating node strength of PCA/ICA factors:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      path_file_check<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_str.csv",sep=""))
      if (!file.exists(path_file_check)){
        df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
        df_strength<-NULL
        list_ca_mri<-list("pca"=df_pca_mri)
        for (method in names(list_ca_mri)){
          df_ca_mri<-list_ca_mri[[method]]
          list_sex<-unique(df_ca_mri$sex)
          list_dim<-unique(df_ca_mri$dim)
          for (sex in list_sex){
            for (dim in list_dim){
              df_ca_mri_subset<-df_ca_mri[df_ca_mri$sex==sex & df_ca_mri$dim==dim,]
              for (idx_row in seq(nrow(data_fc$df_roi))){
                node<-data_fc$df_roi[idx_row,"id"]
                label_node<-data_fc$df_roi[idx_row,"label"]
                strength<-colSums(abs(df_ca_mri_subset[df_ca_mri_subset$from==node | df_ca_mri_subset$to==node,
                                                       sprintf("comp_%03d",seq(dim))]))
                df_strength_add<-data.frame(method=method,sex=sex,dim=dim,comp=seq(dim),node=node,label_node=label_node,strength=strength)
                df_strength<-rbind(df_strength,df_strength_add)
              }
            }
          }
        }
        fwrite(df_strength,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_str.csv",sep="")),row.names=F)
      }
    }
    # End calculating node strength in each factor
    
    # Calculate factor attribution-clinical GLM/GAM
    print(paste("Calculating factor-clinical GLM/GAM:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      if (!file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_aic.csv",sep="")))){
        df_pca_subj<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep=""))))
        df_gamm<-df_anova<-df_aic<-NULL
        
        for (idx_var in names(param$list_tanner)){
          list_covar<-param$list_covar_tanner
          list_covar[["tanner"]]<-param$list_tanner[[idx_var]]
          
          # Prepare clinical data
          list_df_clin<-list()
          # clinical wave = 1 and 2 (longitudinal difference and mean)
          data_clin<-func_clinical_data_long(paths,c(1,2),param$subset_subj,list_covar,rem_na_clin=T,prefix=paste("var-",idx_var,sep=""),print_terminal=F)
          list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
          df_clin<-data_clin$df_clin
          df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
          colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
          df_clin<-func_clinical_data_diffmean(df_clin,list_id_subj,list_covar)
          df_clin<-func_std_clin(df_clin,separate_sex=T)$df_clin
          list_df_clin<-c(list_df_clin,list("long"=df_clin))
          
          # clinical wave = 1 or 2
          subset_subj_temp<-param$subset_subj[wave_mri]
          for (wave_clin in c(1,2)){
            names(subset_subj_temp)<-as.character(wave_clin)
            data_clin<-func_clinical_data_long(paths,wave_clin,subset_subj_temp,list_covar,rem_na_clin=T,prefix=paste("var-",idx_var,sep=""),print_terminal=F)
            df_clin<-func_std_clin(data_clin$df_clin,separate_sex=T)$df_clin
            df_clin$wave<-NULL
            list_df_clin<-c(list_df_clin,list(df_clin))
            names(list_df_clin)[length(list_df_clin)]<-paste("cs",as.character(wave_clin),sep='')
          }
          
          # Calculate GAMM
          for (label_df_clin in names(list_df_clin)){
            df_clin<-list_df_clin[[label_df_clin]]
            for (label_sex in names(param$list_sex)){
              list_sex<-c(1,2)
              eval(parse(text=paste('list_sex<-list_sex[list_sex',param$list_sex[[label_sex]],']',sep='')))
              list_sex<-list(list_sex)
              for (id_comp in seq(param$dim_ca)){
                df_pca_subset<-df_pca_subj[df_pca_subj$sex==label_sex,c("ID_pnTTC",sprintf("comp_%03d",id_comp))]
                colnames(df_pca_subset)<-c("ID_pnTTC","value")
                df_join<-dplyr::inner_join(df_pca_subset,df_clin,by="ID_pnTTC")
                data_gamm<-gamm_core4(df_join,list_mod_in=param$list_mod_tanner[[label_df_clin]],list_sex_in=list_sex,calc_parallel_in=F,test_mod_in=F)
                df_gamm_add<-data_gamm$df_gamm; df_anova_add<-data_gamm$df_anova; df_aic_add<-data_gamm$df_aic
                df_gamm_add$sex<-df_anova_add$sex<-df_aic_add$sex<-NULL
                df_head<-data.frame(wave_clin=label_df_clin,variable=idx_var,method="pca",sex=label_sex,dim=param$dim_ca,comp=id_comp)
                df_gamm<-rbind(df_gamm,cbind(df_head,df_gamm_add)); df_anova<-rbind(df_anova,cbind(df_head,df_anova_add)); df_aic<-rbind(df_aic,cbind(df_head,df_aic_add))
              }
            }
          }
        } # End for Tanenr type
        
        fwrite(df_gamm,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_gamm.csv",sep="")),row.names=F)
        fwrite(df_anova,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_anova.csv",sep="")),row.names=F)
        fwrite(df_aic,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_aic.csv",sep="")),row.names=F)
      } # End if calculated result does not exist
    } # End for MRI wave
  
  } # End for atlas
  
  # Reload and bind all results
  print("Binding results.")
  func_combine_result(paths_,list_atlas_,list_var=list("null"=NULL),list_wave=paste("m",param$list_wave_mri,sep=""),list_type_measure=list(list("measure"="")),
                      list_filename=c("fc_pca_subj","fc_ica_subj","fc_pca_var","fc_ica_var","fc_pca_vaf","fc_ica_vaf",
                                      "fc_pca_var_grp","fc_ica_var_grp","fc_ca_cor","fc_ca_str",
                                      "fc_ca_gamm","fc_ca_anova","fc_ca_aic"))
  print("Finished ca_fc_cs()")
}


#**************************************************
# FC mean and SD calculation ======================
#**************************************************
mean_fc<-function(paths_=paths,list_atlas_=list_atlas,param=param_gamm_fc){
  print("Starting mean_fc().")
  nullobj<-func_createdirs(paths_,str_proc="mean_fc()",copy_log=T,list_param=param)
  list_dir_src<-c("401_fc_acompcor","411_fc_acompcor_gsr","421_fc_aroma","431_fc_aroma_gsr")
  for (dir_src in list_dir_src){
    paths<-func_path(path_exp_=path_exp,dir_in_=dir_src,dir_out_=dir_out,path_exp_full_=path_exp_full)
    for (atlas in list_atlas_){
      data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=c("1","2"),include_grp=F,
                             abs_nfc=F,std_fc=F,div_mean_fc=F)
      df_fc<-data_fc$df_fc
      df_mean_sd<-data.frame()
      list_wave_exist<-sort(unique(df_fc$ses))
      for (wave in list_wave_exist){
        list_id_subj<-sort(unique(df_fc[df_fc$ses==wave,"ID_pnTTC"]))
        for (id_subj in list_id_subj){
          df_fc_subset<-df_fc[df_fc$ses==wave & df_fc$ID_pnTTC==id_subj,]
          mean_fc<-mean(df_fc_subset$z_r)
          sd_fc<-sd(df_fc_subset$z_r)
          df_mean_sd<-rbind(df_mean_sd,data.frame("wave"=wave,"ID_pnTTC"=id_subj,"mean"=mean_fc,"sd"=sd_fc))
        }
      }
      
      
      for (idx_var in names(param$list_tanner)){
        print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_var]][["label"]],sep=""))
        list_covar<-param$list_covar_tanner
        list_covar[["tanner"]]<-param$list_tanner[[idx_var]]
        
        df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                         prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
        df_clin$wave<-as.character(df_clin$wave)
        df_join<-inner_join(df_mean_sd,df_clin,by=c("wave","ID_pnTTC"))
        for (var_clin in c("age","tanner")){
          for (var_fc in c("mean","sd")){
            df_plot<-df_join[,c("wave","ID_pnTTC",var_fc,var_clin)]
            colnames(df_plot)[c(3,4)]<-c("var_fc","var_clin")
            df_plot$clin<-as.numeric.factor(df_plot$clin)
            #plt<-(ggplot(data=df_plot,aes(x="clin",y="mean",group="ID_pnTTC"))
            plt<-(ggplot(data=df_plot,aes(x=var_clin,y=var_fc,group=ID_pnTTC))
                  + geom_point()
                  + geom_line()
                  #+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  #               position=position_dodge(0.4))
                  + ggtitle(paste(dir_src,atlas,idx_var,var_clin,var_fc,sep=' '))
                  + xlab(var_clin)
                  + ylab(var_fc)
                  + theme_light()
                  + theme(plot.title = element_text(hjust = 0.5))
            )
            ggsave(paste("preproc-",dir_src,"_atl-",atlas,"_var-",idx_var,"_clin-",var_clin,"_",var_fc,"_fc.png",sep=""),
                   plot=plt,path=file.path(paths_$output,"output","plot"),height=5,width=7,dpi=300)
          }
        }
      }
    }
  }
}


#**************************************************
# GLMM/GAMM of FC =================================
#**************************************************
gamm_fc_core<-function(paths,data_fc,atlas,param,
                       list_covar,list_mod,list_term,idx_var,
                       calc_parallel,test_mod
){
  # Prepare clinical data and demean
  df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                   prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
  # Select subjects with longitudinal data
  if (param$force_long){
    list_id_subj<-df_clin[df_clin$wave==param$list_wave[1],'ID_pnTTC']
    for (wave in param$list_wave[-1]){
      list_id_subj<-sort(intersect(list_id_subj,df_clin[df_clin$wave==wave,'ID_pnTTC']))
    }
    df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
  }
  # Select subjects with non-decreasing data
  if (!is.null(param$omit_decreasing)){
    list_id_subj<-sort(unique(df_clin$ID_pnTTC))
    list_id_subj_omit<-NULL
    for (var in param$omit_decreasing){
      if (var %in% colnames(df_clin)){
        for (id_subj in list_id_subj){
          value<-as.numeric.factor(df_clin[df_clin$ID_pnTTC==id_subj & df_clin$wave==param$list_wave[1],var])
          for (wave in param$list_wave[-1]){
            value_wave<-as.numeric.factor(df_clin[df_clin$ID_pnTTC==id_subj & df_clin$wave==wave,var])
            if (value_wave<value){
              list_id_subj_omit<-c(list_id_subj_omit,id_subj)
            }
            value<-value_wave
          }
        }
      }
    }
    list_id_subj_omit<-sort(unique(list_id_subj_omit))
    list_id_subj<-list_id_subj[list_id_subj %nin% list_id_subj_omit]
    df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
  }
  # Group Tanner stage data
  if (!is.null(param$group_tanner)){
    df_clin$tanner<-as.numeric.factor(df_clin$tanner)
    for (name_group in names(param$group_tanner)){
      list_tanner<-param$group_tanner[[name_group]]
      df_clin[df_clin$tanner %in% list_tanner,"tanner"]<-name_group
    }
    df_clin$tanner<-factor(df_clin$tanner,levels=names(param$group_tanner))
  }
  
  df_clin<-func_std_clin(df_clin,separate_sex=T)$df_clin
  fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_clin.csv",sep="")),row.names=F)
  
  # Prepare FC data
  df_fc<-data_fc$df_fc; df_fc_grp<-data_fc$df_fc_grp

  label_wave<-"long"
  # Calculate model
  data_gamm<-func_calc_gamm(paths,df_clin,df_fc,df_fc_grp,data_fc,calc_parallel,test_mod,
                            atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
  df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova_grp<-data_gamm$df_anova_grp
  
  if (param$tfnbs){
    if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
    #if (calc_parallel){clust<-makeCluster(floor(detectCores()*1/2))}else{clust<-makeCluster(1)}
    list_sex<-param$list_sex
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                                  "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC",
                                  "param","func_bfs","%nin%","left_join","ggsave"),
                  envir=environment())
    
    # Calculate threshold-free network based statistics
    data_tfnbs<-func_iterate_tfnbs2(paths,clust,df_gamm,df_anova,df_deltah_in=NULL,var_exp_perm_in=NULL,data_fc,plot_result=T,return_nbs=T,
                                    atlas,param,list_mod,list_term,idx_var,label_wave)
    # Permutation test
    func_tfnbs_permutation(paths,clust,data_fc,df_clin,data_tfnbs_in=data_tfnbs,calc_parallel,plot_result=T,
                           atlas,param,list_mod,list_term,idx_var,label_wave)
    stopCluster(clust)
  }else{
    # Threshold and plot graph edges
    data_plot<-func_threshold_gamm(paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                                   atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
    df_plot<-data_plot$df_plot; df_plot_grp<-data_plot$df_plot_grp
    
    # Detect sub-network by breadth-first approach
    data_bfs<-func_detect_subnet(paths,df_plot,df_gamm,data_fc,plot_result=F,
                                  atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
    # Permutation test
    data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,data_bfs,data_fc,calc_parallel,plot_result=T,
                                   atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
  }
}


gamm_fc<-function(paths_=paths,list_atlas_=list_atlas,param=param_gamm_fc){
  
  print("Starting gamm_fc().")
  nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    print(paste("Preparing FC data: ",atlas,sep=""))
    #data_fc<-prep_data_fc(paths_,atlas,param$key_group,abs_nfc=param$abs_nfc)
    data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=c("1","2"),include_grp=T,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
    data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))
    fwrite(data_fc$df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_src_fc.csv",sep="")),row.names=F)
    fwrite(data_fc$df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_src_fc_grp.csv",sep="")),row.names=F)
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_tanner in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gamm_fc_core(paths_,data_fc,atlas,param,list_covar,
                   list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
                   calc_parallel=T,test_mod=F)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(param$list_hormone)){
      print(paste("Atlas: ",atlas,", Hormone type: ",param$list_hormone[[idx_hormone]][["label"]],sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
      gamm_fc_core(paths_,data_fc,atlas,param,list_covar,
                   list_mod=param$list_mod_hormone,list_term=param$list_term_hormone,idx_var=idx_hormone,
                   calc_parallel=T,test_mod=F)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  list_var<-c(param$list_tanner,param$list_hormone)
  func_combine_result(paths_,list_atlas_,list_var,"long",list(list("measure"="")),c("gamm","plot","gamm_anova","gamm_aic","gamm_grp","plot_grp","gamm_anova_grp","gamm_aic_grp","bfs_edge","bfs_node","bfs_size","bfs_pred","perm_max","perm_thr","perm_fwep","perm_sign","tfnbs"))
  print("Finished gamm_fc().")
}


#**************************************************
# GLMM/GAMM of FC mix both sex ====================
#**************************************************
gamm_fc_mix_core<-function(paths,data_fc,atlas,param,
                           list_covar,list_mod,list_term,idx_var,
                           calc_parallel,test_mod){
  # Prepare clinical data and demean
  df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=F,
                                   prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
  df_clin<-df_clin[!(is.na(df_clin$tanner_m) & is.na(df_clin$tanner_f)),]
  # Fill Tanner=NA values with 1
  if (param$fill_na_tanner){
    for (col in c("tanner_m","tanner_f")){
      if (col %in% colnames(df_clin)){
        df_clin[is.na(df_clin[col]),col]<-1
      }
    }
  }
  # Select subjects with longitudinal data
  if (param$force_long){
    list_id_subj<-df_clin[df_clin$wave==param$list_wave[1],'ID_pnTTC']
    for (wave in param$list_wave[-1]){
      list_id_subj<-sort(intersect(list_id_subj,df_clin[df_clin$wave==wave,'ID_pnTTC']))
    }
    df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
  }
  # Select subjects with non-decreasing data
  if (!is.null(param$omit_decreasing)){
    list_id_subj<-sort(unique(df_clin$ID_pnTTC))
    list_id_subj_omit<-NULL
    for (var in param$omit_decreasing){
      if (var %in% colnames(df_clin)){
        for (id_subj in list_id_subj){
          value<-as.numeric.factor(df_clin[df_clin$ID_pnTTC==id_subj & df_clin$wave==param$list_wave[1],var])
          for (wave in param$list_wave[-1]){
            value_wave<-as.numeric.factor(df_clin[df_clin$ID_pnTTC==id_subj & df_clin$wave==wave,var])
            if (value_wave<value){
              list_id_subj_omit<-c(list_id_subj_omit,id_subj)
            }
            value<-value_wave
          }
        }
      }
    }
    list_id_subj_omit<-sort(unique(list_id_subj_omit))
    list_id_subj<-list_id_subj[list_id_subj %nin% list_id_subj_omit]
    df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
  }
  #df_clin<-func_demean_clin(df_clin,separate_sex=F)$df_clin # separate_sex=F in analysis with mixed sex
  df_clin<-func_std_clin(df_clin,separate_sex=F)$df_clin # separate_sex=F in analysis with mixed sex
  fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_clin.csv",sep="")),row.names=F)
  
  # Prepare FC data
  df_fc<-data_fc$df_fc; df_fc_grp<-data_fc$df_fc_grp
  fwrite(df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_fc.csv",sep="")),row.names=F)
  fwrite(df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_fc_grp.csv",sep="")),row.names=F)
  
  label_wave<-"long"
  # Calculate model
  data_gamm<-func_calc_gamm(paths,df_clin,df_fc,df_fc_grp,data_fc,calc_parallel,test_mod,
                            atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
  df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova_grp<-data_gamm$df_anova_grp
  # Threshold and plot graph edges
  data_plot<-func_threshold_gamm(paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                                 atlas,param,"1_2",list_covar,list_mod,list_term,idx_var,label_wave)
  df_plot<-data_plot$df_plot; df_plot_grp<-data_plot$df_plot_grp
  # Detect sub-network by breadth-first approach
  data_bfs<-func_detect_subnet(paths,df_plot,df_gamm,data_fc,plot_result=F,
                                atlas,param,"1_2",list_covar,list_mod,list_term,idx_var,label_wave)
  # Permutation test
  df_clin$sex<-"1_2"
  data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,data_bfs,data_fc,calc_parallel,plot_result=T,
                                 atlas,param,"1_2",list_covar,list_mod,list_term,idx_var,label_wave)
}


gamm_fc_mix<-function(paths_=paths,list_atlas_=list_atlas,param=param_gamm_fc_mix){
  print("Starting gamm_fc_mix().")
  nullobj<-func_createdirs(paths_,str_proc="gamm_fc_mix()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    print(paste("Preparing FC data: ",atlas,sep=""))
    data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=c("1","2"),include_grp=T,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
    data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_tanner in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",idx_tanner,sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner_m"]]<-param$list_tanner[[idx_tanner]][["male"]]
      list_covar[["tanner_f"]]<-param$list_tanner[[idx_tanner]][["female"]]
      gamm_fc_mix_core(paths_,data_fc,atlas,param,list_covar,
                       list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
                       calc_parallel=T,test_mod=F)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(param$list_hormone)){
      print(paste("Atlas: ",atlas,", Hormone type: ",idx_hormone,sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone_m"]]<-param$list_hormone[[idx_hormone]][["male"]]
      list_covar[["hormone_f"]]<-param$list_hormone[[idx_hormone]][["female"]]
      gamm_fc_mix_core(paths_,data_fc,atlas,param,list_covar,
                       list_mod=param$list_mod_hormone,list_term=param$list_term_hormone,idx_var=idx_hormone,
                       calc_parallel=T,test_mod=F)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  list_var<-c(param$list_tanner,param$list_hormone)
  func_combine_result(paths_,list_atlas_,list_var,"long",list(list("measure"="")),c("gamm","plot","gamm_anova","gamm_aic","gamm_grp","plot_grp","gamm_anova_grp","gamm_aic_grp","bfs_edge","bfs_node","bfs_size","bfs_pred","perm_max","perm_thr","perm_fwep"))
  
  print("Finished gamm_fc_mix().")
}



#**************************************************
# Component Analyses of FC ========================
#**************************************************

group_factor<-function(df_ca_mri,dim,df_roi,df_grp,list_sex){
  df_ca_mri<-dplyr::inner_join(df_ca_mri,df_roi[,c("id","group")],by=c("from"="id"))
  df_ca_mri<-dplyr::rename(df_ca_mri,"from_group"="group")
  df_ca_mri<-dplyr::inner_join(df_ca_mri,df_roi[,c("id","group")],by=c("to"="id"))
  df_ca_mri<-dplyr::rename(df_ca_mri,"to_group"="group")
  
  df_ca_mri_grp<-NULL
  for (label_sex in names(list_sex)){
    for (idx_grp1 in seq(nrow(df_grp))){
      for (idx_grp2 in seq(idx_grp1,nrow(df_grp))){
        df_ca_mri_subset<-rbind(df_ca_mri[df_ca_mri$sex==label_sex
                                          & df_ca_mri$from_group==df_grp[idx_grp1,"id"]
                                          & df_ca_mri$to_group==df_grp[idx_grp2,"id"],],
                                df_ca_mri[df_ca_mri$sex==label_sex
                                          & df_ca_mri$from_group==df_grp[idx_grp2,"id"]
                                          & df_ca_mri$to_group==df_grp[idx_grp1,"id"],])
        df_ca_mri_subset<-df_ca_mri_subset[,sprintf("comp_%03d",seq(dim))]
        df_ca_mri_grp<-rbind(df_ca_mri_grp,
                             cbind(sex=label_sex,abs=F,from=df_grp[idx_grp1,"id"],to=df_grp[idx_grp2,"id"],
                                   label_from=df_grp[idx_grp1,"label"],label_to=df_grp[idx_grp2,"label"],
                                   t(colMeans(df_ca_mri_subset))))
        df_ca_mri_subset_abs<-abs(df_ca_mri_subset)
        df_ca_mri_grp<-rbind(df_ca_mri_grp,
                             cbind(sex=label_sex,abs=T,from=df_grp[idx_grp1,"id"],to=df_grp[idx_grp2,"id"],
                                   label_from=df_grp[idx_grp1,"label"],label_to=df_grp[idx_grp2,"label"],
                                   t(colMeans(df_ca_mri_subset_abs))))
        
      }
    }
  }
  df_ca_mri_grp<-as.data.frame(cbind(dim=dim,df_ca_mri_grp))
  return(df_ca_mri_grp)
}

comp_clin_cor<-function(df_comp_subj,df_clin,n_covar,list_sex,atlas,method,
                        wave_mri,wave_clin,idx_var,label_var){
  list_dim<-sort(unique(df_comp_subj$dim))
  # Divide df_clin into covariates for Pearson correlation and Spearman rank correlation
  # criteria is if the covariate has less than 10 unique values
  df_clin_pearson<-df_clin_spearman<-df_clin[,"ID_pnTTC",drop=F]
  for (covar_clin in colnames(df_clin)[2:ncol(df_clin)]){
    n_uniquevalues<-length(unique(df_clin[,covar_clin]))
    if(n_uniquevalues<10){
      df_clin_spearman<-cbind(df_clin_spearman,df_clin[,covar_clin,drop=F])
    }else{
      df_clin_pearson<-cbind(df_clin_pearson,df_clin[,covar_clin,drop=F])
    }
  }
  list_df_clin<-list("pearson"=df_clin_pearson,"spearman"=df_clin_spearman)
  
  df_cor_flat_rbind<-NULL
  for (dim_ca in list_dim){
    for (label_sex in names(list_sex)){
      df_comp_subj_subset<-df_comp_subj[df_comp_subj$dim==dim_ca & df_comp_subj$sex==label_sex,
                                        c("ID_pnTTC",sprintf("comp_%03d",1:dim_ca))]
      df_plot<-NULL
      for (type_corr in names(list_df_clin)){
        df_clin_subset<-list_df_clin[[type_corr]]
        df_join<-inner_join(df_clin_subset,df_comp_subj_subset,by="ID_pnTTC")
        df_join$ID_pnTTC<-NULL
        data_cor<-func_cor(df_join,type=type_corr)
        
        # Pick up correlation between covariates and factors
        df_cor_flat<-data_cor$cor_flat
        df_cor_flat<-df_cor_flat[df_cor_flat$from %in% colnames(df_clin_subset) & df_cor_flat$to %in% colnames(df_comp_subj_subset),]
        df_cor_flat<-df_cor_flat[,c("from","to","r","p")]
        colnames(df_cor_flat)<-c("covar","comp","r","p")
        df_cor_flat$comp<-sapply(sapply(df_cor_flat$comp,substr,start=6,stop=8),as.integer)
        rownames(df_cor_flat)<-NULL
        
        # rbind data
        df_plot<-rbind(df_plot,df_cor_flat)
        df_cor_flat_rbind<-rbind(df_cor_flat_rbind,cbind(sex=label_sex,dim=dim_ca,df_cor_flat))
      }
      
      ## Visualize result
      #df_plot$r<-as.numeric(df_plot$r)
      #df_plot<-df_plot[df_plot$covar!="sex",]
      #df_plot$covar<-str_to_title(df_plot$covar)
      #df_plot$significance<-""
      #df_plot[df_plot$p<0.05,"significance"]<-"*"
      #df_plot[df_plot$p<0.001,"significance"]<-"**"
      #
      #plot<-(ggplot(data=df_plot,aes(x=comp,y=r,fill=covar))
      #       + geom_bar(stat="identity",color="white",width=0.7,position=position_dodge())
      #       + scale_fill_brewer(palette="Accent")
      #       + scale_x_continuous(breaks=seq(dim_ca),limits=c(0.5,max(list_dim)+0.5),expand=c(0.003,0.003))
      #       + scale_y_continuous(breaks=seq(-0.5,0.5,0.1),limits=c(-0.5,0.5))
      #       + geom_hline(yintercept = 0, linetype = 2)
      #       + geom_text(
      #         aes(label = significance, group = covar), 
      #         position = position_dodge(0.7),
      #         vjust = 0, size = 3.5
      #       )
      #       + ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Variable: ",label_var,
      #                       ", Clinical: ",wave_clin, ",MRI: ",wave_mri,
      #                       ", Dimension: ",as.character(dim_ca),
      #                       ", Sex: ",label_sex,sep=""))
      #       + xlab("Factor") + ylab("r or rho") + theme_classic()
      #       + theme(plot.title = element_text(hjust = 0.5),
      #               legend.position="top",legend.justification="center",legend.direction="horizontal",legend.title=element_blank())
      #)
      #name_file<-paste("atl-",atlas,"_method-",method,"_var-",idx_var,"_ses-c",wave_clin,"m",wave_mri,
      #                 "_sex-",label_sex,"_dim-",sprintf("%03d",dim_ca),
      #                 "_fc_ca_cor.png",sep="")
      #ggsave(name_file,plot,path=file.path(paths$output,"output","plot"),height=5,width=10,dpi=200)
    }
  }
  return(list("df_cor_flat"=df_cor_flat_rbind))
}

ca_fc_cs<-function(paths_=paths,list_atlas_=list_atlas,param=param_ca_fc_cs,skip_ca_plot=F){
  print("Starting ca_fc_cs()")
  nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs()",copy_log=T,list_param=param)
  # Increase memory limit for later ICA calculation
  memory.limit(1000000)
  for (atlas in list_atlas_){
    # Load and examine FC data
    data_fc<-prep_data_fc2(paths,atlas,param$key_group,list_wave=param$list_wave_mri,include_grp=F,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    df_fc<-data_fc$df_fc; df_edge<-data_fc$df_edge
    
    # Calculate PCA/ICA factors (if not yet)
    for (wave_mri in param$list_wave_mri){
      path_pca_subj<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep=""))
      #path_ica_subj<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep=""))
      
      #if (!(file.exists(path_pca_subj) & file.exists(path_ica_subj))){
      if (!(file.exists(path_pca_subj))){
        #df_pca_mri<-df_pca_subj<-df_pca_vaf<-df_ica_mri<-df_ica_subj<-df_ica_vaf<-NULL
        df_pca_mri<-df_pca_subj<-df_pca_vaf<-NULL
        for (label_sex in names(param$list_sex)){
          print(paste("Calculating Atlas: ",atlas,", MRI wave: ",wave_mri,", Sex: ",label_sex,sep=""))
          
          # Prepare subject subsetting condition (MRI QC criteria and sex) according to specified MRI wave
          # Existence of clinical variables are not considered here
          subset_subj_temp<-list(c(param$subset_subj[[wave_mri]],list(list("key"="Sex","condition"=param$list_sex[[label_sex]]))))
          names(subset_subj_temp)<-wave_mri
          df_clin<-func_clinical_data_long(paths_,wave_mri,subset_subj_temp,list_covar=NULL,rem_na_clin=F,
                                             prefix=paste("wav-m",wave_mri,"_sex-",label_sex,sep=""),print_terminal=F)$df_clin
          df_clin<-dplyr::rename(df_clin,"ses"="wave")
          
          df_fc_ses<-df_fc[df_fc$ses==wave_mri,]
          df_fc_calc<-df_clin_exist<-NULL
          for (id_subj in df_clin$ID_pnTTC){
            df_fc_subj<-df_fc_ses[df_fc_ses$ID_pnTTC==id_subj,"z_r"]
            df_fc_calc<-rbind(df_fc_calc,df_fc_subj)
            df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ID_pnTTC==id_subj,])
          }
          rownames(df_fc_calc)<-as.character(seq(nrow(df_fc_calc))); colnames(df_fc_calc)<-as.character(seq(ncol(df_fc_calc)))
          colnames(df_clin_exist)<-colnames(df_clin)
          df_clin_exist$ses<-NULL
          gc()
          
          # Calculate PCA of FC
          dim_ca<-max(param$list_dim_ca)
          data_pca<-func_pca(df_src=df_fc_calc,df_var=data_fc$df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca,calc_corr=F)
          
          # Save results
          df_pca_mri<-rbind(df_pca_mri,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_mri))
          df_pca_subj<-rbind(df_pca_subj,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_subj))
          df_pca_vaf<-rbind(df_pca_vaf,cbind(sex=label_sex,dim=dim_ca,data_pca$df_vaf))
          data_pca<-NULL
          gc()
          
          # Calculate ICA of FC
          #for (dim_ca in list_dim_ca_){
          #  data_ica<-func_ica(df_src=df_fc_calc,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca,calc_corr=F)
          #  
          #  # Save results
          #  df_ica_mri<-rbind.fill(df_ica_mri,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_mri))
          #  df_ica_subj<-rbind.fill(df_ica_subj,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_subj))
          #  df_ica_vaf<-rbind.fill(df_ica_vaf,cbind(sex=label_sex,dim=dim_ca,data_ica$df_vaf))
          #  data_ica<-NULL
          #  gc()
          #} # end for over ICA dimensions
        } # end for over sex
        fwrite(df_pca_mri,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep="")),row.names=F)
        fwrite(df_pca_subj,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep="")),row.names=F)
        fwrite(df_pca_vaf,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_vaf.csv",sep="")),row.names=F)
        #fwrite(df_ica_mri,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep="")),row.names=F)
        #fwrite(df_ica_subj,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep="")),row.names=F)
        #fwrite(df_ica_vaf,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_vaf.csv",sep="")),row.names=F)
      } # end if not PCA/ICA factors are already calculated
    } # end for waves
    # end calculating PCA/ICA factors
    
    # Group-wise average of factor-MRI matrix
    print(paste("Calculating group-wise contribution to factors:",atlas,sep=" "))

    for (wave_mri in param$list_wave_mri){
      path_pca_mri_grp<-file.path(paths_$output,"output","temp", paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var_grp.csv",sep=""))
      #path_ica_mri_grp<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var_grp.csv",sep=""))
      #if (!(file.exists(path_pca_mri_grp) & file.exists(path_ica_mri_grp))){
      if (!(file.exists(path_pca_mri_grp))){
        df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
        #df_ica_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep=""))))
        #df_pca_mri_grp<-df_ica_mri_grp<-data.frame()
        df_pca_mri_grp<-NULL
        # PCA
        dim<-max(param$list_dim_ca)
        df_pca_mri_grp<-group_factor(df_pca_mri,dim,data_fc$df_roi,data_fc$df_grp,param$list_sex)
        fwrite(df_pca_mri_grp,path_pca_mri_grp,row.names=F)
        # ICA
        #for (dim in list_dim_ca_){
        #  df_ica_mri_subset<-df_ica_mri[df_ica_mri$dim==dim,]
        #  df_ica_mri_grp<-rbind.fill(df_ica_mri_grp,
        #                             group_factor(df_ica_mri_subset,dim,dict_roi,data_fc$df_grp$id,list_sex_))
        #}
        #fwrite(df_ica_mri_grp,path_ica_mri_grp,row.names=F)
      }
    }
    # end calculating group-wise average
    
    # Generate visual output of MRI factors
    if (skip_ca_plot){
      print(paste("Heatmap plot of PCA/ICA factors skipped:",atlas,sep=" "))
    }else{
      print(paste("Generating heatmap plot of PCA/ICA factors:",atlas,sep=" "))
      for (wave_mri in param$list_wave_mri){
        df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
        df_pca_mri_grp<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var_grp.csv",sep=""))))
        #df_ica_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep=""))))
        #df_ica_mri_grp<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var_grp.csv",sep=""))))
        
        # Visual output of PCA/ICA factors
        for (label_sex in names(param$list_sex)){
          # PCA
          dim_ca<-max(param$list_dim_ca)
          df_pca_mri_subset<-df_pca_mri[df_pca_mri$sex==label_sex & df_pca_mri$dim==dim_ca,]
          df_pca_mri_grp_subset<-df_pca_mri_grp[df_pca_mri_grp$sex==label_sex & df_pca_mri_grp$dim==dim_ca,]
          df_pca_mri_subset$sex<-df_pca_mri_subset$dim<-df_pca_mri_grp_subset$sex<-df_pca_mri_grp_subset$dim<-NULL
          # Visualize factor-FC matrix in heatmap plot
          plot_ca_fc_heatmap(paths_=paths_,data_fc,df_pca_mri_subset,df_pca_mri_grp_subset,atlas=atlas,dim_ca=dim_ca,
                             method="pca",label_sex=label_sex,ses=wave_mri)
          
          # ICA
          #for (dim_ca in list_dim_ca_){
          #  df_ica_mri_subset<-df_ica_mri[df_ica_mri$sex==label_sex & df_ica_mri$dim==dim_ca,]
          #  df_ica_mri_grp_subset<-df_ica_mri_grp[df_ica_mri_grp$sex==label_sex & df_ica_mri_grp$dim==dim_ca,]
          #  df_ica_mri_subset$sex<-df_ica_mri_subset$dim<-df_ica_mri_grp_subset$sex<-df_ica_mri_grp_subset$dim<-NULL
          #  # Visualize factor-FC matrix in heatmap plot
          #  plot_ca_fc_heatmap(paths_=paths_,df_ica_mri_subset,df_ica_mri_grp_subset,atlas=atlas,dim_ca=dim_ca,
          #                     method="ica",label_sex=label_sex,ses=wave_mri)
          #}
        }
      } # End for MRI wave
    }
    # End of generating factor heatmap plot
    
    # Calculate node strength in each factor
    print(paste("Calculating node strength of PCA/ICA factors:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      path_file_check<-file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_str.csv",sep=""))
      if (!file.exists(path_file_check)){
        df_pca_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_var.csv",sep=""))))
        #df_ica_mri<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_var.csv",sep=""))))
        df_strength<-NULL
        #list_ca_mri<-list("pca"=df_pca_mri,"ica"=df_ica_mri)
        list_ca_mri<-list("pca"=df_pca_mri)
        for (method in names(list_ca_mri)){
          df_ca_mri<-list_ca_mri[[method]]
          list_sex<-unique(df_ca_mri$sex)
          list_dim<-unique(df_ca_mri$dim)
          for (sex in list_sex){
            for (dim in list_dim){
              df_ca_mri_subset<-df_ca_mri[df_ca_mri$sex==sex & df_ca_mri$dim==dim,]
              for (idx_row in seq(nrow(data_fc$df_roi))){
                node<-data_fc$df_roi[idx_row,"id"]
                label_node<-data_fc$df_roi[idx_row,"label"]
                strength<-colSums(abs(df_ca_mri_subset[df_ca_mri_subset$from==node | df_ca_mri_subset$to==node,
                                                       sprintf("comp_%03d",seq(dim))]))
                df_strength_add<-data.frame(method=method,sex=sex,dim=dim,comp=seq(dim),node=node,label_node=label_node,strength=strength)
                df_strength<-rbind(df_strength,df_strength_add)
              }
            }
          }
        }
        fwrite(df_strength,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_str.csv",sep="")),row.names=F)
      }
    }
    # End calculating node strength in each factor
    
    # Calculate factor attibution-clinical relation
    print(paste("Calculating factor-clinical correlation:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      if (!file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_cor.csv",sep="")))){
        df_pca_subj<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep=""))))
        #df_ica_subj<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep=""))))
        df_cor<-NULL
        for (wave_clin in param$list_wave_clin){
          # Pickup subsetting condition of MRI wave, and rename it to clinical wave
          subset_subj_temp<-list(param$subset_subj[[wave_mri]])
          names(subset_subj_temp)<-wave_clin
          #1 Tanner stage
          for (idx_tanner in names(param$list_tanner)){
            #print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
            list_covar<-param$list_covar_tanner
            list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
            n_covar<-length(list_covar)
            prefix<-paste("wav-c",wave_clin,"m",wave_mri,"_var-",idx_tanner,sep="")
            df_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                             rem_na_clin=T,prefix=prefix,print_terminal=F)$df_clin
            df_clin$wave<-NULL
            
            # Calculate correlation between component attribution and clinical covariate
            data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                                    n_covar=n_covar,list_sex=param$list_sex,atlas=atlas,method="pca",
                                    wave_mri=wave_mri,wave_clin=wave_clin,
                                    idx_var=idx_tanner,label_var=param$list_tanner[[idx_tanner]]$label)
            df_cor<-rbind(df_cor,cbind(wave_clin=paste("c",wave_clin,sep=""),
                                       variable=idx_tanner,method="pca",data_cor$df_cor_flat))
            #data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
            #                        n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="ica",
            #                        wave_mri=wave_mri,wave_clin=wave_clin,
            #                        idx_var=idx_tanner,label_var=list_tanner_[[idx_tanner]]$label)
            #df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
            #                           variable=idx_tanner,method="ica",data_cor$df_cor_flat))
          } # End for Tanner stages
          
          #2 Hormones
          for (idx_hormone in names(param$list_hormone)){
            list_covar<-param$list_covar_hormone
            list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
            n_covar<-length(list_covar)
            prefix<-paste("wav-c",wave_clin,"m",wave_mri,"_var-",idx_hormone,sep="")
            df_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                             rem_na_clin=T,prefix=prefix,print_terminal=F)$df_clin
            df_clin$wave<-NULL
            
            # Calculate correlation between component attribution and clinical covariate
            data_cor<-comp_clin_cor(df_comp_subj=df_pca_subj,df_clin=df_clin,
                                    n_covar=n_covar,list_sex=param$list_sex,atlas=atlas,method="pca",
                                    wave_mri=wave_mri,wave_clin=wave_clin,
                                    idx_var=idx_hormone,label_var=param$list_hormone[[idx_hormone]]$label)
            df_cor<-rbind(df_cor,cbind(wave_clin=paste("c",wave_clin,sep=""),
                                       variable=idx_hormone,method="pca",data_cor$df_cor_flat))
            #data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
            #                        n_covar=n_covar,list_sex=list_sex_,atlas=atlas,method="ica",
            #                        wave_mri=wave_mri,wave_clin=wave_clin,
            #                        idx_var=idx_hormone,label_var=list_hormone_[[idx_hormone]]$label)
            #df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=paste("c",as.character(wave_clin),"m",as.character(wave_mri),sep=""),
            #                           variable=idx_hormone,method="ica",data_cor$df_cor_flat))
          } # End for hormones
        } # End for clinical wave
        fwrite(df_cor,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_cor.csv",sep="")),row.names=F)
      } # End if calculated result does not exist
    } # End for MRI wave
    # End of calculating factor-clinical correlation
    
    # Integrated visual output of factor-clinical correlation
    print(paste("Generating integrated graph of factor-clinical correlation:",atlas,sep=" "))
    list_var<-c(param$list_tanner,param$list_hormone)
    for (wave_mri in param$list_wave_mri){
      df_cor<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_cor.csv",sep=""))))
      for (idx_var in names(list_var)){
        for (label_sex in names(param$list_sex)){
          list_method<-"pca"
          for (method in list_method){
            list_dim<-sort(unique(df_cor[df_cor$variable==idx_var & df_cor$sex==label_sex & df_cor$method==method,"dim"]))
            for (dim in list_dim){
              df_plot<-NULL
              for (wave_clin in param$list_wave_clin){
                label_wave_clin<-paste("c",wave_clin,sep="")
                df_cor_subset<-df_cor[df_cor$wave_clin==label_wave_clin & df_cor$variable==idx_var
                                      & df_cor$sex==label_sex & df_cor$method==method & df_cor$dim==dim,]
                df_plot<-rbind(df_plot,df_cor_subset)
              }
              fig_ca_clin<-plot_ca_clin(df_plot,list_var,idx_var)
              fig_ca_clin<-(fig_ca_clin
                            + ggtitle(paste("Atlas: ",atlas,", Method: ",method,", Variable: ",idx_var,", MRI wave: ",wave_mri,", Dimension: ",as.character(dim),", Sex: ",label_sex,sep="")))
              ggsave(paste("atl-",atlas,"_method-",method,"_var-",idx_var,"_wav-m",wave_mri,"_sex-",label_sex,"_dim-",sprintf("%03d",dim),"_fc_ca_cor_int.png",sep=""),
                     fig_ca_clin,path=file.path(paths_$output,"output","plot"),height=5,width=10,dpi=200)
            }
          }
        }
      }
    } # End of integrated visual output of factor-clinical correlation
    
    # Calculate factor attribution-clinical GLM/GAM
    print(paste("Calculating factor-clinical GLM/GAM:",atlas,sep=" "))
    for (wave_mri in param$list_wave_mri){
      if (!file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_aic.csv",sep="")))){
        df_pca_subj<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_pca_subj.csv",sep=""))))
        #df_ica_subj<-as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_ses-m",wave_mri,"_fc_ica_subj.csv",sep=""))))
        df_gamm<-df_anova<-df_aic<-NULL
        for (wave_clin in param$list_wave_clin){
          # Pickup subsetting condition of MRI wave, and rename it to clinical wave
          subset_subj_temp<-list(param$subset_subj[[wave_mri]])
          names(subset_subj_temp)<-wave_clin
          #1 Tanner stage
          for (idx_tanner in names(param$list_tanner)){
            list_covar<-param$list_covar_tanner
            list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
            n_covar<-length(list_covar)
            prefix<-paste("wav-c",wave_clin,"m",wave_mri,"_var-",idx_tanner,sep="")
            df_clin<-func_clinical_data_long(paths_,wave_clin,subset_subj_temp,list_covar,
                                             rem_na_clin=T,prefix=prefix,print_terminal=F)$df_clin
            df_clin<-func_demean_clin(df_clin,separate_sex = F)$df_clin
            df_clin$wave<-NULL
            
            for (label_sex in names(param$list_sex)){
              list_sex<-c(1,2)
              eval(parse(text=paste('list_sex<-list_sex[list_sex',param$list_sex[[label_sex]],']',sep='')))
              list_sex<-list(list_sex)
              for (id_comp in seq(max(param$list_dim_ca))){
                df_pca_subset<-df_pca_subj[df_pca_subj$sex==label_sex,c("ID_pnTTC",sprintf("comp_%03d",id_comp))]
                colnames(df_pca_subset)<-c("ID_pnTTC","value")
                df_join<-dplyr::inner_join(df_pca_subset,df_clin,by="ID_pnTTC")
                data_gamm<-gamm_core4(df_join,list_mod_in=param$list_mod_tanner,list_sex_in=list_sex,
                                      calc_parallel_in=F,test_mod_in=F)
                df_gamm_add<-data_gamm$df_gamm; df_anova_add<-data_gamm$df_anova; df_aic_add<-data_gamm$df_aic
                df_gamm_add$sex<-df_anova_add$sex<-df_aic_add$sex<-NULL
                df_head<-data.frame(wave_clin=paste("c",wave_clin,sep=""),variable=idx_tanner,method="pca",
                                    sex=label_sex,dim=max(param$list_dim_ca),comp=id_comp)
                df_gamm<-rbind(df_gamm,cbind(df_head,df_gamm_add)); df_anova<-rbind(df_anova,cbind(df_head,df_anova_add)); df_aic<-rbind(df_aic,cbind(df_head,df_aic_add))

              }
            }
          } # End for Tanner stages
          
        } # End for clinical wave
        fwrite(df_gamm,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_gamm.csv",sep="")),row.names=F)
        fwrite(df_anova,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_anova.csv",sep="")),row.names=F)
        fwrite(df_aic,file.path(paths_$output,"output","temp",paste("atl-",atlas,"_wav-m",wave_mri,"_fc_ca_aic.csv",sep="")),row.names=F)
      } # End if calculated result does not exist
    } # End for MRI wave
    
    
  } # End for atlas
  
  # Reload and bind all results
  print("Binding results.")
  func_combine_result(paths_,list_atlas_,list_var=NULL,list_wave=paste("m",param$list_wave_mri,sep=""),list(list("measure"="")),
                      list_filename=c("fc_pca_subj","fc_ica_subj","fc_pca_var","fc_ica_var","fc_pca_vaf","fc_ica_vaf",
                                      "fc_pca_var_grp","fc_ica_var_grp","fc_ca_cor","fc_ca_str",
                                      "fc_ca_gamm","fc_ca_anova","fc_ca_aic"))
  print("Finished ca_fc_cs()")
}




#**************************************************
# GLM/GAM of FC cross-section =====================
#**************************************************
gam_fc_cs_core<-function(paths,atlas,param,list_sex,
                         list_covar,list_mod,list_term,idx_var,
                         calc_parallel,test_mod
                         ){
  for (label_wave in names(param$list_wave)){
    print(paste("Atlas: ",atlas,", Measure: ",idx_var,", Wave: ",label_wave,sep=""))
    wave_clin<-param$list_wave[[label_wave]]$clin
    wave_mri<-param$list_wave[[label_wave]]$mri
    
    # Prepare clinical data and demean
    # QC subsetting condition must accord with MRI wave, but under the name of clinical wave
    subset_subj<-param$subset_subj[wave_mri]
    names(subset_subj)<-wave_clin
    data_clin<-func_clinical_data_long(paths,wave_clin,subset_subj,list_covar,rem_na_clin=T,prefix=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src",sep=""),print_terminal=F)
    #df_clin<-func_demean_clin(data_clin$df_clin,separate_sex=T)$df_clin
    df_clin<-func_std_clin(data_clin$df_clin,separate_sex=T)$df_clin
    fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_clin.csv",sep="")),row.names=F)
    df_clin$wave<-wave_mri # Need to meet MRI wave for later joining
    
    # Prepare FC data
    print(paste("Preparing FC data: ",atlas,sep=""))
    data_fc<-prep_data_fc2(paths,atlas,param$key_group,list_wave=wave_mri,include_grp=T,
                           abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
    df_fc<-data_fc$df_fc; df_fc_grp<-data_fc$df_fc_grp
    fwrite(df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_fc.csv",sep="")),row.names=F)
    fwrite(df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_fc_grp.csv",sep="")),row.names=F)
    
    # Calculate model
    data_gamm<-func_calc_gamm(paths,df_clin,df_fc,df_fc_grp,data_fc,calc_parallel,test_mod,
                              atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
    df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova_grp<-data_gamm$df_anova_grp
    # Threshold and plot graph edges
    data_plot<-func_threshold_gamm(paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                                   atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
    df_plot<-data_plot$df_plot; df_plot_grp<-data_plot$df_plot_grp
    # Detect sub-network by breadth-first approach
    data_bfs<-func_detect_subnet(paths,df_plot,df_gamm,data_fc,plot_result=F,
                                  atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
    # Permutation test
    data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,data_bfs,data_fc,calc_parallel,plot_result=T,
                                   atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
  }
}

gam_fc_cs<-function(paths_=paths,list_atlas_=list_atlas,param=param_gam_fc_cs){
  print("Starting gam_fc_cs().")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    # Loop over Tanner stage
    for (idx_tanner in names(param$list_tanner)){
      #print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gam_fc_cs_core(paths_,atlas,param,list_sex=list(1,2),list_covar,
                     list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
                     calc_parallel=T,test_mod=F)
    }
    # Loop over hormones
    for (idx_hormone in names(param$list_hormone)){
      #print(paste("Atlas: ",atlas,", Hormone type: ",param$list_hormone[[idx_hormone]][["label"]],sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
      gam_fc_cs_core(paths_,atlas,param,list_sex=list(1,2),list_covar,
                     list_mod=param$list_mod_hormone,list_term=param$list_term_hormone,idx_var=idx_hormone,
                     calc_parallel=T,test_mod=F)
    } 
  } # Finished looping over atlas
  
  print("Combining results.")
  list_var<-c(param$list_tanner,param$list_hormone)
  func_combine_result(paths_,list_atlas_,list_var,names(param$list_wave),list(list("measure"="")),c("gamm","plot","gamm_anova","gamm_aic","gamm_grp","plot_grp","gamm_anova_grp","gamm_aic_grp","bfs_edge","bfs_node","bfs_size","bfs_pred","perm_max","perm_thr","perm_fwep"))
  
  print("Finished gam_fc_cs().")
}


#**************************************************
# gam_fc(), gamm_fc() common functions ============
#**************************************************
func_calc_gamm<-function(paths,param,data_fc,df_clin,atlas,list_covar,list_mod,list_term,
                         idx_var,label_wave,calc_parallel,test_mod){
  
  file_check<-file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_aic_grp.csv",sep=""))
  if (file.exists(file_check)){
    print("Calculated GAMM/ANOVA already exists")
    df_gamm<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm.csv",sep=""))))
    df_anova<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_anova.csv",sep=""))))
    df_gamm_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_grp.csv",sep=""))))
    df_anova_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_anova_grp.csv",sep=""))))
  }else{
    print("Calculating GAMM/ANOVA")
    # Join FC and clinical data
    df_join<-join_fc_clin(data_fc$df_fc,df_clin)
    df_join_grp<-join_fc_clin(data_fc$df_fc_grp,df_clin)
    
    # Prepare parallelization cluster
    if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
    list_sex<-param$list_sex
    list_term_pred<-param$list_term_pred
    clusterExport(clust,varlist=c("list_mod","list_sex","list_term_pred","calc_parallel","test_mod","as.formula","as.numeric.factor",
                                  "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC","expand.grid","%in%","%nin%"),
                  envir=environment())
    
    # Calculate model
    data_gamm<-iterate_gamm5(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
    data_gamm_grp<-iterate_gamm5(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
    stopCluster(clust)
    
    # Add multiple comparison-corrected p values
    df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    
    # Save results
    fwrite(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm.csv",sep="")),row.names = F)
    fwrite(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_anova.csv",sep="")),row.names = F)
    fwrite(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_aic.csv",sep="")),row.names = F)
    fwrite(data_gamm$df_pred,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_pred.csv",sep="")),row.names = F)
    fwrite(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_grp.csv",sep="")),row.names = F)
    fwrite(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_anova_grp.csv",sep="")),row.names = F)
    fwrite(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_aic_grp.csv",sep="")),row.names = F)
    fwrite(data_gamm_grp$df_pred,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_gamm_pred_grp.csv",sep="")),row.names = F)
  } # End if file exists
  
  return(list("df_gamm"=df_gamm,"df_anova"=df_anova,"df_pred"=data_gamm$df_pred,"df_gamm_grp"=df_gamm_grp,"df_anova_grp"=df_anova_grp,"df_pred_grp"=data_gamm_grp$df_pred))
}


func_threshold_gamm<-function(paths,param,data_gamm,data_fc,atlas,list_covar,list_mod,list_term,idx_var,label_wave){
                              #paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                              #atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave){
  if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_plot.csv",sep="")))){
    print("Thresholded GAMM/ANOVA already exists.")
    df_plot<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_plot.csv",sep=""))))
    if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot_grp.csv",sep="")))){
      df_plot_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_plot_grp.csv",sep=""))))
    }else{
      df_plot_grp<-data.frame()
    }
  }else{
    print("Thresholding GAMM/ANOVA results")
    df_gamm<-data_gamm$df_gamm; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova<-data_gamm$df_anova; df_anova_grp<-data_gamm$df_anova_grp
    list_plot<-list()
    df_plot<-df_plot_grp<-data.frame()
    for (idx_mod in names(list_mod)){
      for (idx_term in names(list_term)){
        var_exp<-list_term[[idx_term]][["var_exp"]]
        for (idx_sex in param$list_sex){
          # Subset GAMM result dataframe for plotting
          df_gamm_subset<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp & df_gamm$sex==idx_sex,]
          df_gamm_grp_subset<-df_gamm_grp[df_gamm_grp$model==idx_mod & df_gamm_grp$term==var_exp & df_gamm_grp$sex==idx_sex,]
          if (nrow(df_gamm_subset)==0){
            df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp & df_anova$sex==idx_sex,]
            df_anova_grp_subset<-df_anova_grp[df_anova_grp$model==idx_mod & df_anova_grp$term==var_exp & df_anova_grp$sex==idx_sex,]
            if (nrow(df_anova_subset)>0){
              # In case the term does not exist in df_gamm, plot using df_anova instead
              df_gamm_subset<-df_anova_subset
              df_gamm_grp_subset<-df_anova_grp_subset
            }
          }
          if (nrow(df_gamm_subset)>0){
            plot_gamm<-plot_gam_fc3(df_gamm_subset,df_gamm_grp_subset,data_fc)
            if (idx_sex==1){
              label_sex<-"m"
            }else if (idx_sex==2){
              label_sex<-"f"
            }else{
              label_sex<-"mf"
            }
            plot_gamm<-annotate_figure(plot_gamm,
                                       top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp,", sex: ",label_sex,", p value: all",sep=""),
                                                       color = "black", size = 14))
            list_plot<-c(list_plot,list(list("plot"=plot_gamm,"height"=13,"width"=10,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                             "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term,"_sex-",label_sex,"_pval-all_net.png",sep=""))))
            for (type_p in names(param$list_p)){
              list_thresh_p<-param$list_p[[type_p]]
              for (thresh_p in list_thresh_p){
                df_plot_subset<-df_gamm_subset[df_gamm_subset[[type_p]]<thresh_p,]
                df_plot_grp_subset<-df_gamm_grp_subset[df_gamm_grp_subset[[type_p]]<thresh_p,]
                plot_gamm<-plot_gam_fc3(df_plot_subset,df_plot_grp_subset,data_fc)
                plot_gamm<-annotate_figure(plot_gamm,
                                           top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp,", sex: ",label_sex,", p value: ",type_p,"<",thresh_p,sep=""),
                                                           color = "black", size = 14))
                list_plot<-c(list_plot,list(list("plot"=plot_gamm,"height"=13,"width"=10,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                 "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term,"_sex-",label_sex,"_pval-",type_p,"_",thresh_p,"_net.png",sep=""))))
                df_head<-data.frame(p_type=type_p,p_threshold=thresh_p)
                if (nrow(df_plot_subset)>0){
                  df_plot<-bind_rows(df_plot,cbind(df_head,df_plot_subset))
                }
                if (nrow(df_plot_grp_subset)>0){
                  df_plot_grp<-bind_rows(df_plot_grp,cbind(df_head,df_plot_grp_subset))
                }
              }
            }# end of loop over list_p
          }
        }
      }
    }
    
    if (length(list_plot)>0){
      clust<-makeCluster(floor(detectCores()*3/4))
      plot_parallel(clust,list_plot)
      stopCluster(clust)
    }
    
    # Save results
    if (nrow(df_plot)>0){
      fwrite(df_plot,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_plot.csv",sep="")),row.names = F)
    }
    if (nrow(df_plot_grp)>0){
      fwrite(df_plot_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_plot_grp.csv",sep="")),row.names = F)
    }
  }
  return(list("df_plot"=df_plot,"df_plot_grp"=df_plot_grp))
}

# parallelization
func_iterate_tfnbs2<-function(paths,clust,df_gamm,df_anova,df_deltah_in,var_exp_perm_in,data_fc,plot_result=F,return_nbs=T,
                              atlas,param,list_mod,list_term,idx_var,label_wave){
  list_src_tfnbs<-list()
  for (idx_mod in param$param_nbs$list_mod){
    for (set_term in param$param_nbs$list_term){
      idx_term_perm<-set_term$term_perm
      var_exp_perm<-list_term[[idx_term_perm]]$var_exp
      for (idx_term_detect in set_term$term_detect){
        var_exp_detect<-list_term[[idx_term_detect]]$var_exp
        if (!is.null(var_exp_detect)){
          for (idx_sex in param$list_sex){
            if (idx_sex==1){label_sex<-"m"}else if (idx_sex==2){label_sex<-"f"}else{label_sex<-"mf"}
            df_stat<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp_detect & df_gamm$sex==idx_sex,]
            if (nrow(df_stat)==0){
              df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp_detect & df_anova$sex==idx_sex,]
              if (nrow(df_anova_subset)>0){
                # In case the term does not exist in df_gamm, plot using df_anova instead
                df_stat<-df_anova_subset
              }
            }
            if (nrow(df_stat)>0){
              if (is.na(df_stat[1,"F"])){ # GAMM result
                list_df_stat<-list("both"=df_stat,"pos"=df_stat[df_stat$estimate>0,],"neg"=df_stat[df_stat$estimate<0,])
              }else{ # ANCOVA result
                list_df_stat<-list("both"=df_stat,"pos"=data.frame(),"neg"=data.frame())
              }
            }else{
              list_df_stat<-list("both"=data.frame(),"pos"=data.frame(),"neg"=data.frame())
            }
            for (type_sign in names(list_df_stat)){
              #print(paste(idx_mod,var_exp_detect,label_sex,type_sign))
              df_stat_temp<-list_df_stat[[type_sign]]
              if (nrow(df_stat_temp)>0){
                df_head<-data.frame(model=idx_mod,term_perm=idx_term_perm,term_detect=idx_term_detect,sex=idx_sex,sign=type_sign)
                df_head_label<-data.frame(term_perm=var_exp_perm,term_detect=var_exp_detect,sex=label_sex)
                if (is.null(df_deltah_in)){
                  delta_h_in<-NULL
                  list_src_tfnbs<-c(list_src_tfnbs,
                                    list(list("df_stat"=df_stat_temp,"delta_h_in"=delta_h_in,"df_head"=df_head,"df_head_label"=df_head_label)))
                }else{
                  # if delta_h_in is not NULL (= calculating permutation),
                  # calculate tfNBS only for the set of parameters on delta_h_in
                  if (var_exp_perm==var_exp_perm_in){
                    delta_h_in<-df_deltah_in[df_deltah_in$model==idx_mod & df_deltah_in$term_perm==idx_term_perm & df_deltah_in$term_detect==idx_term_detect
                                             & df_deltah_in$sex==idx_sex & df_deltah_in$sign==type_sign,"delta_h"]
                    list_src_tfnbs<-c(list_src_tfnbs,
                                      list(list("df_stat"=df_stat_temp,"delta_h_in"=delta_h_in,"df_head"=df_head,"df_head_label"=df_head_label)))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  # parallel computing
  list_dst_tfnbs<-parLapply(clust,list_src_tfnbs,func_tfnbs3)
  
  # combine results
  df_tfnbs<-data.frame(); df_max<-data.frame(); df_deltah<-data.frame(); list_plot<-list()
  for (dst_tfnbs in list_dst_tfnbs){
    df_head<-dst_tfnbs$df_head
    df_head_label<-dst_tfnbs$df_head_label
    if (return_nbs){
      df_tfnbs<-rbind(df_tfnbs,cbind(df_head,dst_tfnbs$df_tfnbs))
    }
    df_max<-rbind(df_max,data.frame(df_head,max_nbs=dst_tfnbs$max_nbs))
    df_deltah<-rbind(df_deltah,data.frame(df_head,delta_h=dst_tfnbs$delta_h))
    if (plot_result){
      plot_nbs<-plot_tfnbs(dst_tfnbs$df_tfnbs,data_fc)
      plot_nbs<-(plot_nbs+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",df_head_label$term_detect,", sex: ",df_head_label$sex,", sign: ",df_head$sign,sep="")))
      list_plot<-c(list_plot,list(list("plot"=plot_nbs,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                       "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",df_head$model,"_trm-",df_head$term_detect,"_sex-",df_head_label$sex,"_sgn-",df_head$sign,"_tfnbs.png",sep=""))))
    }
  }
  
  if (length(list_plot)>0){
    nullobj<-parLapply(clust,list_plot,plot_parallel_core)
  }
  if (return_nbs){
    fwrite(df_tfnbs,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_tfnbs.csv",sep="")),row.names = F)
  }
  return(list("df_tfnbs"=df_tfnbs,"df_max"=df_max,"df_deltah"=df_deltah))
}


func_iterate_tfnbs<-function(paths,clust,df_gamm,df_anova,df_deltah_in,data_fc,plot_result=F,return_nbs=T,
                             atlas,param,list_mod,list_term,idx_var,label_wave){
  list_plot<-list()
  df_tfnbs<-data.frame(); df_max<-data.frame(); df_deltah<-data.frame()
  for (idx_mod in param$param_nbs$list_mod){
    for (set_term in param$param_nbs$list_term){
      for (idx_term_detect in set_term$term_detect){
        var_exp_detect<-list_term[[idx_term_detect]]$var_exp
        if (!is.null(var_exp_detect)){
          for (idx_sex in param$list_sex){
            if (idx_sex==1){label_sex<-"m"}else if (idx_sex==2){label_sex<-"f"}else{label_sex<-"mf"}
            df_stat<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp_detect & df_gamm$sex==idx_sex,]
            if (nrow(df_stat)==0){
              df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp_detect & df_anova$sex==idx_sex,]
              if (nrow(df_anova_subset)>0){
                # In case the term does not exist in df_gamm, plot using df_anova instead
                df_stat<-df_anova_subset
              }
            }
            if (nrow(df_stat)>0){
              if (is.na(df_stat[1,"F"])){ # GAMM result
                list_df_stat<-list("both"=df_stat,"pos"=df_stat[df_stat$estimate>0,],"neg"=df_stat[df_stat$estimate<0,])
              }else{ # ANCOVA result
                list_df_stat<-list("both"=df_stat,"pos"=data.frame(),"neg"=data.frame())
              }
            }else{
              list_df_stat<-list("both"=data.frame(),"pos"=data.frame(),"neg"=data.frame())
            }
            for (type_sign in names(list_df_stat)){
              #print(paste(idx_mod,var_exp_detect,label_sex,type_sign))
              df_stat_temp<-list_df_stat[[type_sign]]
              if (nrow(df_stat_temp)>0){
                if (is.null(df_deltah_in)){
                  delta_h_in<-NULL
                }else{
                  delta_h_in<-df_deltah_in[df_deltah_in$model==idx_mod & df_deltah_in$term==var_exp_detect
                                           & df_deltah_in$sex==idx_sex & df_deltah_in$sign==type_sign,"delta_h"]
                }
                data_tfnbs<-func_tfnbs2(df_stat_temp,delta_h_in,param,clust)
                df_head<-data.frame(model=idx_mod,term=var_exp_detect,sex=idx_sex,sign=type_sign)
                if (return_nbs){
                  df_tfnbs<-rbind(df_tfnbs,cbind(df_head,data_tfnbs$df_tfnbs))
                }
                df_max<-rbind(df_max,data.frame(df_head,max_nbs=data_tfnbs$max_nbs))
                df_deltah<-rbind(df_deltah,data.frame(df_head,delta_h=data_tfnbs$delta_h))
                if (plot_result){
                  plot_nbs<-plot_tfnbs(data_tfnbs$df_tfnbs,data_fc)
                  plot_nbs<-(plot_nbs+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp_detect,", sex: ",label_sex,", sign: ",type_sign,sep="")))
                  list_plot<-c(list_plot,list(list("plot"=plot_nbs,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                   "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_sgn-",type_sign,"_tfnbs.png",sep=""))))
                }
              }
            }
          }
        }
      }
    }
  }
  if (length(list_plot)>0){
    clust<-makeCluster(floor(detectCores()*3/4))
    plot_parallel(clust,list_plot)
    stopCluster(clust)
  }
  if (return_nbs){
    fwrite(df_tfnbs,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_tfnbs.csv",sep="")),row.names = F)
  }
  return(list("df_tfnbs"=df_tfnbs,"df_max"=df_max,"df_deltah"=df_deltah))
}

func_tfnbs_permutation<-function(paths,clust,data_fc,df_clin,data_tfnbs_in,calc_parallel,plot_result=T,
                                 atlas,param,list_mod,list_term,idx_var,label_wave){
  
  if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_fwep.csv",sep="")))){
    print("Calculated permutation exists")
  }else{
    print("Calculating permutation")
    set.seed(0)
    pb<-txtProgressBar(min=0,max=param$param_nbs$n_perm,style=3,width=50)
    df_max_nbs<-data.frame();df_deltah<-data.frame()
    for (idx_perm in seq(param$param_nbs$n_perm)){
      for (set_term in param$param_nbs$list_term){
        var_exp_perm<-list_term[[set_term$term_perm]]$var_exp
        if (!is.null(var_exp_perm)){
          # Sex-wise permutation of term (expvar) of interst
          df_clin_perm<-NULL
          for (idx_sex in param$list_sex){
            df_clin_perm_add<-df_clin[df_clin$sex==idx_sex,]
            df_clin_perm_add[,var_exp_perm]<-sample(df_clin_perm_add[,var_exp_perm])
            df_clin_perm<-rbind(df_clin_perm,df_clin_perm_add)
          }
          # Join FC and permuted clinical data
          df_join<-join_fc_clin(data_fc$df_fc,df_clin_perm)
          
          # Calculate model
          data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=F)
          df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova
          # Calculate threshold-free network-based statistics
          data_tfnbs<-func_iterate_tfnbs2(paths,clust,df_gamm,df_anova,df_deltah_in=data_tfnbs_in$df_deltah,var_exp_perm_in=var_exp_perm,data_fc,plot_result=F,return_nbs=F,
                                         atlas,param,list_mod,list_term,idx_var,label_wave)
          df_max_nbs<-rbind(df_max_nbs,data.frame(id_perm=idx_perm,data_tfnbs$df_max))
          df_deltah<-rbind(df_deltah,data.frame(id_perm=idx_perm,data_tfnbs$df_deltah))
        }
      }
      setTxtProgressBar(pb,idx_perm)
    }
    close(pb)
    fwrite(df_max_nbs,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_max.csv",sep="")),row.names = F)
    fwrite(df_deltah,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_deltah.csv",sep="")),row.names = F)
    
    # Summarize permutation result
    list_plot<-list()
    df_thresh_nbs<-df_fwep<-df_sign<-data.frame()
    for (idx_mod in param$param_nbs$list_mod){
      for (set_term in param$param_nbs$list_term){
        idx_term_perm<-set_term$term_perm
        var_exp_perm<-list_term[[idx_term_perm]]$var_exp
        for (idx_term_detect in set_term$term_detect){
          var_exp_detect<-list_term[[idx_term_detect]][["var_exp"]]
          for (idx_sex in param$list_sex){
            if (idx_sex==1){
              label_sex<-"m";title_sex<-"male";color_plt<-"steelblue2"
            }else if (idx_sex==2){
              label_sex<-"f";title_sex<-"female";color_plt<-"lightcoral"
            }else{
              label_sex<-"mf";title_sex<-"both";color_plt<-"grey60"
            }
            
            for (type_sign in c("both","pos","neg")){
              list_max_nbs<-df_max_nbs[df_max_nbs$model==idx_mod & df_max_nbs$term_perm==idx_term_perm & df_max_nbs$term_detect==idx_term_detect & df_max_nbs$sex==idx_sex & df_max_nbs$sign==type_sign, "max_nbs"]
              list_max_nbs<-sort(list_max_nbs)
              if (length(list_max_nbs)>0){
                df_tfnbs<-data_tfnbs_in$df_tfnbs
                df_tfnbs_fwep<-df_tfnbs[df_tfnbs$model==idx_mod & df_tfnbs$term_perm==idx_term_perm & df_tfnbs$term_detect==idx_term_detect & df_tfnbs$sex==idx_sex
                                        & df_tfnbs$sign==type_sign,]
                for (idx_row in seq(nrow(df_tfnbs_fwep))){
                  df_tfnbs_fwep[idx_row,"fwep"]<-sum(df_tfnbs_fwep[idx_row,"nbs"]<list_max_nbs)/param$param_nbs$n_perm
                }
                df_tfnbs_sign<-df_tfnbs_fwep[df_tfnbs_fwep$fwep<=param$param_nbs$p_perm_threshold,]
                df_head<-data.frame(model=idx_mod,term_perm=idx_term_perm,term_detect=idx_term_detect,sex=idx_sex,sign=type_sign)
                df_fwep<-rbind(df_fwep,df_tfnbs_fwep)
                if (nrow(df_tfnbs_sign)>0){
                  df_sign<-rbind(df_sign,df_tfnbs_sign)
                }
                thr_nbs<-list_max_nbs[ceiling(length(list_max_nbs)*(1-param$param_nbs$p_perm_threshold))]
                df_thresh_nbs<-rbind(df_thresh_nbs,data.frame(df_head,thresh_nbs=thr_nbs))
                if (plot_result){
                  plot_nbs<-plot_tfnbs(df_tfnbs_sign,data_fc)
                  plot_nbs<-plot_nbs+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp_detect,", sex: ",label_sex,", sign: ",type_sign,sep=""))
                  list_plot<-c(list_plot,list(list("plot"=plot_nbs,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                   "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_sgn-",type_sign,"_tfnbs_sign.png",sep=""))))
                  #plot_perm<-plot_permutation2(list_max=list_max_nbs,thr_nbs=thr_nbs,color_plt)
                  plot_perm<-plot_permutation2(list_max=list_max_nbs,color_plt)
                  plot_perm<-plot_perm+ ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,
                                                      "\nexpvar: ",var_exp_detect,", sex: ",label_sex,", sign: ",type_sign,sep=""))
                  list_plot<-c(list_plot,list(list("plot"=plot_perm,"height"=5,"width"=7,"dpi"=300,"path"=file.path(paths$output,"output","plot"),
                                                   "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_sgn-",type_sign,"_perm.png",sep=""))))
                }
              }
            }
          }
        }
      }
    }
    if (length(list_plot)>0){
      nullobj<-parLapply(clust,list_plot,plot_parallel_core)
    }
    fwrite(df_thresh_nbs,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_thr.csv",sep="")),row.names = F)
    fwrite(df_sign,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_sign.csv",sep="")),row.names = F)
    fwrite(df_fwep,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_fwep.csv",sep="")),row.names = F)
  }
}


func_detect_subnet<-function(paths,param,data_plot,data_gamm,data_fc,
                             atlas,list_covar,list_mod,list_term,idx_var,label_wave,plot_result=F){
  df_plot<-data_plot$df_plot
  df_gamm<-data_gamm$df_gamm
  if (nrow(df_plot)>0){
    print("Detecting subnetworks")
    list_plot<-list()
    df_net<-df_node<-df_size_net<-df_pred_ancova<-data.frame()
    #list_output<-list()
    for (idx_mod in names(list_mod)){
      for (set_term in param$param_nbs$list_term){
        for (idx_term_detect in set_term$term_detect){
          var_exp_detect<-list_term[[idx_term_detect]]$var_exp
          if (!is.null(var_exp_detect)){
            for (idx_sex in param$list_sex){
              if (idx_sex==1){label_sex<-"m"}else if (idx_sex==2){label_sex<-"f"}else{label_sex<-"mf"}
              for (p_cdt in param$list_p$p){
                df_sign<-df_plot[df_plot$p_type=="p" & df_plot$p_threshold==p_cdt
                                 & df_plot$model==idx_mod & df_plot$term==var_exp_detect & df_plot$sex==idx_sex,]
                if (nrow(df_sign)>0){
                  if (is.na(df_sign[1,"F"])){ # GAMM result
                    list_df_sign<-list("both"=df_sign,"pos"=df_sign[df_sign$estimate>0,],"neg"=df_sign[df_sign$estimate<0,])
                  }else{ # ANCOVA result
                    list_df_sign<-list("both"=df_sign,"pos"=data.frame(),"neg"=data.frame())
                  }
                }else{
                  list_df_sign<-list("both"=data.frame(),"pos"=data.frame(),"neg"=data.frame())
                }
                
                for (type_sign in names(list_df_sign)){
                  df_sign_temp<-list_df_sign[[type_sign]]
                  if (nrow(df_sign_temp)>0){
                    data_bfs<-func_bfs(df_sign_temp)
                    if(length(data_bfs$list_network)>0){
                      for (idx_net in seq(length(data_bfs$list_network))){
                        network<-data_bfs$list_network[[idx_net]]
                        df_head<-data.frame(model=idx_mod,term=var_exp_detect,sex=idx_sex,sign=type_sign)
                        df_net<-rbind(df_net,data.frame(sign=type_sign,id_net=idx_net,network$df_edge))
                        df_node_add<-inner_join(data.frame(p_threshold=p_cdt,id_net=idx_net,network$df_node),data_fc$df_roi,by=c("node"="id"))
                        df_node_add<-dplyr::rename(df_node_add,"label_node"="label","group_node"="group")
                        df_node<-rbind(df_node,cbind(df_head,df_node_add))
                        df_size_net<-rbind(df_size_net,cbind(df_head,data.frame(p_threshold=p_cdt,id_net=idx_net,size=network$size_net)))
                        if (plot_result){
                          plot_subnet<-plot_net(df_edge=network$df_edge,df_node=network$df_node,df_roi=data_fc$df_roi)
                          plot_subnet<-(plot_subnet+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", sign: ",type_sign,", #",as.character(idx_net),sep="")))
                          list_plot<-c(list_plot,list(list("plot"=plot_subnet,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                           "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_sgn-",type_sign,"_idx-",as.character(idx_net),"_subnet.png",sep=""))))
                        }  
                        if(idx_term_detect %in% names(param$param_ancova_pred)){
                          df_pred_ancova_add<-func_pred_ancova(df_edge=network$df_edge,df_gamm=df_gamm,data_fc=data_fc,param_ancova_pred=param$param_ancova_pred,idx_term_detect,var_exp_detect)
                          df_pred_ancova<-rbind(df_pred_ancova,data.frame(sign=type_sign,p_threshold=p_cdt,id_net=idx_net,df_pred_ancova_add))
                          if (plot_result){
                            plot_pred<-(plot_pred_ancova(df_pred_ancova_add)+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,"\nexpvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", sign: ",type_sign,", #",as.character(idx_net),sep=""))+ xlab(list_term[[idx_term_detect]][["title"]]))
                            list_plot<-c(list_plot,list(list("plot"=plot_pred,"height"=5,"width"=5,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                             "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_sgn-",type_sign,"_idx-",as.character(idx_net),"_pred.png",sep=""))))
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (length(list_plot)>0){
      clust<-makeCluster(floor(detectCores()*3/4))
      plot_parallel(clust,list_plot)
      stopCluster(clust)
    }
    fwrite(df_net,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_bfs_edge.csv",sep="")),row.names = F)
    fwrite(df_node,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_bfs_node.csv",sep="")),row.names = F)
    fwrite(df_size_net,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_bfs_size.csv",sep="")),row.names = F)
    if (nrow(df_pred_ancova)>0){
      fwrite(df_pred_ancova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_bfs_pred.csv",sep="")),row.names = F)
    }
  }
  return(list("df_net"=df_net,"df_node"=df_node,"df_size_net"=df_size_net,"df_pred_ancova"=df_pred_ancova))
}

func_nbs_permutation<-function(paths,param,df_clin,data_bfs,data_fc, atlas,list_covar,list_mod,list_term,idx_var,label_wave,calc_parallel,plot_result=T){
  #paths,df_fc,df_clin,data_bfs,data_fc,calc_parallel,plot_result=T,
                               #atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave){
  
  if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_fwep.csv",sep="")))){
    print("Calculated permutation exists")
  }else{
    print("Calculating permutation")
    
    # Prepare parallelization cluster
    test_mod<-F
    if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
    list_sex<-param$list_sex
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                                  "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC"),
                  envir=environment())
    set.seed(0)
    pb<-txtProgressBar(min=0,max=param$param_nbs$n_perm,style=3,width=50)
    df_max_size<-data.frame()
    for (idx_perm in seq(param$param_nbs$n_perm)){
      for (set_term in param$param_nbs$list_term){
        var_exp_perm<-list_term[[set_term$term_perm]]$var_exp
        if (!is.null(var_exp_perm)){
          # Sex-wise permutation of term (expvar) of interst
          df_clin_perm<-NULL
          for (idx_sex in param$list_sex){
            df_clin_perm_add<-df_clin[df_clin$sex==idx_sex,]
            df_clin_perm_add[,var_exp_perm]<-sample(df_clin_perm_add[,var_exp_perm])
            df_clin_perm<-rbind(df_clin_perm,df_clin_perm_add)
          }
          # Join FC and permuted clinical data
          df_join<-join_fc_clin(data_fc$df_fc,df_clin_perm)
          
          # Calculate model
          data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
          df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova
          for (idx_mod in names(list_mod)){
            list_term_detect<-set_term$term_detect
            for (idx_term_detect in list_term_detect){
              var_exp_detect<-list_term[[idx_term_detect]][["var_exp"]]
              for (idx_sex in param$list_sex){
                # Subset GAMM result dataframe for plotting
                df_gamm_subset<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp_detect & df_gamm$sex==idx_sex,]
                if (nrow(df_gamm_subset)==0){
                  df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp_detect & df_anova$sex==idx_sex,]
                  if (nrow(df_anova_subset)>0){
                    # In case the term does not exist in df_gamm, plot using df_anova instead
                    df_gamm_subset<-df_anova_subset
                  }
                }
                if (nrow(df_gamm_subset)>0){ # If the model/expvar/sex exist either in df_gamm or df_anova
                  for (p_cdt in param$list_p$p){
                    df_sign<-df_gamm_subset[df_gamm_subset$p<p_cdt,]
                    if (nrow(df_sign)>0){
                      if (is.na(df_sign[1,"F"])){ # GAMM result
                        list_df_sign<-list("both"=df_sign,"pos"=df_sign[df_sign$estimate>0,],"neg"=df_sign[df_sign$estimate<0,])
                      }else{ # ANCOVA result
                        list_df_sign<-list("both"=df_sign,"pos"=data.frame(),"neg"=data.frame())
                      }
                    }else{
                      list_df_sign<-list("both"=data.frame(),"pos"=data.frame(),"neg"=data.frame())
                    }
                    
                    for (type_sign in names(list_df_sign)){
                      df_sign_temp<-list_df_sign[[type_sign]]
                      if (nrow(df_sign_temp)>0){
                        max_size<-func_bfs(df_sign_temp)$max_size
                      }else{
                        max_size<-0
                      }
                      df_max_size<-rbind(df_max_size,data.frame(id_perm=idx_perm,model=idx_mod,term=var_exp_detect,sex=idx_sex,p_threshold=p_cdt,sign=type_sign,max_size=max_size))
                    }
                  }
                }
              }
            }
          }
        }
      }
      setTxtProgressBar(pb,idx_perm)
    }
    stopCluster(clust)
    close(pb)
    fwrite(df_max_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_max.csv",sep="")),row.names = F)
    
    # Summarize permutation result
    df_net<-data_bfs$df_net; df_node<-data_bfs$df_node; df_size_net<-data_bfs$df_size_net; df_pred_ancova<-data_bfs$df_pred_ancova
    list_plot<-list()
    df_threshold_size<-df_fwep<-NULL
    for (idx_mod in names(list_mod)){
      for (set_term in param$param_nbs$list_term){
        for (idx_term_detect in set_term$term_detect){
          var_exp_detect<-list_term[[idx_term_detect]][["var_exp"]]
          for (idx_sex in param$list_sex){
            if (idx_sex==1){
              label_sex<-"m";title_sex<-"male"
            }else if (idx_sex==2){
              label_sex<-"f";title_sex<-"female"
            }else{
              label_sex<-"mf";title_sex<-"both"
            }
            for (p_cdt in param$list_p$p){
              for (type_sign in c("both","pos","neg")){
                list_max_size<-df_max_size[df_max_size$model==idx_mod & df_max_size$term==var_exp_detect
                                           & df_max_size$p_threshold==p_cdt & df_max_size$sex==idx_sex & df_max_size$sign==type_sign,
                                           "max_size"]
                if (length(list_max_size)>0){
                  list_max_size<-sort(list_max_size)
                  df_size_net_subset<-df_size_net[df_size_net$model==idx_mod & df_size_net$term==var_exp_detect
                                                  & df_size_net$p_threshold==p_cdt & df_size_net$sex==idx_sex & df_size_net$sign==type_sign,]
                  if(nrow(df_size_net_subset)>0){
                    for (idx_row in seq(nrow(df_size_net_subset))){
                      p_fwe<-sum(list_max_size>df_size_net_subset[idx_row,"size"])/param$param_nbs$n_perm
                      df_size_net_subset[idx_row,"p_fwe"]<-p_fwe
                      if (p_fwe<param$param_nbs$p_perm_threshold){
                        if (plot_result){
                          idx_net<-as.character(df_size_net_subset[idx_row,"id_net"])
                          df_net_subset<-df_net[df_net$model==idx_mod & df_net$term==var_exp_detect & df_net$p_threshold==p_cdt & df_net$sex==idx_sex & df_net$sign==type_sign & df_net$id_net==idx_net,]
                          df_node_subset<-df_node[df_node$model==idx_mod & df_node$term==var_exp_detect & df_node$p_threshold==p_cdt & df_node$sex==idx_sex & df_node$sign==type_sign & df_node$id_net==idx_net,
                                                  c("node","degree","label_node","group_node")]
                          plot_subnet<-plot_net(df_edge=df_net_subset,df_node=df_node_subset,df_roi=data_fc$df_roi)
                          plot_subnet<-(plot_subnet+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", sign: ",type_sign,", #",as.character(idx_net),sep="")))
                          list_plot<-c(list_plot,list(list("plot"=plot_subnet,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                           "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_sgn-",type_sign,"_idx-",as.character(idx_net),"_subnet.png",sep=""))))
                          if(idx_term_detect %in% names(param$param_ancova_pred)){
                            df_pred_ancova_subset<-df_pred_ancova[df_pred_ancova$model==idx_mod & df_pred_ancova$term==var_exp_detect & df_pred_ancova$p_threshold==p_cdt & df_pred_ancova$sex==idx_sex & df_pred_ancova$sign==type_sign & df_pred_ancova$id_net==idx_net,]
                            plot_pred<-(plot_pred_ancova(df_pred_ancova_subset)+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,"\nexpvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", sign: ",type_sign,", #",as.character(idx_net),sep=""))+ xlab(list_term[[idx_term_detect]][["title"]]))
                            list_plot<-c(list_plot,list(list("plot"=plot_pred,"height"=5,"width"=5,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                             "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_sgn-",type_sign,"_idx-",as.character(idx_net),"_pred.png",sep=""))))
                          }
                        }  
                      }
                    }
                    df_fwep<-rbind(df_fwep,df_size_net_subset)
                  }
                  thr_size_nbs<-list_max_size[ceiling(length(list_max_size)*(1-param$param_nbs$p_perm_threshold))]
                  df_threshold_size<-rbind(df_threshold_size,data.frame(model=idx_mod,term=var_exp_detect,sex=idx_sex,p_threshold=p_cdt,sign=type_sign,
                                                                        thr_size=thr_size_nbs))
                  
                  title_plot<-list_term[[idx_term_detect]][["title"]]
                  plot_permutation(paths,list_max=list_max_size,thr_size_nbs,
                                   atlas,var=idx_var,wave=label_wave,idx_mod,idx_term_detect,label_sex,title_plot,title_sex,p_cdt,type_sign)
                }
              }
            }
          }
        }
      }
    }
    if (length(list_plot)>0){
      clust<-makeCluster(floor(detectCores()*3/4))
      plot_parallel(clust,list_plot)
      stopCluster(clust)
    }
    fwrite(df_threshold_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_thr.csv",sep="")),row.names = F)
    fwrite(df_fwep,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_perm_fwep.csv",sep="")),row.names = F)
  }
}

func_pred_ancova<-function(df_edge,df_gamm,data_fc,param_ancova_pred,idx_term,var_exp){
  df_edge<-df_edge[,c("from","to","sex","model")]
  #df_edge$sex<-as.character(df_edge$sex)
  df_gamm<-df_gamm[df_gamm$term %in% param_ancova_pred[[idx_term]]$term,c("from","to","sex","model","term","estimate")]
  df_gamm<-inner_join(df_edge,df_gamm,by=c("from","to","sex","model"))
  df_gamm<-inner_join(df_gamm,data_fc$df_edge,by=c("from","to"))
  df_plot<-NULL
  for (id_edge in unique(df_gamm$id_edge)){
    df_plot_add<-df_gamm[df_gamm$id_edge==id_edge,]
    term_base<-param_ancova_pred[[idx_term]][1,"term"]
    df_plot_add[df_plot_add$term!=term_base,"estimate"]<-df_plot_add[df_plot_add$term!=term_base,"estimate"]+df_plot_add[df_plot_add$term==term_base,"estimate"]
    df_plot<-rbind(df_plot,df_plot_add)
  }
  df_plot<-inner_join(df_plot,param_ancova_pred[[idx_term]],by="term")
  df_plot$term<-var_exp
  df_plot$label_edge<-paste(df_plot$label_from,df_plot$label_to,sep=" - ")
  df_plot<-df_plot[,c("from","to","label_from","label_to","label_edge","sex","model","term","level","estimate")]
  
  return(df_plot)
}


#**************************************************
# Sex difference of FC ============================
#**************************************************

func_nbs<-function(paths,atlas,wave,df_fc,df_clin,list_mod,calc_slope,list_plot,list_sex,
                   df_roi,df_edge,df_grp,thr_p_cdt,n_perm,thr_p_perm,calc_parallel,test_mod=F){
  print(paste("Calculating model, atlas: ",atlas,", wave: ",wave,sep=""))
  if (calc_parallel){
    clust<-makeCluster(floor(detectCores()*3/4))
  }else{
    clust<-makeCluster(1)
  }
  clusterExport(clust,
                varlist=c("list_mod","list_sex","calc_parallel","test_mod","sort","gam","as.formula","summary.gam",
                          "anova.gam","as.numeric.factor"),
                envir=environment())
  data_nbs<-func_nbs_core(clust=clust,df_fc=df_fc,df_clin=df_clin,
                          df_roi=df_roi,df_edge=df_edge,list_mod=list_mod,
                          thr_p_cdt=thr_p_cdt,list_plot=list_plot,
                          progressbar=F,output_gamm=F,calc_slope=calc_slope,test_mod=test_mod)
  if(test_mod){
    return(data_nbs)
  }else{
    data_nbs<-data_nbs$data_nbs
    
    # Permutation test
    print(paste("Calculating permutation, atlas: ",atlas,", wave: ",wave,sep=""))
    set.seed(0)
    list_max<-list()
    pb<-txtProgressBar(min=0,max=n_perm,style=3,width=50)
    for (idx_perm in seq(n_perm)){
      df_clin_perm<-df_clin
      df_clin_perm$sex<-sample(df_clin_perm$sex)
      #print(as.character(idx_perm))
      data_nbs_perm<-func_nbs_core(clust=clust,df_fc=df_fc,df_clin=df_clin_perm,
                                   df_roi=df_roi,df_edge=df_edge,list_mod=list_mod,
                                   thr_p_cdt=thr_p_cdt,list_plot=list_plot,
                                   progressbar=F,output_gamm=F,calc_slope=calc_slope,test_mod=F)$data_nbs
      for (model in names(data_nbs_perm)){
        for (plot in names(data_nbs_perm[[model]])){
          if(idx_perm==1){
            list_max_sex<-list("m"=data_nbs_perm[[model]][[plot]][["m"]][["max_size"]],
                               "f"=data_nbs_perm[[model]][[plot]][["f"]][["max_size"]])
            list_max[[model]][[plot]]<-list_max_sex
          }else{
            for (sex in c("m","f")){
              list_max[[model]][[plot]][[sex]]<-c(list_max[[model]][[plot]][[sex]],
                                                  data_nbs_perm[[model]][[plot]][[sex]][["max_size"]])
            }
          }
        }
      }
      setTxtProgressBar(pb,idx_perm)
    }
    close(pb)
    
    # Summarize permutation and threshold subgraphs
    print(paste("Preparing output, atlas: ",atlas,", wave: ",wave,sep=""))
    df_net<-df_size_net<-df_perm<-df_thr_size<-NULL
    list_output<-list()
    for (model in names(data_nbs)){
      for (plot in names(data_nbs[[model]])){
        for (sex in c("m","f")){
          data_nbs_subset<-data_nbs[[model]][[plot]][[sex]]
          list_max_subset<-list_max[[model]][[plot]][[sex]]
          list_max_subset_sort<-sort(list_max_subset)
          thr_size_perm<-list_max_subset_sort[ceiling(length(list_max_subset_sort)*(1-thr_p_perm))]
          
          if (sex=="m"){
            title_sex<-"M>F"
            color_plt<-"steelblue2"
          }else{
            title_sex<-"F>M"
            color_plt<-"lightcoral"
          }
          title_plot<-list_plot[[plot]][["title"]]
          #list_output<-c(list_output,
          #               list(plot_permutation(paths,list_max=list_max_subset_sort,thr_size_perm,
          #                                     atlas,wave,model,plot,sex,title_plot,title_sex,color_plt)))
          plot_permutation(paths,list_max=list_max_subset_sort,thr_size_perm,
                           atlas,var="sex",wave,model,plot,sex,title_plot,title_sex,color_plt)
          df_head<-data.frame(atlas=atlas,wave=wave,mod=model,plot=plot,sex=sex)
          list_network_sign<-list()
          if(length(data_nbs_subset$list_network)>0){
            for (idx_net in seq(length(data_nbs_subset$list_network))){
              network<-data_nbs_subset$list_network[[idx_net]]
              list_output<-c(list_output,
                             list(plot_sex_diff_fc(paths,network$df_edge,atlas,wave,df_roi,df_grp,
                                                   model,plot,sex,title_plot,title_sex,idx_net)))
              df_net<-rbind(df_net,cbind(df_head, data.frame(id_net=idx_net,network$df_edge)))
              df_size_net<-rbind(df_size_net,cbind(df_head,data.frame(id_net=idx_net,size=network$size_net)))
              if (network$size_net>=thr_size_perm){
                list_network_sign<-c(list(network))
              }
            }
          }
          data_nbs[[model]][[plot]][[sex]][["list_network_sign_perm"]]<-list_network_sign
          data_nbs[[model]][[plot]][[sex]][["list_max_size_perm"]]<-list_max_subset
          data_nbs[[model]][[plot]][[sex]][["thr_size_perm"]]<-thr_size_perm
          df_thr_size<-rbind(df_thr_size,
                             cbind(df_head,data.frame(thr_size=thr_size_perm)))
          df_perm<-rbind(df_perm,
                         cbind(df_head,data.frame(id_perm=seq(length(list_max_subset)),
                                                  max_size_net=list_max_subset)))
        }
      }
    }
    plot_parallel(clust,list_output)
    stopCluster(clust)
    write.csv(df_net,file.path(paths$output,"output","temp",
                               paste("atl-",atlas,"_wave-",wave,"_net.csv",sep="")),row.names=F)
    write.csv(df_size_net,file.path(paths$output,"output","temp",
                                    paste("atl-",atlas,"_wave-",wave,"_size_net.csv",sep="")),row.names=F)
    write.csv(df_thr_size,file.path(paths$output,"output","temp",
                                    paste("atl-",atlas,"_wave-",wave,"_thr_perm.csv",sep="")),row.names=F)
    write.csv(df_perm,file.path(paths$output,"output","temp",
                                paste("atl-",atlas,"_wave-",wave,"_perm.csv",sep="")),row.names=F)
  }
}

sex_diff_fc<-function(paths_=paths,list_atlas_=list_atlas,key_group_='group_3',
                      subset_subj_=sex_diff_fc_subset_subj,list_covar_=sex_diff_fc_list_covar,
                      list_mod_diff_=sex_diff_fc_list_mod_diff,
                      list_mod_long_=sex_diff_fc_list_mod_long,
                      list_mod_cs_=sex_diff_fc_list_mod_cs,
                      list_plot_=sex_diff_fc_list_plot,
                      thr_p_cdt_=sex_diff_fc_thr_p_cdt,thr_p_perm_=sex_diff_fc_thr_p_perm,
                      n_perm_=sex_diff_fc_n_perm){
  print("Starting sex_diff_fc()")
  nullobj<-func_createdirs(paths_,str_proc="sex_diff_fc()",copy_log=T)
  # Increase memory limit
  memory.limit(1000000)
  for (atlas in list_atlas_){
    #if (!file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_wave-",wave,"_perm.csv",sep="")))){
    if (!file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_wave-diff_perm.csv",sep="")))){
      print(paste("Preparing data: ",atlas,sep=""))
      data_fc<-prep_data_fc(paths_,atlas,key_group_,include_diff=T,include_grp=F)
      data_clin<-func_clinical_data_long(paths_,list_wave=c("1","2"),subset_subj_,list_covar_,rem_na_clin=T,
                                         prefix=paste("atl-",atlas,sep=""),print_terminal=F)
      
      # Longitudinal change analysis
      df_fc_diff<-data_fc$df_fc[data_fc$df_fc$ses=="2-1",]
      # List of subjects meeting QC criteria and has non-NA covariates in both waves
      list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
      df_clin_diff<-data_clin$df_clin
      colnames(df_clin_diff)[colnames(df_clin_diff)=="wave"]<-"ses"
      df_clin_diff<-func_clinical_data_diffmean(df_clin_diff,list_id_subj,list_covar_)
      df_clin_diff$wave<-"2-1"
      #df_clin_diff<-func_demean_clin(df_clin_diff,thr_cont=10)$df_clin
      mean_mean_age<-mean(df_clin_diff$mean_age)
      df_clin_diff$mean_age<-df_clin_diff$mean_age-mean_mean_age
      
      # Network-based statistics
      func_nbs(paths=paths_,atlas=atlas,wave="diff",
               df_fc=df_fc_diff,df_clin=df_clin_diff,
               list_mod=list_mod_diff_,calc_slope=T,list_plot=list_plot_,list_sex=list(c(1,2)),
               df_roi=data_fc$df_roi,df_edge=data_fc$df_edge,df_grp=data_fc$df_grp,
               thr_p_cdt=thr_p_cdt_,n_perm=n_perm_,thr_p_perm=thr_p_perm_,
               calc_parallel=T,test_mod=F)

      # Generalized linear(/additive) mixed model analysis
      #df_fc_long<-data_fc$df_fc[data_fc$df_fc$ses %in% c(1,2),]
      #df_clin_long<-data_clin$df_clin

      
      # Cross-sectional analysis
      
    }
  }
  
  print("Binding results")
  df_net<-df_size_net<-df_perm<-df_thr_size<-NULL
  for (atlas in list_atlas_){
    for (wave in c("diff")){
      df_net<-rbind(df_net,
                    as.data.frame(fread(file.path(paths$output,"output","temp",
                                                  paste("atl-",atlas,"_wave-",wave,"_net.csv",sep="")))))
      df_size_net<-rbind(df_size_net,
                         as.data.frame(fread(file.path(paths$output,"output","temp",
                                                       paste("atl-",atlas,"_wave-",wave,"_size_net.csv",sep="")))))
      df_thr_size<-rbind(df_thr_size,
                         as.data.frame(fread(file.path(paths$output,"output","temp",
                                                       paste("atl-",atlas,"_wave-",wave,"_thr_perm.csv",sep="")))))
      df_perm<-rbind(df_perm,
                     as.data.frame(fread(file.path(paths$output,"output","temp",
                                                   paste("atl-",atlas,"_wave-",wave,"_perm.csv",sep="")))))
    }
  }
  write.csv(df_net,file.path(paths$output,"output","result","net.csv"),row.names=F)
  write.csv(df_size_net,file.path(paths$output,"output","result","size_net.csv"),row.names=F)
  write.csv(df_thr_size,file.path(paths$output,"output","result","thr_perm.csv"),row.names=F)
  write.csv(df_perm,file.path(paths$output,"output","result","perm.csv"),row.names=F)
  print("Finished sex_diff_fc()")
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


#**************************************************
# OBSOLETE ========================================
#**************************************************
#gamm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
#                  list_wave_=list_wave,list_atlas_=list_atlas,
#                  #list_measure_=list_measure,list_str_group_=list_str_group,
#                  list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
#                  list_type_p_=list_type_p,thr_p_=thr_p
#){
#  print("Starting gamm_fc().")
#  nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T)
#  dict_roi <- func_dict_roi(paths_)
#  
#  # Load and subset clinical data according to specified subsetting condition and covariate availability
#  print('Loading clinical data.')
#  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
#                                     list_covar=list_covar_,rem_na_clin=T)
#  df_clin<-data_clin$df_clin
#  
#  
#  for (atlas in list_atlas_){
#    
#    #****************************
#    # ROI-wise FC GAMM calculation
#    #****************************
#    # Load ROI-wise FC data
#    print(paste('Loading FC data, atlas:',atlas,sep=' '))
#    df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
#    df_join<-join_fc_clin(df_fc,df_clin)
#    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_src.csv",sep="")),
#              row.names=F)
#    
#    # Calculate and save ROI-wise GAMM of FC
#    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
#    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
#    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
#    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
#    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
#    write.csv(data_gamm$df_out_gamm,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm.csv",sep="")),row.names = F)
#    write.csv(data_gamm$df_out_aic,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_aic.csv",sep="")),row.names = F)
#    
#    # Calculate multiple comparison-corrected p values
#    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
#                             list_mod,list_plot,calc_seed_level=T)
#    write.csv(df_plot_gamm,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_plt.csv",sep="")),row.names = F)
#    
#    # Graphical output of ROI-wise GAMM of FC
#    plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
#                list_type_p_,thr_p,paths_)
#    
#    #****************************
#    # Group-wise FC GAMM calculation
#    #****************************
#    # Load group-wise FC data
#    print(paste('Loading group FC data, atlas:',atlas,sep=' '))
#    df_fc_grp<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc_grp.csv',sep='')))
#    df_join_grp<-join_fc_clin(df_fc_grp,df_clin)
#    write.csv(df_join_grp,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_src.csv",sep="")),
#              row.names=F)
#    
#    # Calculate and save group-wise GAMM of FC
#    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
#    list_roi_grp<-sort(unique(c(as.character(df_join_grp$from),as.character(df_join_grp$to))))
#    #df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
#    df_roi_grp<-data.frame(id=list_roi_grp,label=capitalize(list_roi_grp),group="group")
#    data_gamm_grp<-iterate_gamm(df_join_grp,df_roi_grp,list_mod_)
#    write.csv(data_gamm_grp$df_out_gamm,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm.csv",sep="")),row.names = F)
#    write.csv(data_gamm_grp$df_out_aic,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_aic.csv",sep="")),row.names = F)
#    
#    # Calculate multiple comparison-corrected p values
#    df_plot_gamm_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_roi_grp,analysis="grp",atlas,
#                                 list_mod,list_plot,calc_seed_level=T)
#    write.csv(df_plot_gamm_grp,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_plt.csv",sep="")),row.names = F)
#    
#    # Graphical output of group-wise GAMM of FC
#    plot_gam_fc(df_plot_gamm_grp,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
#                list_type_p_,thr_p,paths_)
#    
#    #****************************
#    # Multi-scale FC GAMM calculation
#    #****************************
#    # Subset ROI-wise GAMM result to include only within-group connections
#    df_gamm_ms<-NULL
#    for (group in list_roi_grp){
#      list_roi_within_grp<-as.character(df_roi[df_roi$group==group,"id"])
#      df_gamm_ms_add<-data_gamm$df_out_gamm[which(is.element(as.character(data_gamm$df_out_gamm[,"from"]),list_roi_within_grp)
#                                                  & is.element(as.character(data_gamm$df_out_gamm[,"to"]),list_roi_within_grp)),]
#      df_gamm_ms_add<-cbind(group=group,df_gamm_ms_add)
#      df_gamm_ms<-rbind(df_gamm_ms,df_gamm_ms_add)
#    }
#    
#    # Combine within-group ROI-wise GAMM results and between-group GAMM results
#    df_gamm_ms<-rbind(df_gamm_ms,cbind(group="group",data_gamm_grp$df_out_gamm))
#    
#    # Calculate multiple comparison-corrected p values
#    df_plot_gamm_ms<-add_mltcmp(df_gamm_ms,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
#                                calc_seed_level=F)
#    write.csv(df_plot_gamm_ms,
#              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-ms_gamm_plt.csv",sep="")),row.names = F)
#    
#    # Split data into ROI-wise and group-wise GAMM results, graphical output
#    for (group in list_roi_grp){
#      df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group==group,-1]
#      df_roi_split<-df_roi[df_roi$group==group,]
#      label_analysis<-paste("ms_grp-",group,sep="")
#      plot_gamm_fc(df_plot_gamm_ms_split,df_roi_split,analysis=label_analysis,atlas,list_mod,list_plot,
#                   list_type_p_,thr_p,paths_)
#    }
#    df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group=="group",-1]
#    plot_gam_fc(df_plot_gamm_ms_split,df_roi_grp,analysis="ms_grp-group",atlas,list_mod,list_plot,
#                list_type_p_,thr_p,paths_)
#    
#  }
#  print('Finished gamm_fc().')
#}

#**************************************************
# Component analyses of FC ========================
#**************************************************
# OBSOLETE
#ca_fc<-function(paths_=paths,list_atlas_=list_atlas,list_wave_=list_wave,
#                list_covar_=list_covar,subset_subj_=subset_subj,list_dim_ca_=list_dim_ca,
#                plot_result=F){
#  print("Starting ca_fc().")
#  memory.limit(200000)
#  nullobj<-func_createdirs(paths_,str_proc="ca_fc()",copy_log=T)
#  
#  # Load and subset clinical data according to specified subsetting condition and covariate availability
#  print('Loading clinical data.')
#  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=F)
#  df_clin<-data_clin$df_clin
#  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
#  
#  for (atlas in list_atlas_){
#    # Load and examine FC data
#    print(paste("Loding FC of atlas: ",atlas,sep=""))
#    
#    # Create graph edge dataframe and node list
#    #df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
#    df_conn<-fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
#    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
#    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
#    n_edge<-dim(df_edge)[1]
#    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
#    n_node<-length(list_node)
#    
#    # Create list of subjects who meet subsetting condition and whose MRI data exist
#    list_ses_exist <- sort(unique(df_conn$ses))
#    list_id_subj_exist<-list()
#    for (ses in list_ses_exist){
#      df_conn_ses<-df_conn[df_conn$ses==ses,]
#      id_subj_exist<-unique(df_conn_ses$ID_pnTTC)
#      id_subj_subset<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
#      id_subj_exist<-intersect(id_subj_exist,id_subj_subset)
#      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist)
#    }
#    
#    # Cbind FC data (Fisher-z transform of FC) as input for PCA function
#    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
#    df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
#    colnames(df_clin_exist)<-colnames(df_clin)
#    for (ses in list_ses_exist){
#      for (id_subj in list_id_subj_exist[[ses]]){
#        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
#        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
#        df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
#        df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ses==ses & df_clin$ID_pnTTC==id_subj,])
#      }
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
#    write.csv(data_pca$df_comp_mri,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variable.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_comp_subj,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_subject.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_variance,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variance.csv",sep="")),row.names=F)
#    write.csv(data_pca$df_comp_clin_flat,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_correlation.csv",sep="")),row.names=F)
#    
#    # Plot PCA results
#    if (plot_result){
#      print("Plotting PCA of FC.")
#      list_plot_pca<-plot_ca(df_src=data_pca$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_pca$dim)
#      for (i_dim in names(list_plot_pca)){
#        for (name_covar in names(list_plot_pca[[i_dim]])){
#          plot<-list_plot_pca[[i_dim]][[name_covar]]
#          plot<-(plot
#                 + ggtitle("PCA of FC"))
#          ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_pca.eps",sep=""),plot=plot,device=cairo_ps,
#                 path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
#        }
#      }
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
#    write.csv(df_comp_mri_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_variable.csv",sep="")),row.names=F)
#    write.csv(df_comp_subj_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_subject.csv",sep="")),row.names=F)
#    write.csv(df_variance_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_variance.csv",sep="")),row.names=F)
#    write.csv(df_comp_clin_flat_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_ica_correlation.csv",sep="")),row.names=F)
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