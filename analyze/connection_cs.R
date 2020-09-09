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
#path_exp_full<-NULL
path_exp_full<-"/media/atiroms/HDD_04/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"421_fc_aroma"
dir_out<-"423_fc_gam_aroma"

#dir_in<-"421_fc_aroma"
#dir_out<-"425_fc_ca_aroma"

#list_dim_ca<-c(10,20,40)
#list_dim_ca<-c(5,10,20,40)
#list_dim_ca<-c(5,10)
list_dim_ca<-10
ratio_vis<-0.01

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100x7","schaefer200x7","schaefer400x7",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
#list_atlas<-c("aal116","power264","shen268")
list_atlas<-"aal116"
#list_atlas<-"schaefer400x7"

list_type_p=c("p","p_bh","seed_p_bh")
thr_p <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(data.table)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# Iterate ca_fc() over clinical variables =========
# and waves =======================================
#**************************************************

comp_clin_cor<-function(df_comp_subj,df_clin,n_covar,list_sex){
  list_dim<-sort(unique(df_comp_subj$dim))
  df_cor_rbind<-df_cor_flat_rbind<-NULL
  for (dim in list_dim){
    for (label_sex in names(list_sex)){
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
      df_cor_flat<-df_cor_flat[df_cor_flat$from %in% colnames(df_clin) & df_cor_flat$to %in% colnames(df_comp_subj_subset),]
      df_cor_flat<-df_cor_flat[,c("from","to","r","p")]
      colnames(df_cor_flat)<-c("covar","comp","r","p")
      df_cor_flat$comp<-sapply(sapply(df_cor_flat$comp,substr,start=6,stop=8),as.integer)
      df_cor_flat<-cbind(dim=dim,df_cor_flat)
      df_cor_rbind<-rbind(df_cor_rbind,cbind(sex=label_sex,df_cor))
      df_cor_flat_rbind<-rbind(df_cor_flat_rbind,cbind(sex=label_sex,df_cor_flat))
    }
  }
  return(list("df_cor"=df_cor_rbind,"df_cor_flat"=df_cor_flat_rbind))
}

ca_fc_cs_multi<-function(paths_=paths,list_waves_=ca_fc_list_waves,subset_subj_=ca_fc_subset_subj,
                         list_sex_=ca_fc_list_sex,list_atlas_=list_atlas,
                         list_covar_tanner_=ca_fc_list_covar_tanner,list_tanner_=ca_fc_list_tanner,
                         list_covar_hormone_=ca_fc_list_covar_hormone,list_hormone_=ca_fc_list_hormone,
                         list_dim_ca_=list_dim_ca,ratio_vis_=ratio_vis){
  print("Starting ca_fc_cs_multi()")
  nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs_multi()",copy_log=T)
  df_cor<-df_ca_subj_bind<-df_ca_var_bind<-df_ca_vaf_bind<-NULL
  for (atlas in list_atlas_){
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
    for (waves in names(list_waves_)){
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
                      ", Sex: ",label_sex,", preaparing data.",sep=""))
          
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
          
          # Visualize component-FC matrix in circular plot
          plot_ca_fc(paths_=paths_,data_pca$df_comp_mri,atlas=atlas,dim_ca=dim_ca,
                     ratio_vis=ratio_vis,method="pca",label_sex=label_sex,ses=wave_mri)
          
          # Save results
          df_pca_mri<-rbind(df_pca_mri,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_mri))
          df_pca_subj<-rbind(df_pca_subj,cbind(sex=label_sex,dim=dim_ca,data_pca$df_comp_subj))
          df_pca_vaf<-rbind(df_pca_vaf,cbind(sex=label_sex,dim=dim_ca,data_pca$df_vaf))
          df_ca_var_bind<-rbind(df_ca_var_bind,cbind(atlas=atlas,method="pca",df_pca_mri))
          df_ca_subj_bind<-rbind(df_ca_subj_bind,cbind(atlas=atlas,method="pca",df_pca_subj))
          df_ca_vaf_bind<-rbind(df_ca_vaf_bind,cbind(atlas=atlas,method="pca",df_pca_vaf))
        
          data_pca<-NULL
          gc()
          
          # Calculate ICA of FC
          for (dim_ca in list_dim_ca_){
            data_ica<-func_ica(df_src=df_conn_calc,df_var=df_edge,df_indiv=df_clin_exist,dim_ca=dim_ca,calc_corr=F)
            
            # Visualize component-FC matrix in circular plot
            plot_ca_fc(paths_=paths_,data_ica$df_comp_mri,atlas=atlas,dim_ca=dim_ca,
                       ratio_vis=ratio_vis,method="ica",label_sex=label_sex,ses=wave_mri)
            
            # Save results
            df_ica_mri<-rbind.fill(df_ica_mri,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_mri))
            df_ica_subj<-rbind.fill(df_ica_subj,cbind(sex=label_sex,dim=dim_ca,data_ica$df_comp_subj))
            df_ica_vaf<-rbind.fill(df_ica_vaf,cbind(sex=label_sex,dim=dim_ca,data_ica$df_vaf))
          }
          df_ca_var_bind<-rbind(df_ca_var_bind,cbind(atlas=atlas,method="ica",df_ica_mri))
          df_ca_subj_bind<-rbind(df_ca_subj_bind,cbind(atlas=atlas,method="ica",df_ica_subj))
          df_ca_vaf_bind<-rbind.fill(df_ca_vaf_bind,cbind(atlas=atlas,method="ica",df_ica_vaf))
  
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
      
      # Calculate component-clinical correlation
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
                                n_covar=n_covar,list_sex=list_sex_)
        df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_tanner,
                                   method="pca",data_cor$df_cor_flat))
        data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                                n_covar=n_covar,list_sex=list_sex_)
        df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_tanner,
                                   method="ica",data_cor$df_cor_flat))
      } # finished looping over Tanner stages
      
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
                                n_covar=n_covar,list_sex=list_sex_)
        df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_hormone,
                                   method="pca",data_cor$df_cor_flat))
        data_cor<-comp_clin_cor(df_comp_subj=df_ica_subj,df_clin=df_clin,
                                n_covar=n_covar,list_sex=list_sex_)
        df_cor<-rbind(df_cor,cbind(atlas=atlas,ses=waves,variable=idx_hormone,
                                   method="ica",data_cor$df_cor_flat))
      } # finished looping over Hormones
    } # finished looping over waves
  } # finished looping over atlases
  
  # Add ROI label to component-MRI variable matrix
  df_roi<-func_dict_roi(paths=paths_)
  df_roi<-df_roi[,c("id","label")]
  df_ca_var_bind<-left_join(df_ca_var_bind,df_roi,by=c("from"="id"))
  df_ca_var_bind<-rename(df_ca_var_bind,c("label"="from_label"))
  df_ca_var_bind<-left_join(df_ca_var_bind,df_roi,by=c("to"="id"))
  df_ca_var_bind<-rename(df_ca_var_bind,c("label"="to_label"))
  
  # Save results
  write.csv(df_ca_var_bind,file.path(paths_$output,"output",paste("fc_ca_var.csv",sep="")),row.names=F)
  write.csv(df_ca_subj_bind,file.path(paths_$output,"output",paste("fc_ca_subj.csv",sep="")),row.names=F)
  write.csv(df_ca_vaf_bind,file.path(paths_$output,"output",paste("fc_ca_vaf.csv",sep="")),row.names=F)
  write.csv(df_cor,file.path(paths_$output,"output",paste("fc_ca_cor.csv",sep="")),row.names=F)
  
  print("Finished ca_fc_cs_multi()")
}


#**************************************************
# Component analyses of FC in cross section =======
#**************************************************
# OBSOLETE
ca_fc_cs<-function(paths_=paths,list_atlas_=list_atlas,wave_clin_=wave_clin,wave_mri_=wave_mri,
                   list_covar_=list_covar,subset_subj_=subset_subj,list_dim_ca_=list_dim_ca,
                   plot_result=F,suffix_=""){
  
  print("Starting ca_fc().")
  #memory.limit(200000)
  nullobj<-func_createdirs(paths_,str_proc="ca_fc()",copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,wave_clin_,subset_subj_,list_covar_,
                                     rem_na_clin=F,suffix=suffix_)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  for (atlas in list_atlas_){
    # Load and examine FC data
    print(paste("Loading FC of atlas: ",atlas,sep=""))
    
    # Create graph edge dataframe and node list
    #df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_conn<-fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_conn<-df_conn[df_conn$ses==wave_mri_,]
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
    n_edge<-dim(df_edge)[1]
    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
    n_node<-length(list_node)
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    list_subj_mri<-unique(df_conn[df_conn$ses==wave_mri_,]$ID_pnTTC)
    list_subj_clin<-unique(df_clin[df_clin$ses==wave_clin_,]$ID_pnTTC)
    list_subj_calc<-intersect(list_subj_mri,list_subj_clin)
    print(paste("MRI data absent in",
                as.character(length(list_subj_clin)-length(list_subj_calc)),
                "subjects.",sep=" "))
    
    # Cbind FC data (Fisher-z transform of FC) as input for PCA function
    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
    df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
    colnames(df_clin_exist)<-colnames(df_clin)
    for (id_subj in list_subj_calc){
      df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
      df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
      df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ID_pnTTC==id_subj,])
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
    write.csv(data_pca$df_comp_mri,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_variable.csv",sep="")),row.names=F)
    write.csv(data_pca$df_comp_subj,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_subject.csv",sep="")),row.names=F)
    write.csv(data_pca$df_vaf,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_variance.csv",sep="")),row.names=F)
    write.csv(data_pca$df_comp_clin_flat,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_pca_correlation.csv",sep="")),row.names=F)
    
    # Plot PCA results
    if (plot_result){
      print("Plotting PCA of FC.")
      #list_plot_pca<-plot_ca(df_src=data_pca$df_comp_subj,list_name_covar=names(list_covar_),n_dim=data_pca$dim)
      #for (i_dim in names(list_plot_pca)){
      #  for (name_covar in names(list_plot_pca[[i_dim]])){
      #    plot<-list_plot_pca[[i_dim]][[name_covar]]
      #    plot<-(plot
      #           + ggtitle("PCA of FC"))
      #    ggsave(paste("atl-",atlas,"_comp-",sprintf("%03d",as.numeric(i_dim)),"-",sprintf("%03d",as.numeric(i_dim)+1),"_cov-",name_covar,"_pca.eps",sep=""),plot=plot,device=cairo_ps,
      #           path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
      #  }
      #}
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
    write.csv(df_comp_mri_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_variable.csv",sep="")),row.names=F)
    write.csv(df_comp_subj_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_subject.csv",sep="")),row.names=F)
    write.csv(df_variance_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_variance.csv",sep="")),row.names=F)
    write.csv(df_comp_clin_flat_rbind,file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_ica_correlation.csv",sep="")),row.names=F)
    
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
# Iterate gam_fc_cs() over clinical variables =====
# and waves =======================================
#**************************************************

gam_fc_cs_multi<-function(paths_=paths,list_atlas_=list_atlas,
                          list_waves_=gam_fc_list_waves,subset_subj_=gam_fc_subset_subj,
                          list_covar_tanner_=gam_fc_list_covar_tanner,list_tanner_=gam_fc_list_tanner,
                          list_mod_tanner_=gam_fc_list_mod_tanner,list_plot_tanner_=gam_fc_list_plot_tanner,
                          list_covar_hormone_=gam_fc_list_covar_hormone,list_hormone_=gam_fc_list_hormone,
                          list_mod_hormone_=gam_fc_list_mod_hormone,list_plot_hormone_=gam_fc_list_plot_hormone,
                          list_type_p_=list_type_p,thr_p_=thr_p){
  print("Starting gam_fc_cs_multi()")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)
  
  for (waves in names(list_waves_)){
    wave_clin<-list_waves_[[waves]]$wave_clin
    wave_mri<-list_waves_[[waves]]$wave_mri
    print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))
    
    # Prepare subject subsetting condition (MRI QC criteria) according to specified waves
    subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
    subset_subj_temp<-list(subset_subj_temp)
    names(subset_subj_temp)<-wave_clin

    #1 Tanner stage
    for (idx_tanner in names(list_tanner_)){
      print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
      list_covar<-list_covar_tanner_
      list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
      suffix<-paste("ses-",waves,"_var-",idx_tanner,sep="")
      
      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
                         wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas_,
                         list_mod_=list_mod_tanner_,list_plot_=list_plot_tanner_,
                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
                         suffix_=suffix)
    } # finished looping over Tanner stages
    
    
    #2 Hormones
    for (idx_hormone in names(list_hormone_)){
      print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
      list_covar<-list_covar_hormone_
      list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
      suffix<-paste("ses-",waves,"_var-",idx_hormone,sep="")
      
      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
                         wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                         list_mod_=list_mod_hormone_,list_plot_=list_plot_hormone_,
                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
                         suffix_=suffix)
    } # finished looping over Hormones
    
  } # finished looping over waves
  print("Finished gam_fc_cs_multi()")
}



#**************************************************
# Additive/Linear model of FC in cross-section ====
#**************************************************

join_fc_clin<-function(df_fc,df_clin,wave_clin,wave_mri){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  df_fc<-df_fc[df_fc$ses==wave_mri,]
  df_fc$ses<-NULL
  #colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c("ID_pnTTC","from","to","value")]
  
  df_clin<-df_clin[df_clin$wave==wave_clin,]
  df_clin$wave<-NULL
  
  # Join clinical and FC data frames
  print('Joining clinical and FC data.')
  #df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  df_join<-inner_join(df_fc,df_clin,by='ID_pnTTC')
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  return(df_join)
}

gam_fc_cs<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                    wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                    list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
                    list_type_p_=list_type_p,thr_p_=thr_p,suffix_=suffix
                    ){
  print("Starting gam_fc_cs().")
  #nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T,prefix=suffix_,print_terminal=F)
  df_clin<-data_clin$df_clin
  
  for (atlas in list_atlas_){
    
    # Load ROI-wise FC data
    print(paste('Loading FC data, atlas:',atlas,sep=' '))
    #df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_fc<-as.data.frame(fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep=''))))
    df_join<-join_fc_clin(df_fc,df_clin,wave_clin_,wave_mri_)
    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save ROI-wise GAMM of FC
    print(paste('Calculating GAM, atlas: ',atlas,sep=''))
    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
    write.csv(data_gamm$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod_,list_plot_,calc_seed_level=T)
    write.csv(df_plot_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam_plt.csv",sep="")),row.names = F)
    
    # Graphical output of ROI-wise GAMM of FC
    plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod_,list_plot_,
                list_type_p_,thr_p_,paths_,suffix_)
    
  }
  print('Finished gam_fc_cs().')
}
