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
#path_exp_full<-"/media/atiroms/SSD_02/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"421_fc_aroma"
#dir_out<-"423.1_fc_gam_aroma_test5"
dir_out<-"424_fc_gamm_aroma_test5"
#list_atlas<-c("aal116","gordon333","ho112","power264",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
#list_atlas<-c("aal116")
list_atlas<-c("ho112")
#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100x7","schaefer200x7","schaefer400x7",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")


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
gam_fc_diff_core<-function(paths,data_fc,atlas,param,list_sex,
                           list_covar,list_mod,list_term,idx_var,
                           calc_parallel,test_mod
                           ){
  # Prepare clinical data and demean
  data_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                     prefix=paste("var-",idx_var,sep=""),print_terminal=F)
  list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
  df_clin_diff<-data_clin$df_clin
  colnames(df_clin_diff)[colnames(df_clin_diff)=="wave"]<-"ses"
  df_clin_diff<-func_clinical_data_diffmean(df_clin_diff,list_id_subj,list_covar)
  df_clin_diff<-func_demean_clin(df_clin_diff,thr_cont=6,separate_sex=T)$df_clin # thr_cont=3 to demean Tanner, =6 not to
  df_clin_diff$wave<-"2-1"
  fwrite(df_clin_diff,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_clin_diff.csv",sep="")),row.names=F)
  
  # Prepare FC data
  df_fc_diff<-data_fc$df_fc[data_fc$df_fc$ses=="2-1",]
  df_fc_grp_diff<-data_fc$df_fc_grp[data_fc$df_fc_grp$ses=="2-1",]
  fwrite(df_fc_diff,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc_diff.csv",sep="")),row.names=F)
  fwrite(df_fc_grp_diff,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc_grp_diff.csv",sep="")),row.names=F)
  
  # Calculate model
  file_check<-file.path(paths$output,"output","temp",
                        paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep=""))
  if (file.exists(file_check)){
    print("Calculated result already exists.")
    df_gamm<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))
    df_anova<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))
    df_gamm_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))
    df_anova_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))
  }else{
    # Join FC and clinical data
    df_join_diff<-join_fc_clin(df_fc_diff,df_clin_diff)
    df_join_grp_diff<-join_fc_clin(df_fc_grp_diff,df_clin_diff)
    
    # Prepare parallelization cluster
    if (calc_parallel){
      clust<-makeCluster(floor(detectCores()*3/4))
    }else{
      clust<-makeCluster(1)
    }
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                                  "as.formula","as.numeric.factor",
                                  "lm","lmer","gam",
                                  "summary","anova","summary.gam","anova.gam","AIC"),
                  envir=environment())
    
    # Calculate model
    data_gamm<-iterate_gamm4(clust,df_join_diff,data_fc$df_edge,progressbar=F,test_mod=test_mod)
    data_gamm_grp<-iterate_gamm4(clust,df_join_grp_diff,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
    stopCluster(clust)
    
    # Add multiple comparison-corrected p values
    df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    
    # Save results
    fwrite(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)
    fwrite(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
    fwrite(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
    fwrite(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
    fwrite(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep="")),row.names = F)
    fwrite(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
  } # End if file exists
  
  # Threshold and plot graph edges
  if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")))){
    df_plot<-read.csv(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")))
  }else{
    df_plot<-df_plot_grp<-data.frame()
    for (idx_mod in names(list_mod)){
      for (idx_term in names(list_term)){
        var_exp<-list_term[[idx_term]][["var_exp"]]
        for (idx_sex in list_sex){
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
            if (idx_sex==1){label_sex<-"m"}else{label_sex<-"f"}
            plot_gamm<-annotate_figure(plot_gamm,
                                       top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                                             ", expvar: ",var_exp,", sex: ",label_sex,", p value: all",sep=""),
                                                       color = "black", size = 14))
            ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,
                         "_sex-",label_sex,"_pval-all_net.png",sep=""),
                   plot=plot_gamm,path=file.path(paths$output,"output","plot"),height=13,width=10,dpi=600)
            for (p in param$list_p){
              df_plot_subset<-df_gamm_subset[df_gamm_subset[[p$type]]<p$threshold,]
              df_plot_grp_subset<-df_gamm_grp_subset[df_gamm_grp_subset[[p$type]]<p$threshold,]
              plot_gamm<-plot_gam_fc3(df_plot_subset,df_plot_grp_subset,data_fc)
              plot_gamm<-annotate_figure(plot_gamm,
                                         top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                                               ", expvar: ",var_exp,", sex: ",label_sex,", p value: ",p$type,"<",p$threshold,sep=""),
                                                         color = "black", size = 14))
              ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,
                           "_sex-",label_sex,"_pval-",p$type,"_",p$threshold,
                           "_net.png",sep=""),
                     plot=plot_gamm,path=file.path(paths$output,"output","plot"),height=13,width=10,dpi=600)
              df_head<-data.frame(p_type=p$type,p_threshold=p$threshold)
              if (nrow(df_plot_subset)>0){
                df_plot<-bind_rows(df_plot,cbind(df_head,df_plot_subset))
              }
              if (nrow(df_plot_grp_subset)>0){
                df_plot_grp<-bind_rows(df_plot_grp,cbind(df_head,df_plot_grp_subset))
              }
            }# end of loop over list_p
          }
        }
      }
    }
    # Save results
    if (nrow(df_plot)>0){
      fwrite(df_plot,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")),row.names = F)
    }
    if (nrow(df_plot_grp)>0){
      fwrite(df_plot_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot_grp.csv",sep="")),row.names = F)
    }
  }
  
  # Detect sub-network by breadth-first approach
  if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")))){
    df_net<-df_node<-df_size_net<-NULL
    #list_output<-list()
    for (idx_mod in param$param_nbs$list_mod){
      for (idx_term in param$param_nbs$list_term){
        var_exp<-list_term[[idx_term]][["var_exp"]]
        for (idx_sex in list_sex){
          df_sign<-df_plot[df_plot$p_type=="p" & df_plot$p_threshold==param$param_nbs$p_cdt_threshold
                           & df_plot$model==idx_mod & df_plot$term==var_exp & df_plot$sex==idx_sex,]
          data_bfs<-func_bfs(df_sign)
          if (idx_sex==1){label_sex<-"m"}else{label_sex<-"f"}
          if(length(data_bfs$list_network)>0){
            for (idx_net in seq(length(data_bfs$list_network))){
              network<-data_bfs$list_network[[idx_net]]
              plot_subnet<-plot_circular2(df_edge=network$df_edge,df_node=network$df_node,df_roi=data_fc$df_roi,rule_order="degree")
              plot_subnet<-(plot_subnet
                            +ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                           ", expvar: ",var_exp,", sex: ",label_sex,", p value: p<",param$param_nbs$p_cdt_threshold,
                                           ", #",as.character(idx_net),sep="")))
              ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,
                           "_sex-",label_sex,"_pval-p_",param$param_nbs$p_cdt_threshold,
                           "_idx-",as.character(idx_net),"_subnet.png",sep=""),
                     plot=plot_subnet,path=file.path(paths$output,"output","plot"),height=10,width=10,dpi=600)
              df_head<-data.frame(model=idx_mod,term=var_exp,sex=idx_sex)
              df_net<-rbind(df_net,data.frame(id_net=idx_net,network$df_edge))
              df_node_add<-inner_join(data.frame(id_net=idx_net,network$df_node),data_fc$df_roi,by=c("node"="id"))
              df_node_add<-dplyr::rename(df_node_add,"label_node"="label","group_node"="group")
              df_node<-rbind(df_node,cbind(df_head,df_node_add))
              df_size_net<-rbind(df_size_net,cbind(df_head,data.frame(id_net=idx_net,size=network$size_net)))
            }
          }
        }
      }
    }
    fwrite(df_net,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_edge.csv",sep="")),row.names = F)
    fwrite(df_node,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_node.csv",sep="")),row.names = F)
    fwrite(df_size_net,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_size.csv",sep="")),row.names = F)
  }
  
  # Permutation for NBS
  # Prepare parallelization cluster
  if (calc_parallel){
    clust<-makeCluster(floor(detectCores()*3/4))
  }else{
    clust<-makeCluster(1)
  }
  clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                                "as.formula","as.numeric.factor",
                                "lm","lmer","gam",
                                "summary","anova","summary.gam","anova.gam","AIC"),
                envir=environment())
  set.seed(0)
  pb<-txtProgressBar(min=0,max=param$param_nbs$n_perm,style=3,width=50)
  df_max_size<-data.frame()
  for (idx_perm in seq(param$param_nbs$n_perm)){
    for (idx_term in param$param_nbs$list_term){
      var_exp<-list_term[[idx_term]][["var_exp"]]
      # Sex-wise permutation of term (expvar) of interst
      df_clin_perm<-NULL
      for (idx_sex in list_sex){
        df_clin_perm_add<-df_clin_diff[df_clin_diff$sex==idx_sex,]
        df_clin_perm_add[,var_exp]<-sample(df_clin_perm_add[,var_exp])
        df_clin_perm<-rbind(df_clin_perm,df_clin_perm_add)
      }
      # Join FC and permuted clinical data
      df_join_diff<-join_fc_clin(df_fc_diff,df_clin_perm)
      
      # Calculate model
      data_gamm<-iterate_gamm4(clust,df_join_diff,data_fc$df_edge,progressbar=F,test_mod=test_mod)
      df_gamm<-data_gamm$df_gamm
      df_anova<-data_gamm$df_anova
      for (idx_mod in param$param_nbs$list_mod){
        for (idx_sex in list_sex){
          # Subset GAMM result dataframe for plotting
          df_gamm_subset<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp & df_gamm$sex==idx_sex,]
          if (nrow(df_gamm_subset)==0){
            df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp & df_anova$sex==idx_sex,]
            if (nrow(df_anova_subset)>0){
              # In case the term does not exist in df_gamm, plot using df_anova instead
              df_gamm_subset<-df_anova_subset
            }
          }
          if (nrow(df_gamm_subset)>0){ # If the model/expvar/sex exist either in df_gamm or df_anova
            df_sign<-df_gamm_subset[df_gamm_subset$p<param$param_nbs$p_cdt_threshold,]
            max_size<-func_bfs(df_sign)$max_size
            df_max_size<-rbind(df_max_size,data.frame(id_perm=idx_perm,model=idx_mod,term=var_exp,sex=idx_sex,max_size=max_size))
          }
        }
      }
    }
    setTxtProgressBar(pb,idx_perm)
  }
  stopCluster(clust)
  close(pb)
  fwrite(df_max_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_perm_max.csv",sep="")),row.names = F)
  
  # Summarize permutation result
  df_threshold_size<-NULL
  for (idx_mod in param$param_nbs$list_mod){
    for (idx_term in param$param_nbs$list_term){
      var_exp<-list_term[[idx_term]][["var_exp"]]
      for (idx_sex in list_sex){
        list_max_size<-df_max_size[df_max_size$model==idx_mod & df_max_size$term==var_exp
                                   & df_max_size$sex==idx_sex,"max_size"]
        list_max_size<-sort(list_max_size)
        thr_size_nbs<-list_max_size[ceiling(length(list_max_size)*(1-param$param_nbs$p_perm_threshold))]
        df_threshold_size<-rbind(df_threshold_size,data.frame(model=idx_mod,term=var_exp,sex=idx_sex,
                                                              thr_size=thr_size_nbs))
        if (idx_sex==1){
          label_sex<-"m";title_sex<-"male";color_plt<-"steelblue2"
        }else{
          label_sex<-"f";title_sex<-"female";color_plt<-"lightcoral"
        }
        title_plot<-list_term[[idx_term]][["title"]]
        plot_permutation(paths,list_max=list_max_size,thr_size_nbs,
                         atlas,var=idx_var,wave="2-1",idx_mod,idx_term,label_sex,title_plot,title_sex,color_plt)
      }
    }
  }
  fwrite(df_threshold_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_perm_thr.csv",sep="")),row.names = F)
}

gam_fc_diff<-function(paths_=paths,list_atlas_=list_atlas,param=param_gam_fc_diff){
  print("Starting gam_fc_diff().")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc_diff()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    print(paste("Preparing FC data: ",atlas,sep=""))
    data_fc<-prep_data_fc(paths_,atlas,param$key_group,include_diff=T,abs_nfc=param$abs_nfc)
    data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
    data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_tanner in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gam_fc_diff_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
                       param$list_mod_tanner,param$list_term_tanner,idx_tanner,
                       calc_parallel=T,test_mod=F)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(param$list_hormone)){
      print(paste("Atlas: ",atlas,", Hormone type: ",param$list_hormone[[idx_hormone]][["label"]],sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
      gam_fc_diff_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
                       param$list_mod_hormone,param$list_term_hormone,idx_hormone,
                       calc_parallel=T,test_mod=F)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  df_gamm<-df_aic<-df_plot<-df_anova<-df_gamm_grp<-df_aic_grp<-df_plot_grp<-df_anova_grp<-df_net<-df_node<-df_size_net<-df_perm_max<-df_perm_thr<-data.frame()
  list_var<-c(param$list_tanner,param$list_hormone)
  for (atlas in list_atlas_){
    for (idx_var in names(list_var)){
      df_head<-data.frame(atlas=atlas,variable=idx_var)
      # ROI GAM results
      df_gamm<-rbind(df_gamm,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))))
      if(file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot.csv",sep="")))){
        df_plot<-bind_rows(df_plot,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot.csv",sep=""))))))
      }
      df_anova<-rbind(df_anova,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))))
      df_aic<-rbind(df_aic,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep=""))))))
      # Group GAM results
      df_gamm_grp<-rbind(df_gamm_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))))
      if(file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot_grp.csv",sep="")))){
        df_plot_grp<-bind_rows(df_plot_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot_grp.csv",sep=""))))))
      }
      df_anova_grp<-rbind(df_anova_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))))
      df_aic_grp<-rbind(df_aic_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep=""))))))
      # Network-based statistics results
      df_net<-rbind(df_net,cbind(df_head,as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_edge.csv",sep=""))))))
      df_node<-rbind(df_node,cbind(df_head,as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_node.csv",sep=""))))))
      df_size_net<-rbind(df_size_net,cbind(df_head,as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_size.csv",sep=""))))))
      # Permutation results
      df_perm_max<-rbind(df_perm_max,cbind(df_head,as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_perm_max.csv",sep=""))))))
      df_perm_thr<-rbind(df_perm_thr,cbind(df_head,as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_perm_thr.csv",sep=""))))))
    }
  }
  fwrite(df_gamm,file.path(paths_$output,"output","result","gamm.csv"),row.names = F)
  if (nrow(df_plot)>0){
    fwrite(df_plot,file.path(paths_$output,"output","result","plot.csv"),row.names = F)
  }
  fwrite(df_anova,file.path(paths_$output,"output","result","gamm_anova.csv"),row.names = F)
  fwrite(df_aic,file.path(paths_$output,"output","result","gamm_aic.csv"),row.names = F)
  fwrite(df_gamm_grp,file.path(paths_$output,"output","result","gamm_grp.csv"),row.names = F)
  if (nrow(df_plot_grp)>0){
    fwrite(df_plot_grp,file.path(paths_$output,"output","result","plot_grp.csv"),row.names = F)
  }
  fwrite(df_anova_grp,file.path(paths_$output,"output","result","gamm_anova_grp.csv"),row.names = F)
  fwrite(df_aic_grp,file.path(paths_$output,"output","result","gamm_aic_grp.csv"),row.names = F)
  fwrite(df_net,file.path(paths_$output,"output","result","bfs_edge.csv"),row.names = F)
  fwrite(df_node,file.path(paths_$output,"output","result","bfs_node.csv"),row.names = F)
  fwrite(df_size_net,file.path(paths_$output,"output","result","bfs_size.csv"),row.names = F)
  fwrite(df_perm_max,file.path(paths_$output,"output","result","perm_max.csv"),row.names = F)
  fwrite(df_perm_thr,file.path(paths_$output,"output","result","perm_thr.csv"),row.names = F)
  print("Finished gam_fc_diff().")
}

#**************************************************
# GLMM/GAMM of FC =================================
#**************************************************
gamm_fc_core<-function(paths,data_fc,atlas,param,list_sex,
                       list_covar,list_mod,list_term,idx_var,
                       calc_parallel,test_mod
                       ){
  file_check<-file.path(paths$output,"output","temp",
                        paste("atl-",atlas,"_var-",idx_var,"_gamm_grp_aic.csv",sep=""))
  if (file.exists(file_check)){
    print("Calculated result already exists.")
    df_gamm<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))
    df_gamm_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))
    df_anova<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))
    df_anova_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))
  }else{
    # Prepare clinical data and demean
    df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                     prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
    df_clin<-func_demean_clin(df_clin,thr_cont=6,separate_sex=T)$df_clin # thr_cont=4 to demean Tanner, =5 not to
    df_join<-join_fc_clin(data_fc$df_fc,df_clin)
    df_join_grp<-join_fc_clin(data_fc$df_fc_grp,df_clin)
    
    # Calculate model
    if (calc_parallel){
      clust<-makeCluster(floor(detectCores()*3/4))
    }else{
      clust<-makeCluster(1)
    }
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                                  "as.formula","as.numeric.factor",
                                  "lm","lmer","gam",
                                  "summary","anova","summary.gam","anova.gam","AIC"),
                  envir=environment())
    data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
    data_gamm_grp<-iterate_gamm4(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
    stopCluster(clust)
    
    # Add multiple comparison-corrected p values
    df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
    df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
    
    # Save results
    write.csv(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
    write.csv(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)
    write.csv(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
    write.csv(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
    write.csv(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep="")),row.names = F)
  }
  
  df_plot<-df_plot_grp<-data.frame()
  for (idx_mod in names(list_mod)){
    for (idx_term in names(list_term)){
      var_exp<-list_term[[idx_term]][["var_exp"]]
      for (idx_sex in list_sex){
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
          if (idx_sex==1){label_sex<-"m"}else{label_sex<-"f"}
          plot_gamm<-annotate_figure(plot_gamm,
                                     top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                                           ", expvar: ",var_exp,", sex: ",label_sex,", p value: all",sep=""),
                                                     color = "black", size = 14))
          ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,
                       "_sex-",label_sex,"_pval-all_net.png",sep=""),
                 plot=plot_gamm,path=file.path(paths$output,"output","plot"),height=13,width=10,dpi=600)
          for (p in param$list_p){
            df_plot_subset<-df_gamm_subset[df_gamm_subset[[p$type]]<p$threshold,]
            df_plot_grp_subset<-df_gamm_grp_subset[df_gamm_grp_subset[[p$type]]<p$threshold,]
            df_plot_subset$F<-as.numeric(df_plot_subset$F)
            df_plot_grp_subset$F<-as.numeric(df_plot_grp_subset$F)
            plot_gamm<-plot_gam_fc3(df_plot_subset,df_plot_grp_subset,data_fc)
            plot_gamm<-annotate_figure(plot_gamm,
                                       top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                                             ", expvar: ",var_exp,", sex: ",label_sex,", p value: ",p$type,"<",p$threshold,sep=""),
                                                       color = "black", size = 14))
            ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,
                         "_sex-",label_sex,"_pval-",p$type,"_",p$threshold,
                         "_net.png",sep=""),
                   plot=plot_gamm,path=file.path(paths$output,"output","plot"),height=13,width=10,dpi=600)
            df_head<-data.frame(p_type=p$type,p_threshold=p$threshold)
            if (nrow(df_plot_subset)>0){
              df_plot<-bind_rows(df_plot,cbind(df_head,df_plot_subset))
            }
            if (nrow(df_plot_grp_subset)>0){
              df_plot_grp<-bind_rows(df_plot_grp,cbind(df_head,df_plot_grp_subset))
            }
          }# end of loop over list_p
        }
      }
    }
  }
  # Save results
  if (nrow(df_plot)>0){
    write.csv(df_plot,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")),row.names = F)
  }
  if (nrow(df_plot_grp)>0){
    write.csv(df_plot_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot_grp.csv",sep="")),row.names = F)
  }
}

gamm_fc<-function(paths_=paths,list_atlas_=list_atlas,param=param_gamm_fc){
  print("Starting gamm_fc().")
  nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T,list_param=param)
  memory.limit(1000000)
  
  # Loop over atlases
  for (atlas in list_atlas_){
    print(paste("Preparing FC data: ",atlas,sep=""))
    data_fc<-prep_data_fc(paths_,atlas,param$key_group,abs_nfc=param$abs_nfc)
    data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
    data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))
    
    # Loop over clinical variables
    #1 Tanner stage
    for (idx_tanner in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gamm_fc_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
                   param$list_mod_tanner,param$list_term_tanner,idx_tanner,
                   calc_parallel=T,test_mod=F)
    } # Finished looping over Tanner stages
    
    #2 Hormones
    for (idx_hormone in names(param$list_hormone)){
      print(paste("Atlas: ",atlas,", Hormone type: ",param$list_hormone[[idx_hormone]][["label"]],sep=""))
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
      gamm_fc_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
                   param$list_mod_hormone,param$list_term_hormone,idx_hormone,
                   calc_parallel=T,test_mod=F)
    } # Finished looping over Hormones
  } # Finished looping over atlas
  
  print("Combining results.")
  df_gamm<-df_aic<-df_plot<-df_anova<-df_gamm_grp<-df_aic_grp<-df_plot_grp<-df_anova_grp<-data.frame()
  list_var<-c(param$list_tanner,param$list_hormone)
  for (atlas in list_atlas_){
    for (idx_var in names(list_var)){
      df_head<-data.frame(atlas=atlas,variable=idx_var)
      df_gamm<-rbind(df_gamm,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))))
      if(file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot.csv",sep="")))){
        df_plot<-rbind(df_plot,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot.csv",sep=""))))))
      }
      df_anova<-rbind(df_anova,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))))
      df_aic<-rbind(df_aic,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep=""))))))
      df_gamm_grp<-rbind(df_gamm_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))))
      if(file.exists(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot_grp.csv",sep="")))){
        df_plot_grp<-rbind(df_plot_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var, "_plot_grp.csv",sep=""))))))
      }
      df_anova_grp<-rbind(df_anova_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))))
      df_aic_grp<-rbind(df_aic_grp,cbind(df_head,as.data.frame(fread(file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep=""))))))
    }
  }
  write.csv(df_gamm,file.path(paths_$output,"output","result","gamm.csv"),row.names = F)
  write.csv(df_plot,file.path(paths_$output,"output","result","plot.csv"),row.names = F)
  write.csv(df_anova,file.path(paths_$output,"output","result","gamm_anova.csv"),row.names = F)
  write.csv(df_aic,file.path(paths_$output,"output","result","gamm_aic.csv"),row.names = F)
  write.csv(df_gamm_grp,file.path(paths_$output,"output","result","gamm_grp.csv"),row.names = F)
  write.csv(df_plot_grp,file.path(paths_$output,"output","result","plot_grp.csv"),row.names = F)
  write.csv(df_anova_grp,file.path(paths_$output,"output","result","gamm_anova_grp.csv"),row.names = F)
  write.csv(df_aic_grp,file.path(paths_$output,"output","result","gamm_aic_grp.csv"),row.names = F)
  
  print("Finished gamm_fc().")
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
