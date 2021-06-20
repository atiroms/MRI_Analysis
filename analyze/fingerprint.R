#**************************************************
# Description =====================================
#**************************************************
# R script to analyze fingerprint data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for gta_bin() and gta_weight()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/str_FS"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
#path_exp_full <-"/media/atiroms/HDD_05/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full <-NULL

#dir_in<-"422_fp_aroma"
dir_in<-"422.1_fp_aroma_test2"
dir_out<-"429_qc_fp_aroma_test2"
#dir_out<-"426.1_fp_gam_aroma_test1"
#dir_out<-"428.1_fp_var_aroma_test2"
#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100x7","schaefer200x7","schaefer400x7",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
list_atlas<-c("aal116","gordon333","ho112","power264",
              "schaefer100x17","schaefer200x17","schaefer400x17",
              "shen268")
#list_atlas<-"ho112"
#list_atlas<-"aal116"
#list_atlas<-c("ho112","power264")


#**************************************************
# Libraries =======================================
#**************************************************
library(easypackages)
libraries("dplyr","mgcv","multcomp","parallel","DescTools","ggplot2","data.table")


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# GLM of Fingerprint change =======================
#**************************************************
qc_fp<-function(paths_=paths,list_atlas_=list_atlas){
  print("Starting qc_fp().")
  nullobj<-func_createdirs(paths_,str_proc="qc_fp()",copy_log=T)
  df_concat<-data.frame()
  for (atlas in list_atlas_){
    df_mean_fp<-as.data.frame(fread(file.path(paths_$input,"output","result",paste("atl-",atlas,"_mean_fp.csv",sep=""))))
    df_mean_fp<-df_mean_fp[df_mean_fp$group_1=="whole" & df_mean_fp$group_2=="whole",c("ses","ID_pnTTC","mean_z_r","qc")]
    mean_meanfp<-mean(df_mean_fp$mean_z_r)
    sd_meanfp<-sd(df_mean_fp$mean_z_r)
    df_mean_fp$z<-(df_mean_fp$mean_z_r-mean_meanfp)/sd_meanfp
    df_concat<-rbind(df_concat,data.frame(atlas=atlas,df_mean_fp))
    hist(df_mean_fp$z,breaks=30)
    hist(df_mean_fp[df_mean_fp$qc==1,"z"],breaks=30)
    
    df_out<-data.frame("ID_pnTTC"=seq(max(df_mean_fp$ID_pnTTC)),"W1_ZFP"=NA,"W2_ZFP"=NA)
    for (idx_row in seq(nrow(df_mean_fp))){
      ses<-df_mean_fp[idx_row,"ses"]
      id_subj<-df_mean_fp[idx_row,"ID_pnTTC"]
      zfp<-df_mean_fp[idx_row,"z"]
      if (ses==1){
        df_out[df_out$ID_pnTTC==id_subj,"W1_ZFP"]<-zfp
      }else{
        df_out[df_out$ID_pnTTC==id_subj,"W2_ZFP"]<-zfp
      }
    }
    df_out[is.na(df_out$W1_ZFP),"W1_ZFP"]<-df_out[is.na(df_out$W2_ZFP),"W2_ZFP"]<-"NA"
    fwrite(df_out,file.path(paths$output,"output","result",paste("atl-",atlas,"_zfp.csv",sep="")),row.names = F)
  }
  fwrite(df_concat,file.path(paths$output,"output","result","zfp.csv"),row.names=F)
}

#**************************************************
# GLM of Fingerprint change =======================
#**************************************************
gam_fp_core<-function(paths,df_fp,atlas,param,list_covar,list_mod,list_term,idx_var,
                      calc_parallel,test_mod){
  # Prepare clinical data, calculate diff and mean, and standardize
  data_clin<-func_clinical_data_long(paths,c(1,2),param$subset_subj,list_covar,rem_na_clin=T,prefix=paste("var-",idx_var,sep=""),print_terminal=F)
  list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
  df_clin<-data_clin$df_clin
  df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  df_clin<-func_clinical_data_diffmean(df_clin,list_id_subj,list_covar)
  #df_clin<-data.frame(wave="2-1",func_demean_clin(df_clin,separate_sex=T)$df_clin)
  df_clin<-data.frame(wave="2-1",func_std_clin(df_clin,separate_sex=T)$df_clin)
  
  # Prepare fingerprint data
  df_fp<-df_fp[df_fp$group_1=="whole" & df_fp$group_2=="whole",]
  df_fp$group_1<-df_fp$gorup_2<-NULL
  df_fp_long<-data.frame()
  for (id_subj in list_id_subj){
    df_fp_subset<-df_fp[df_fp$from_ses==1 & df_fp$from_ID_pnTTC==id_subj
                        & df_fp$to_ses==2 & df_fp$to_ID_pnTTC==id_subj,]
    df_fp_subset<-data.frame("ID_pnTTC"=id_subj,df_fp_subset[,"z_r"])
    colnames(df_fp_subset)<-c("ID_pnTTC","value")
    df_fp_long<-rbind(df_fp_long,df_fp_subset)
  }
  
  df_join<-dplyr::inner_join(df_fp_long,df_clin,by="ID_pnTTC")
  data_gamm<-gamm_core4(df_join,list_mod_in=param$list_mod_tanner,list_sex_in=param$list_sex,calc_parallel_in=F,test_mod_in=F)
  fwrite(data_gamm$df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
  fwrite(data_gamm$df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)
  fwrite(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
}


gam_fp<-function(paths_=paths,list_atlas_=list_atlas,param=param_gam_fp){
  print("Starting gam_fp().")
  nullobj<-func_createdirs(paths_,str_proc="gam_fp()",copy_log=T,list_param=param)
  
  for (atlas in list_atlas_){
    print(paste("Loading FP, atlas:",atlas,sep=" "))
    df_fp<-as.data.frame(fread(file.path(paths_$input,"output","result",paste("atl-",atlas,"_fp.csv",sep=""))))
    
    for (idx_tanner in names(param$list_tanner)){
      print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gam_fp_core(paths_,df_fp,atlas,param,list_covar,
                  list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
                  calc_parallel=T,test_mod=F)
    }
  }
  
  print("Combining results.")
  list_var<-param$list_tanner
  func_combine_result(paths_,list_atlas_,list_var,"",list(list("measure"="")),c("gamm","gamm_anova","gamm_aic"))
  print("Finished gam_fc().")
}

#**************************************************
# Variance attribution ============================
#**************************************************
func_stratify<-function(df_fp){
  df_zr_group<-df_stat_zr_group<-NULL
  # Group similarity
  list_zr<-df_fp$z_r; strat<-"group"
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Group + sex similarity
  list_zr<-df_fp[df_fp$from_sex==df_fp$to_sex,"z_r"]; strat<-"group+sex"
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Group + wave similarity
  list_zr<-df_fp[df_fp$from_ses==df_fp$to_ses,"z_r"]; strat<-"group+wave"
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Group (+ sex) + individual similarity
  list_zr<-df_fp[df_fp$from_ID_pnTTC==df_fp$to_ID_pnTTC,"z_r"]; strat<-"group+individual"
  if (length(list_zr)==0){
    list_zr<-NA
  }
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Group + Tanner similarity
  list_zr<-df_fp[df_fp$from_tanner==df_fp$to_tanner,"z_r"]; strat<-"group+tanner"
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Group + individual + Tanner similarity
  list_zr<-df_fp[df_fp$from_ID_pnTTC==df_fp$to_ID_pnTTC & df_fp$from_tanner==df_fp$to_tanner,"z_r"]; strat<-"group+individual+tanner"
  if (length(list_zr)==0){
    list_zr<-NA
  }
  df_zr_group<-rbind(df_zr_group,data.frame(strat=strat,z_r=list_zr))
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=strat,mean_z_r=mean(list_zr),n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
  # Calculate relative mean_z_r
  #df_stat_zr_group$rel_mean_z_r<-df_stat_zr_group$mean_z_r/df_stat_zr_group[df_stat_zr_group$strat=="group+individual+tanner","mean_z_r"]
  df_stat_zr_group$rel_mean_z_r<-df_stat_zr_group$mean_z_r/max(df_stat_zr_group$mean_z_r,na.rm=T)
  df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=c("individual","wave","sex","tanner","individual+tanner"),mean_z_r=NA,n_z_r=NA,sd_z_r=NA,sem_z_r=NA,
                                                      rel_mean_z_r=c(df_stat_zr_group[df_stat_zr_group$strat=="group+individual","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                                     df_stat_zr_group[df_stat_zr_group$strat=="group+wave","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                                     df_stat_zr_group[df_stat_zr_group$strat=="group+sex","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                                     df_stat_zr_group[df_stat_zr_group$strat=="group+tanner","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                                     df_stat_zr_group[df_stat_zr_group$strat=="group+individual+tanner","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group+tanner","rel_mean_z_r"])))
  return(list("df_zr"=df_zr_group,"df_stat_zr"=df_stat_zr_group))
}

variance_fp<-function(paths_=paths,list_atlas_=list_atlas,param=param_variance_fp){
  print("Starting variance_fp().")
  nullobj<-func_createdirs(paths_,str_proc="variance_fp()",copy_log=T,list_param=param)
  
  df_zr<-df_stat_zr<-NULL
  for (atlas in list_atlas_){
    # Load fingerprint data
    df_fp<-as.data.frame(fread(file.path(paths_$input,"output","result",paste("atl-",atlas,"_fp.csv",sep=""))))
    
    # Create list of groups
    list_group<-unique(c(as.character(df_fp$group_1),as.character(df_fp$group_2)))
    if ("whole" %in% list_group){
      list_group<-c("whole",list_group[list_group!="whole"])
    }
    
    for (label_wave in names(param$list_wave)){
      print(paste("Calculating, atlas: ",atlas,", wave: ",label_wave,sep=""))
      wave_clin<-param$list_wave[[label_wave]]$clin
      wave_mri<-param$list_wave[[label_wave]]$mri
      # Load and subset clinical data according to specified subsetting condition and covariate availability
      if (label_wave=="long"){
        df_clin<-func_clinical_data_long(paths,wave_clin,param$subset_sub,list_covar=param$list_covar,rem_na_clin=T,prefix=paste("wav-",label_wave,"_src",sep=""),print_terminal=F)$df_clin
      }else{
        # QC subsetting condition must accord with MRI wave, but under the name of clinical wave
        subset_subj<-param$subset_subj[wave_mri]
        names(subset_subj)<-wave_clin
        df_clin<-func_clinical_data_long(paths,wave_clin,subset_subj,list_covar=param$list_covar,rem_na_clin=T,prefix=paste("wav-",label_wave,"_src",sep=""),print_terminal=F)$df_clin
        # "wave" column must accord with MRI wave for later joining with FP data
        df_clin$wave<-wave_mri
      }
      df_clin<-dplyr::rename(df_clin,"ses"="wave")
      
      for (idx_group_1 in seq(length(list_group))){
        for (idx_group_2 in seq(idx_group_1,length(list_group))){
          group_1<-list_group[idx_group_1]; group_2<-list_group[idx_group_2]
          
          # Prepare dataframe for variance calculation
          df_fp_subset<-df_fp[df_fp$group_1==group_1 & df_fp$group_2==group_2,]
          df_fp_subset<-inner_join(df_fp_subset,df_clin,by=c("from_ID_pnTTC"="ID_pnTTC","from_ses"="ses"))
          df_fp_subset<-dplyr::rename(df_fp_subset,"from_sex"="sex","from_tanner"="tanner")
          df_fp_subset<-inner_join(df_fp_subset,df_clin,by=c("to_ID_pnTTC"="ID_pnTTC","to_ses"="ses"))
          df_fp_subset<-dplyr::rename(df_fp_subset,"to_sex"="sex","to_tanner"="tanner")
          
          df_head<-data.frame(atlas=atlas,wave=label_wave,group_1=group_1,group_2=group_2)
          # Male + Female
          data_stratify<-func_stratify(df_fp_subset)
          #df_zr<-rbind(df_zr,cbind(df_head,sex="all",data_stratify$df_zr))
          df_stat_zr<-rbind(df_stat_zr,cbind(df_head,sex="all",data_stratify$df_stat_zr))
          # Male
          data_stratify<-func_stratify(df_fp_subset[df_fp_subset$from_sex==1 & df_fp_subset$to_sex==1,])
          #df_zr<-rbind(df_zr,cbind(df_head,sex="male",data_stratify$df_zr))
          df_stat_zr<-rbind(df_stat_zr,cbind(df_head,sex="male",data_stratify$df_stat_zr))
          # Female
          data_stratify<-func_stratify(df_fp_subset[df_fp_subset$from_sex==2 & df_fp_subset$to_sex==2,])
          #df_zr<-rbind(df_zr,cbind(df_head,sex="female",data_stratify$df_zr))
          df_stat_zr<-rbind(df_stat_zr,cbind(df_head,sex="female",data_stratify$df_stat_zr))
        }
      } # End of loop over groups
    } # End of loop over waves
  } # End of loop over atlas
  
  # Save results
  #fwrite(df_zr,file.path(paths_$output,"output","result","zr.csv"),row.names=F)
  fwrite(df_stat_zr,file.path(paths_$output,"output","result","stat_zr.csv"),row.names=F)
  
  # Graphical output
  df_stat_zr$group_pair<-paste(as.character(df_stat_zr$group_1),as.character(df_stat_zr$group_2),sep="-")
  df_stat_zr$label_group_pair<-paste(str_to_title(gsub("_"," ",as.character(df_stat_zr$group_1))),
                                     str_to_title(gsub("_"," ",as.character(df_stat_zr$group_2))),sep=" - ")
  for (atlas in list_atlas_){
    list_group_pair<-unique(df_stat_zr[df_stat_zr$atlas==atlas,"group_pair"])
    list_label_group_pair<-unique(df_stat_zr[df_stat_zr$atlas==atlas,"label_group_pair"])
    for (label_wave in names(param$list_wave)){
      print(paste("Preparing plot, atlas: ",atlas,", wave: ",label_wave,sep=""))
      wave_clin<-param$list_wave[[label_wave]]$clin
      wave_mri<-param$list_wave[[label_wave]]$mri
      
      # Male + Female absolute
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="all",levels=c("group","group+sex","group+wave","group+individual"),
                    palette="Greys",list_group_pair,list_label_group_pair,type="abs")
      # Male + Female relative
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="all",levels=c("group","sex","individual"),
                    palette="Greys",list_group_pair,list_label_group_pair,type="rel")
      # Male absolute
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="male",levels=c("group","group+wave","group+individual","group+tanner","group+individual+tanner"),
                    palette="Blues",list_group_pair,list_label_group_pair,type="abs")
      # Male relative
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="male",levels=c("group","tanner","individual+tanner"),
                    palette="Blues",list_group_pair,list_label_group_pair,type="rel")
      # Female absolute
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="female",levels=c("group","group+wave","group+individual","group+tanner","group+individual+tanner"),
                    palette="Reds",list_group_pair,list_label_group_pair,type="abs") 
      # Female relative
      plot_variance(paths_,df_stat_zr,atlas=atlas,wave=label_wave,sex="female",levels=c("group","tanner","individual+tanner"),
                    palette="Reds",list_group_pair,list_label_group_pair,type="rel")
    }
  }
  print("Finished variance_fp()")
}


#**************************************************
# GLM and ANCOVA of Fingerprint change ============
#**************************************************

glm_core<-function(data_src){
  df_src<-data_src$df_src
  atlas<-data_src$atlas
  measure<-data_src$measure
  group_1<-data_src$group_1
  group_2<-data_src$group_2
  suffix<-data_src$suffix
  
  #print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,", GLM/GAM.",  sep=""))
  df_out_aic_add<-df_out_lm_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    list_plot<-list()
    list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
    for (idx_sex in list_sex){
      df_src_sex<-df_src[df_src$sex==idx_sex,]
      mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML")
      p_table<-summary.gam(mod)$p.table
      if (is.null(summary.gam(mod)$s.table)){
        df_out_lm_add_add<-data.frame(sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                      estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                      t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
        
      }else{
        s_table<-summary.gam(mod)$s.table
        df_out_lm_add_add<-rbind(data.frame(sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                            estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                            t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                                 data.frame(sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                            estimate=NA,se=NA,F=s_table[,'F'],
                                            t=NA,p=s_table[,'p-value']))
      }
      df_out_lm_add<-rbind(df_out_lm_add,df_out_lm_add_add)
      df_out_aic_add<-rbind(df_out_aic_add,
                            data.frame(sex=idx_sex,model=idx_mod,aic=mod$aic,aic_best_among_models=0))
      
      # Graphical output of GLM results
      for (idx_graph in names(list_graph_)){
        if (list_graph_[[idx_graph]][["x_axis"]] %in% colnames(mod$model)){
          # Add sex-wise lines/plots to existent plot, initialize if absent
          plot<-plot_gamm(plot_in=list_plot[[idx_graph]],mod_gamm=mod,
                          df_in=df_src_sex,
                          spec_graph=list_graph_[[idx_graph]])
          list_plot[[idx_graph]]<-plot
          
          # Output
          if (idx_sex==list_sex[length(list_sex)]){
            # Prepare x axis label
            axis_x<-list_graph_[[idx_graph]][["x_axis"]]
            for (idx_prefix in list(c("",""),c("ses1_"," 1st wave"),c("ses2_"," 2nd wave"),
                                    c("diff_"," difference"),c("mean_"," mean"))){
              for (idx_covar in names(list_covar_)){
                if (axis_x==paste(idx_prefix[1],idx_covar,sep="")){
                  label_x<-paste(list_covar_[[idx_covar]][["label"]],idx_prefix[2],sep='')
                }
              }
            }
            
            plot<-(plot
                   + ggtitle(paste(list_graph_[[idx_graph]][["title"]],atlas,measure,group_1,group_2,idx_mod,sep=" "))
                   + xlab(label_x)
                   + ylab("Fingerprint correlation")
                   + theme(legend.position = "none"))
            filename_plot<-paste("atl-",atlas,"_msr-",measure,"_",suffix,"_grp1-",group_1,"_grp2-",group_2,"_mod-",idx_mod,
                                 "_plt-",idx_graph,"_fp_glm.eps",sep="")
            ggsave(filename_plot,plot=plot,device=cairo_ps,
                   path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
          }
        }
      }
    }
  }
  
  # Compare AICs of GLM models
  df_out_aic_add_sex_rbind<-data.frame()
  for (idx_sex in list_sex){
    df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
    df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                       'aic_best_among_models']<-1
    df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
  }
  
  df_out_lm_add<-cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,df_out_lm_add)
  df_out_aic_add_sex_rbind<-cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                  df_out_aic_add_sex_rbind)
  
  return(list("df_out_lm_add"=df_out_lm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind))
}

ancova_core<-function(data_input){
  atlas<-data_input$atlas
  measure<-data_input$measure
  group_1=data_input$group_network[1]
  group_2=data_input$group_network[2]
  group_tanner<-data_input$group_tanner
  group_tanner_content<-data_input$group_tanner_content
  idx_sex<-data_input$group_sex
  df_src_ancova<-data_input$df_src_ancova
  
  # Calculate ANCOVA
  if (idx_sex=="all"){
    mod_ancova<-aov(value~long_tanner+diff_age+sex,data=df_src_ancova)
    color_plot<-"seagreen"
  }else if (idx_sex=="male"){
    mod_ancova<-aov(value~long_tanner+diff_age,data=df_src_ancova)
    color_plot<-"steelblue2"
  }else if (idx_sex=="female"){
    mod_ancova<-aov(value~long_tanner+diff_age,data=df_src_ancova)
    color_plot<-"lightcoral"
  }
  df_ancova<-summary(mod_ancova)[[1]]
  df_out_ancova<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,sex=idx_sex,test="ANCOVA",
                            term=rownames(df_ancova),label=NA,
                            p=df_ancova[,'Pr(>F)'],t=NA,F=df_ancova[,'F value'],
                            fit=NA,diff=NA,sigma=NA)
  
  # Extract fitting of ANCOVA
  fit_ancova<-mod_ancova$coefficients
  list_term<-names(mod_ancova$model)
  list_term<-c("(Intercept)",list_term[list_term!="value"])
  #list_term<-list_term[list_term!="value"]
  df_out_fit<-NULL
  for (term in list_term){
    if (term %in% names(mod_ancova$xlevels)){
      list_label<-mod_ancova$xlevels[[term]]
    }else{
      list_label<-NA
    }
    
    for (label in list_label){
      if (is.na(label)){
        term_label<-term
      }else{
        term_label<-paste(term,label,sep='')
      }
      if (term_label %in% names(mod_ancova$coefficients)){
        fit<-mod_ancova$coefficients[[term_label]]
      }else{
        fit<-0
      }
      df_out_fit<-rbind(df_out_fit,
                        data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,sex=idx_sex,test="fit",
                                   term=term,label=label,
                                   p=NA,t=NA,F=NA,
                                   fit=fit,diff=NA,sigma=NA))
    }
  }
  
  # Graphical output of ANCOVA tanner result
  mean_fit<-df_out_fit[df_out_fit$term=="(Intercept)","fit"]
  for (term in list_term[list_term!="(Intercept)" & list_term!="long_tanner"]){
    list_label<-df_out_fit[df_out_fit$term==term,"label"]
    if (is.na(list_label[1])){
      mean_fit<-mean_fit+mean(df_src_ancova[,term]*df_out_fit[df_out_fit$term==term,"fit"])
    }else{
      list_label<-sort(unique(list_label))
      df_calc_mean<-data.frame(df_out_fit[df_out_fit$term==term & df_out_fit$label %in% list_label,c("label","fit")])
      df_src_ancova_mean<-df_src_ancova
      colnames(df_src_ancova_mean)[colnames(df_src_ancova_mean)==term]<-"label"
      df_src_ancova_mean$label<-as.character(df_src_ancova_mean$label)
      df_calc_mean<-left_join(df_src_ancova_mean,df_calc_mean,by="label")
      mean_fit<-mean_fit+mean(df_calc_mean$fit)
    }
  }
  
  df_ancova_plot<-data.frame(matrix(nrow=length(group_tanner_content[["1"]]),
                                    ncol=length(group_tanner_content[["2"]])))
  rownames(df_ancova_plot)<-names(group_tanner_content[["1"]])
  colnames(df_ancova_plot)<-names(group_tanner_content[["2"]])
  for(group_tanner_1 in names(group_tanner_content[["1"]])){
    for(group_tanner_2 in names(group_tanner_content[["2"]])){
      mean_fit_tanner<-df_out_fit[df_out_fit$term=="long_tanner" & df_out_fit$label==paste(group_tanner_1,group_tanner_2,sep="_"),"fit"]
      if (length(mean_fit_tanner)>0){
        df_ancova_plot[group_tanner_1,group_tanner_2]<-mean_fit+mean_fit_tanner
      }else{
        df_ancova_plot[group_tanner_1,group_tanner_2]<-NA
      }
    }
  }
  plot_ancova<-plot_cor_heatmap(input=df_ancova_plot)
  suppressMessages(plot_ancova<-(plot_ancova
                                 + scale_fill_gradient(low="white",high=color_plot,name="r")
                                 + ggtitle(paste("FP Cor ANCOVA,",atlas,measure,group_1,group_2,group_tanner,idx_sex,sep=" "))
                                 + xlab("2nd wave Tanner stage")
                                 + ylab("1st wave Tanner stage")
                                 + theme(plot.title = element_text(hjust = 0.5),
                                         axis.text.x = element_text(size=12,angle = 0,vjust=0,hjust=0.5),
                                         axis.text.y = element_text(size=12))))
  
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp1-",group_1,"_grp2-",group_2,"_str-",group_tanner,"_sex-",idx_sex,"_fp_ancova.eps",sep=""),plot=plot_ancova,device=cairo_ps,
         path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
  
  # Calculate Tukey-Kramer
  tk<-summary(glht(mod_ancova, linfct = mcp('long_tanner' = 'Tukey')))$test
  df_out_posthoc<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,
                             sex=idx_sex,test="Tukey-Kramer",
                             term="long_tanner",label=names(tk$coefficients),
                             p=tk$pvalues[1:length(tk$coefficients)],t=tk$tstat,F=NA,
                             fit=NA,diff=tk$coefficients,sigma=tk$sigma)
  df_out<-rbind(df_out_ancova,df_out_fit,df_out_posthoc)
  return(df_out)
}

heatmap_gammfp<-function(paths_,df_out_lm,list_atlas,list_mod,suffix){
  
  for (atlas in list_atlas){
    df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas,]
    list_measure<-as.character(sort(unique(df_out_lm_subset$measure)))
    for (measure in list_measure){
      df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas & df_out_lm$measure==measure,]
      for (model in names(list_mod)){
        for (idx_sex in c(1,2)){
          df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                      & df_out_lm$measure==measure
                                      & df_out_lm$sex==idx_sex
                                      & df_out_lm$model==model,]
          list_term<-sort(unique(as.character(df_out_lm_subset$term)))
          list_term<-list_term[list_term!="(Intercept)"]
          for (term in list_term){
            df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                        & df_out_lm$measure==measure
                                        & df_out_lm$sex==idx_sex
                                        & df_out_lm$model==model
                                        & df_out_lm$term==term,]
            list_group<-sort(unique(c(as.character(df_out_lm_subset$group_1),
                                      as.character(df_out_lm_subset$group_2))))
            list_group<-list_group[list_group!="whole"]
            n_group<-length(list_group)
            
            if(n_group>0){
              mat_static<-data.frame(matrix(nrow=n_group,ncol=n_group))
              mat_pval<-data.frame(matrix(nrow=n_group,ncol=n_group))
              colnames(mat_static)<-rownames(mat_static)<-colnames(mat_pval)<-rownames(mat_pval)<-list_group
              for (idx_group_1 in seq(n_group)){
                for (idx_group_2 in seq(idx_group_1,n_group)){
                  group_1<-list_group[idx_group_1]
                  group_2<-list_group[idx_group_2]
                  df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                              & df_out_lm$measure==measure
                                              & df_out_lm$sex==idx_sex
                                              & df_out_lm$model==model
                                              & df_out_lm$term==term
                                              & df_out_lm$group_1==group_1
                                              & df_out_lm$group_2==group_2,c("estimate","F","p")]
                  if (nrow(df_out_lm_subset)>0){
                    name_static<-colnames(df_out_lm_subset[,c("estimate","F")])[!is.na(df_out_lm_subset[1,c("estimate","F")])]
                    mat_static[group_1,group_2]<-mat_static[group_2,group_1]<-df_out_lm_subset[[1,name_static]]
                    #mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-df_out_lm_subset[[1,"p"]]
                    pval<-df_out_lm_subset[[1,"p"]]
                    if (pval<0.001){
                      mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-"**"
                    }else if (pval<0.05){
                      mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-"*"
                    }else{
                      mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-""
                    }
                  }else{
                    mat_static[group_1,group_2]<-mat_static[group_2,group_1]<-NA
                    mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-""
                  }
                }
              }
              
              #df_pval<-rownames_to_column(mat_pval,"row")
              #df_pval<-gather(df_pval,key=column,value=p,2:ncol(df_pval))
              if (idx_sex==1){
                label_sex<-"male"
              }else{
                label_sex<-"female"
              }
              #mat_pval<-round(mat_pval,3)
              plot_stat<-plot_cor_heatmap(mat_static,mat_pval)
              suppressMessages(plot_stat<-(plot_stat
                                           + scale_fill_gradientn(colors = matlab.like2(100),
                                                                  lim=c(-max(max(mat_static),-min(mat_static)),max(max(mat_static),-min(mat_static))),
                                                                  name=name_static)
                                           + ggtitle(paste("GLM-FP,",atlas,measure,label_sex,model,term,sep=" "))
                                           + theme(plot.title = element_text(hjust = 0.5),
                                                   axis.title=element_blank())))
              ggsave(paste("atl-",atlas,"_msr-",measure,"_",suffix,"_sex-",label_sex,"_mod-",model,
                           "_plt-",term,"_fp_glm_heatmap.eps",sep=""),plot=plot_stat,device=cairo_ps,
                     path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
            }
          }
        }
      }
    }
  }
}

model_fp_multi<-function(paths_=paths,subset_subj_=model_fp_subset_subj,list_atlas_=list_atlas,
                         list_wave_=list_wave,skip_ancova=T,
                         list_covar_tanner_=model_fp_list_covar_tanner,list_tanner_=model_fp_list_tanner,
                         list_mod_tanner_=model_fp_list_mod_tanner,list_graph_tanner_=model_fp_list_graph_tanner,
                         list_strat_tanner_=model_fp_list_strat_tanner,
                         list_covar_hormone_=model_fp_list_covar_hormone,list_hormone_=model_fp_list_hormone,
                         list_mod_hormone_=model_fp_list_mod_hormone,list_graph_hormone_=model_fp_list_graph_hormone
                         ){
  print("Starting model_fp_multi().")
  nullobj<-func_createdirs(paths_,str_proc="model_fp_multi()")
  
  df_out_lm<-df_out_aic<-NULL
  # Loop over clinical variables
  #1 Tanner stage
  for (idx_tanner in names(list_tanner_)){
    print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
    list_covar<-list_covar_tanner_
    list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
    data_model<-model_fp(paths_,list_atlas_,list_wave_,list_covar_=list_covar,
                         list_mod_=list_mod_tanner_,list_graph_=list_graph_tanner_,
                         subset_subj_,list_strat_tanner_,skip_ancova,suffix=paste("var-",idx_tanner,sep=""))
    df_out_lm<-rbind(df_out_lm,cbind(variable=idx_tanner,data_model$df_out_lm))
    df_out_aic<-rbind(df_out_aic,cbind(variable=idx_tanner,data_model$df_out_aic))
  } # Finished looping over Tanner stages
  
  #2 Hormones
  for (idx_hormone in names(list_hormone_)){
    print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
    list_covar<-list_covar_hormone_
    list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
    data_model<-model_fp(paths_,list_atlas_,list_wave_,list_covar_=list_covar,
                         list_mod_=list_mod_hormone_,list_graph_=list_graph_hormone_,
                         subset_subj_,NULL,skip_ancova,suffix=paste("var-",idx_hormone,sep=""))
    df_out_lm<-rbind(df_out_lm,cbind(variable=idx_hormone,data_model$df_out_lm))
    df_out_aic<-rbind(df_out_aic,cbind(variable=idx_hormone,data_model$df_out_aic))
  } # Finished looping over Hormones
  
  write.csv(df_out_lm,file.path(paths_$output,"output","fp_glm.csv"),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output","fp_glm_aic.csv"),row.names = F)
  print("Finished model_fp_multi().")
}

model_fp<-function(paths_=paths,
                   list_atlas_=list_atlas,list_wave_=list_wave,list_covar_=list_covar,
                   list_mod_=list_mod,list_graph_=list_graph,subset_subj_=subset_subj,
                   list_strat_tanner_=list_strat_tanner,
                   skip_ancova=T,suffix=""
                   ){

  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T,prefix=suffix,print_terminal=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  list_src_glm<-list_src_ancova<-df_src_glm_bind<-NULL
  for (atlas in list_atlas_){
    # Load fingerprint data]
    print(paste("Loading FP, atlas:",atlas,sep=" "))
    df_fp<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep=""))))
    
    list_measure<-as.character(sort(unique(df_fp$measure)))

    for (measure in list_measure){
      print(paste("Preparing dataset, atlas: ",atlas,", measure: ",measure,sep=""))
      df_fp_meas<-df_fp[df_fp$measure==measure,]
    
      # Create list of subjects who meet subsetting condition and whose MRI data exist
      list_ses_exist <- sort(unique(c(df_fp_meas$from_ses,df_fp_meas$to_ses)))
      list_id_subj_exist<-list()
      for (ses in list_ses_exist){
        id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                         df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
        id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
        id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
        list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
      }
      
      # Identify subjects with longitudinal data
      list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                               list_id_subj_exist[["2"]]))
      n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
      
      # Create list of existing ROI subgroups
      list_group<-sort(unique(c(as.character(df_fp_meas$group_1),as.character(df_fp_meas$group_2))))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
      #list_group<-"whole      # Disable subgroup-wise analysis
      
      n_group<-length(list_group)
      for (idx_group_1 in seq(n_group)){
        group_1<-list_group[idx_group_1]
        df_fp_meas_grp1<-df_fp_meas[df_fp_meas$group_1==group_1,]
        for (idx_group_2 in seq(idx_group_1,n_group)){
          group_2<-list_group[idx_group_2]
          df_fp_meas_grp2<-df_fp_meas_grp1[df_fp_meas_grp1$group_2==group_2,]
          df_cor_fp<-data.frame(ID_pnTTC=list_id_subj_exist_twice)
          for (id_subj in list_id_subj_exist_twice){
            if (any(df_fp_meas_grp2$from_ID_pnTTC==id_subj & df_fp_meas_grp2$to_ID_pnTTC==id_subj)){
              df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-df_fp_meas_grp2[df_fp_meas_grp2$from_ID_pnTTC==id_subj
                                                                              & df_fp_meas_grp2$to_ID_pnTTC==id_subj,
                                                                              "z_r"]
            }else{
              df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-NA
            }
          }
          
          # Subset those without longitudinal fp correlation
          list_id_subj_nonna<-df_cor_fp[!is.na(df_cor_fp$value),"ID_pnTTC"]
          df_cor_fp<-df_cor_fp[df_cor_fp$ID_pnTTC %in% list_id_subj_nonna,]
          n_id_subj_exist_twice<-length(list_id_subj_nonna)
          #print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,
          #            ", ",as.character(n_id_subj_exist_twice)," subjects with longitudinal non-NA data.",sep=""))
          
          if (n_id_subj_exist_twice>0){
            # Create dataframe for GLM analysis
            df_src_glm<-func_clinical_data_diffmean(df_src=df_clin,
                                                list_id_subj=list_id_subj_nonna,
                                                list_covar=list_covar_)
            df_src_glm<-inner_join(df_src_glm,df_cor_fp,by="ID_pnTTC")
            df_src_glm$ID_pnTTC<-as.factor(df_src_glm$ID_pnTTC)
            df_src_glm$sex<-as.factor(df_src_glm$sex)
            
            df_src_glm_bind<-rbind(df_src_glm_bind,
                                   cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                         df_src_glm))
            
            list_src_glm<-c(list_src_glm,
                            list(list("df_src"=df_src_glm,"atlas"=atlas,"measure"=measure,
                                      "group_1"=group_1,"group_2"=group_2,"suffix"=suffix)))
          }
          
          # Calculate GLM
          #out_glm<-glm_core(df_src=df_join_grp,atlas,measure,group_1,group_2,
          #                  list_mod_,list_graph_,list_covar_,paths_)
          #df_out_lm<-rbind(df_out_lm,out_glm$df_out_lm_add)
          #df_out_aic<-rbind(df_out_aic,out_glm$df_out_aic_add)
          
          # Prepare ANCOVA dataset for later parallel computing
          if (!skip_ancova){
            print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,", ANCOVA preparation.",  sep=""))
            # Create list of input dataframes for parallel ANCOVA calculation
            for (group_tanner in names(list_strat_tanner_)){
              # group by longitudinal Tanner stage
              df_join_grp_tanner<-df_join_grp
              for (ses in c(1,2)){
                list_strat_tanner_ses<-names(list_strat_tanner_[[group_tanner]][[as.character(ses)]])
                for (label_tanner in list_strat_tanner_ses){
                  list_strat_tanner_ses_group<-list_strat_tanner_[[group_tanner]][[as.character(ses)]][[label_tanner]]
                  #print(list_strat_tanner_ses_group)
                  df_join_grp_tanner[df_join_grp_tanner[[paste('ses',as.character(ses),'_tanner',sep='')]] %in% list_strat_tanner_ses_group,
                                     paste('ses',as.character(ses),'_tanner_label',sep='')]<-label_tanner
                }
              }
              df_join_grp_tanner$long_tanner<-paste(as.character(df_join_grp_tanner$ses1_tanner_label),
                                                    as.character(df_join_grp_tanner$ses2_tanner_label),sep="_")
              df_join_grp_tanner$long_tanner<-as.factor(df_join_grp_tanner$long_tanner)
              
              list_sex<-list("all"=c(1,2),"male"=1,"female"=2)
              
              for (idx_sex in names(list_sex)){
                df_join_grp_tanner_sex<-df_join_grp_tanner[df_join_grp_tanner$sex %in% list_sex[[idx_sex]],]
                list_src_ancova<-c(list_src_ancova,
                                   list(list("atlas"=atlas,"measure"=measure,"group_network"=c(group_1,group_2),
                                             "group_tanner"=group_tanner,"group_tanner_content"=list_strat_tanner_[[group_tanner]],
                                             "group_sex"=idx_sex,"df_src_ancova"=df_join_grp_tanner_sex)))
              }
            }
          }
        }
      } # Finished looping over ROI subgroups
    } # Finished looping over measure
  } # Finished looping over atlas
  
  write.csv(df_src_glm_bind,file.path(paths_$output,"output",
                                      paste(suffix,"_fp_glm_src.csv",sep="")),row.names = F)

  # Parallel GLM calculation
  n_cluster<-floor(detectCores()*3/4)
  print(paste("Calculating GLM/GAM in parallel, threads: ",as.character(n_cluster),sep=""))
  clust<-makeCluster(n_cluster)
  clusterExport(clust,
                varlist=c("list_mod_","list_graph_","list_covar_","paths_","as.numeric.factor",
                          "gam","as.formula","summary.gam","plot_gamm","ggplot",
                          "predict.gam","geom_line","aes","geom_ribbon","geom_point",
                          "theme_light","element_text","ggtitle",
                          "xlab","ylab","theme","ggsave"),
                envir=environment())
  list_dst_glm<-pblapply(list_src_glm,glm_core,cl=clust)
  stopCluster(clust)
  
  df_out_lm<-df_out_aic<-NULL
  for (dst_glm in list_dst_glm){
    df_out_lm<-rbind(df_out_lm,dst_glm$df_out_lm_add)
    df_out_aic<-rbind(df_out_aic,dst_glm$df_out_aic_add)
  }
  
  # GLM/GAM Data saving
  rownames(df_out_lm)<-rownames(df_out_aic)<-NULL
  write.csv(df_out_lm, file.path(paths_$output,"output",paste(suffix,"_fp_glm.csv",sep="")),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output",paste(suffix,"_fp_glm_aic.csv",sep="")),row.names = F)
  
  
  # group-wise GLM/GAM heatmap visualization
  heatmap_gammfp(paths_,df_out_lm,list_atlas_,list_mod_,suffix=suffix)
  
  # Parallel ANCOVA calculation
  if (!skip_ancova){
    print("Calculating ANCOVA in parallel.")
    n_cluster<-min(floor(detectCores()*3/4),length(list_src_ancova))
    clust<-makeCluster(n_cluster)
    clusterExport(clust,
                  varlist=c("aov","summary","glht","mcp","left_join","paths_",
                            "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave","scale_fill_gradient","xlab","ylab"),
                  envir=environment())
    list_df_ancova<-pblapply(list_src_ancova,ancova_core,cl=clust)
    stopCluster(clust)
    df_out_ancova<-NULL
    for (df_ancova in list_df_ancova){
      df_out_ancova<-rbind(df_out_ancova,df_ancova)
    }
    # Data saving
    rownames(df_out_ancova)<-NULL
    write.csv(df_out_ancova,file.path(paths_$output,"output","fp_ancova.csv"),row.names=F)
  }

  return(list("df_out_lm"=df_out_lm,"df_out_aic"=df_out_aic,"df_src"=df_src_glm_bind))
}


#**************************************************
# Fingerprint identification ======================
#**************************************************
identify_fp<-function(paths_=paths,list_atlas_=list_atlas,list_wave_=list_wave,
                      #list_covar_=list_covar,
                      subset_subj_=subset_subj,n_permutation_=n_permutation
                      ){
  print("Starting identify_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar=NULL,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  df_out_combined<-NULL
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep="")))
    
    list_measure<-sort(unique(df_fp$measure))
    
    for (measure in list_measure){
      print(paste("Atlas: ",atlas," Measure: ",measure,sep=""))
      df_fp_meas<-df_fp[df_fp$measure==measure,]
      
      # Create list of subjects who meet subsetting condition and whose MRI data exist
      list_ses_exist <- sort(unique(c(df_fp_meas$from_ses,df_fp_meas$to_ses)))
      list_id_subj_exist<-list()
      for (ses in list_ses_exist){
        id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                         df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
        id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
        id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
        list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
      }
      
      # Extract those with longitudinal data
      list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                               list_id_subj_exist[["2"]]))
      n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
      print(paste(as.character(n_id_subj_exist_twice)," subjects with two sessions.",sep=""))
      df_fp_exist_twice<-df_fp_meas[(df_fp_meas$from_ID_pnTTC %in% list_id_subj_exist_twice) & (df_fp_meas$to_ID_pnTTC %in% list_id_subj_exist_twice),]
      df_fp_exist_twice<-df_fp_exist_twice[(df_fp_exist_twice$from_ses==1 & df_fp_exist_twice$to_ses==2),]
      
      # Output subset with longitudinal data
      write.csv(df_fp_exist_twice,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_input_subset.csv",sep="")),row.names=F)
      
      list_group<-sort(unique(as.character(df_fp_exist_twice$group)))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
      
      # Correlation matrix graphical output
      for (group in list_group){
        df_fp_exist_twice_plot<-df_fp_exist_twice[df_fp_exist_twice$group==group,c('from_ID_pnTTC','to_ID_pnTTC','r')]
        df_fp_exist_twice_plot<-spread(df_fp_exist_twice_plot,key=to_ID_pnTTC,value=r)
        colnames(df_fp_exist_twice_plot)[-1]<-rownames(df_fp_exist_twice_plot)<-sprintf("%05d",df_fp_exist_twice_plot$from_ID_pnTTC)
        df_fp_exist_twice_plot<-df_fp_exist_twice_plot[-1]
        plot_fp_exist_twice<-plot_cor_heatmap(input=df_fp_exist_twice_plot)
        suppressMessages(plot_fp_exist_twice<-(plot_fp_exist_twice
                                               + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                               + ggtitle(paste("Long. FP Cor,",atlas,measure,group,sep=" "))
                                               + xlab("2nd wave")
                                               + ylab("1st wave")
                                               + theme(plot.title = element_text(hjust = 0.5))))
        ggsave(paste("atl-",atlas,"_msr-",measure,"_grp-",group,"_fp_id.eps",sep=""),plot=plot_fp_exist_twice,device=cairo_ps,
                     path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
      }
      
      # Calculate fingerprint identification
      df_ident<-df_perm<-NULL
      
      for (group in list_group){
        df_ident_grp<-data.frame("group"=group,"target"=list_id_subj_exist_twice)
        df_perm_grp<-data.frame("group"=group,"id_perm"=seq(1,n_permutation_))
        df_fp_exist_twice_group<-df_fp_exist_twice[df_fp_exist_twice$group==group,]
        for(ses in c(1,2)){
          if(ses==1){
            df_fp_pool<-df_fp_exist_twice_group[c("from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
          }else{
            df_fp_pool<-df_fp_exist_twice_group[c("to_ses","to_ID_pnTTC","from_ses","from_ID_pnTTC","r")]
          }
          colnames(df_fp_pool)<-c("target_ses","target_ID_pnTTC","pool_ses","pool_ID_pnTTC","r")
          for (id_subj in list_id_subj_exist_twice){
            df_fp_subset<-df_fp_pool[df_fp_pool$target_ID_pnTTC==id_subj,]
            list_id_subj_ordered<-df_fp_subset[order(df_fp_subset$r,decreasing = TRUE,na.last=NA),'pool_ID_pnTTC']
            if (id_subj %in% list_id_subj_ordered){
              rank_similarity<-which(list_id_subj_ordered==id_subj)
            }else{
              rank_similarity<-NA
            }
            df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_rank",sep='')]<-rank_similarity
            if (is.na(rank_similarity)){
              df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-0
            }else{
              if (rank_similarity==1){
                df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-1
              }else{
                df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-0
              }
            }
          }
          
          # Permutation test calculation
          for (i in seq(1,n_permutation_)){
            # Create dataframe for random shuffling
            df_rand<-data.frame(pool_ID_pnTTC=list_id_subj_exist_twice,rand_ID_pnTTC=sample(list_id_subj_exist_twice,length(list_id_subj_exist_twice)))
            
            # Randomize subject ID of pool session accordint to df_rand
            df_fp_rand<-left_join(df_fp_pool,df_rand,by="pool_ID_pnTTC")
            n_identified<-0
            for (id_subj in list_id_subj_exist_twice){
              df_fp_rand_subj<-df_fp_rand[df_fp_rand$target_ID_pnTTC==id_subj,]
              list_id_subj_ordered<-df_fp_rand_subj[order(df_fp_rand_subj$r,decreasing = TRUE,na.last=NA),'rand_ID_pnTTC']
              if (id_subj %in% list_id_subj_ordered){
                rank_similarity<-which(list_id_subj_ordered==id_subj)
              }else{
                rank_similarity<-NA
              }
              if (!is.na(rank_similarity)){
                if (rank_similarity==1){
                  n_identified<-n_identified+1
                }
              }
            }
            df_perm_grp[i,paste(as.character(ses),'_n_ident',sep='')]<-n_identified
            #print(paste("Iteration: ",as.character(i),", subjects identified: ",as.character(n_identified),sep=""))
          }
        }
        df_ident<-rbind(df_ident,df_ident_grp)
        df_perm<-rbind(df_perm,df_perm_grp)
        
        n_subj<-length(list_id_subj_exist_twice)
        n_id_1_tar<-sum(df_ident_grp["1_targeted_identification"])
        n_id_2_tar<-sum(df_ident_grp["2_targeted_identification"])
        n_id<-n_id_1_tar+n_id_2_tar
        prop_id_1_tar<-n_id_1_tar/n_subj
        prop_id_2_tar<-n_id_2_tar/n_subj
        prop_id<-n_id/(n_subj*2)
        p_perm_1_tar<-sum(df_perm_grp["1_n_ident"]>sum(df_ident_grp["1_targeted_identification"]))/n_permutation_
        p_perm_2_tar<-sum(df_perm_grp["2_n_ident"]>sum(df_ident_grp["2_targeted_identification"]))/n_permutation_
        p_perm<-(p_perm_1_tar+p_perm_2_tar)/2
        
        df_out_combined<-rbind(df_out_combined,
                               data.frame(atlas=atlas,measure=measure,group=group,n_subj=n_subj,n_identified=n_id,proportion_identified=prop_id,
                                          n_identified_1_targeted=n_id_1_tar,proportion_identifeid_1_targeted=prop_id_1_tar,
                                          n_identified_2_targeted=n_id_2_tar,proportion_identified_2_targeted=prop_id_2_tar,
                                          p_permutation=p_perm,
                                          p_permutation_1_targeted=p_perm_1_tar,p_permutation_2_targeted=p_perm_2_tar))
      }
      write.csv(df_ident,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_id.csv",sep="")),row.names=F)
      write.csv(df_perm,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_perm.csv",sep="")),row.names=F)
    }
  }
  write.csv(df_out_combined,file.path(paths_$output,"output","fp_id_summary.csv"),row.names=F)
  print("Finished identify_fp().")
}

