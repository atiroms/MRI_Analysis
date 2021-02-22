#source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')

####

paths_=paths
list_atlas_=list_atlas
param<-param_gamm_fc

####

print("Starting gamm_fc().")
nullobj<-func_createdirs(paths_,str_proc="gamm_fc_()",copy_log=T)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,param$key_group)
#data_fc$df_edge$label_from<-data_fc$df_edge$label_to<-NULL
data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
#data_fc$df_edge_grp$label_from<-data_fc$df_edge_grp$label_to<-NULL
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

####

paths<-paths_
list_wave<-param$list_wave
list_sex<-list(1,2)
subset_subj<-param$subset_subj
list_mod<-param$list_mod_tanner
list_plot<-param$list_plot_tanner
idx_var<-idx_tanner
list_p<-param$list_p
calc_parallel=F
test_mod=F
#test_mod=T

####

clust<-makeCluster(1)
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","sort","gam","as.formula","summary.gam",
                              "anova.gam","as.numeric.factor"),
              envir=environment())

# Prepare clinical data and demean
df_clin<-func_clinical_data_long(paths,list_wave,subset_subj,list_covar,rem_na_clin=T,
                                 prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
df_clin<-func_demean_clin(df_clin,thr_cont=6,separate_sex=T)$df_clin # thr_cont=4 to demean Tanner, =5 not to

# ROI-wise GAMM of FC
df_join<-join_fc_clin(data_fc$df_fc,df_clin)
data_gamm<-iterate_gamm3(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_plot,calc_seed_level=F))
df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_plot,calc_seed_level=F))
write.csv(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
write.csv(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
write.csv(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)

# Group-wise GAMM of FC
df_join_grp<-join_fc_clin(data_fc$df_fc_grp,df_clin)
data_gamm_grp<-iterate_gamm3(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_plot,calc_seed_level=F))
df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_plot,calc_seed_level=F))
write.csv(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
write.csv(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep="")),row.names = F)

stopCluster(clust)

df_plot<-df_plot_grp<-NULL
for (idx_mod in names(list_mod)){
  for (idx_plot in names(list_plot)){
    var_exp<-list_plot[[idx_plot]][["var_exp"]]
    for (idx_sex in list_sex){
      # Subset GAMM result dataframe for plotting
      if (idx_sex==1){
        label_sex<-"m"
      }else{
        label_sex<-"f"
      }
      df_gamm_subset<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp & df_gamm$sex==idx_sex,]
      df_gamm_grp_subset<-df_gamm_grp[df_gamm_grp$model==idx_mod & df_gamm_grp$term==var_exp & df_gamm_grp$sex==idx_sex,]
      if (nrow(df_gamm_subset)>0){
        for (p in list_p){
          df_plot_subset<-df_gamm_subset[df_gamm_subset[[p$type]]<p$threshold,]
          df_plot_grp_subset<-df_gamm_grp_subset[df_gamm_grp_subset[[p$type]]<p$threshold,]
          plot_gamm<-plot_gam_fc3(df_plot_subset,df_plot_grp_subset,data_fc)
          plot_gamm<-annotate_figure(plot_gamm,
                                     top = text_grob(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,
                                                           ", expvar: ",var_exp,", sex: ",label_sex,", p value: ",p$type,"<",p$threshold,sep=""),
                                                     color = "black", size = 14))
          ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_plt-",var_exp,
                       "_sex-",label_sex,"_pval-",p$type,"_",p$threshold,
                       "_gamm.png",sep=""),
                 plot=plot_gamm,path=file.path(paths$output,"output","plot"),height=13,width=10,dpi=600)
          df_head<-data.frame(p_type=p$type,p_threshold=p$threshold)
          if (nrow(df_plot_subset)>0){
            df_plot<-rbind(df_plot,cbind(df_head,df_plot_subset))
          }
          if (nrow(df_plot_grp_subset)>0){
            df_plot_grp<-rbind(df_plot_grp,cbind(df_head,df_plot_grp_subset))
          }
        }
      }
    }
  }
}

# Save results
write.csv(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
write.csv(df_plot,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot.csv",sep="")),row.names = F)
write.csv(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
write.csv(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
write.csv(df_plot_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot_grp.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)