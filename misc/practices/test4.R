#source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')
source('D:/atiro/GitHub/MRI_Analysis/analyze/connection.R')

####

paths_=paths
list_atlas_=list_atlas
param<-param_gam_fc_diff

####

print("Starting gam_fc_diff().")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_diff()",copy_log=T,list_param=param)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave="2-1",include_grp=T,abs_nfc=param$abs_nfc)
data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

#gam_fc_diff_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
#                 param$list_mod_tanner,param$list_term_tanner,idx_tanner,
#                 calc_parallel=T,test_mod=F)

####

paths<-paths_
list_sex<-list(1,2)
list_mod<-param$list_mod_tanner
list_term<-param$list_term_tanner
idx_var<-idx_tanner
#calc_parallel=F
calc_parallel=T
test_mod=F
#test_mod=T

####
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

####
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

####

# Threshold and plot graph edges
if (file.exists(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep="")))){
  df_plot<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_plot.csv",sep=""))))
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

####


# Detect sub-network by breadth-first approach
if (nrow(df_plot)>0){
  df_net<-df_node<-df_size_net<-df_pred_ancova<-NULL
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
            #plot_subnet<-plot_circular2(df_edge=network$df_edge,df_node=network$df_node,df_roi=data_fc$df_roi,rule_order="degree")
            plot_subnet<-plot_net(df_edge=network$df_edge,df_node=network$df_node,df_roi=data_fc$df_roi)
            plot_subnet<-(plot_subnet+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,", expvar: ",var_exp,", sex: ",label_sex,", p value: p<",param$param_nbs$p_cdt_threshold,", #",as.character(idx_net),sep="")))
            ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,"_sex-",label_sex,"_pval-p_",param$param_nbs$p_cdt_threshold,"_idx-",as.character(idx_net),"_subnet.png",sep=""),
                   plot=plot_subnet,path=file.path(paths$output,"output","plot"),height=10,width=10,dpi=600)
            #filename_plot<-file.path(paths$output,"output","plot",paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,"_sex-",label_sex,"_pval-p_",param$param_nbs$p_cdt_threshold,"_idx-",as.character(idx_net),"_subnet.png",sep=""))
            #title_plot<-paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,", expvar: ",var_exp,", sex: ",label_sex,", p value: p<",param$param_nbs$p_cdt_threshold,", #",as.character(idx_net),sep="")
            #plot_subnet<-plot_lgl(df_edge=network$df_edge,df_node=network$df_node,df_roi=data_fc$df_roi,filename_plot,title_plot)
            if(idx_term %in% names(param$param_ancova_pred)){
              data_pred_ancova<-plot_pred_ancova(df_edge=network$df_edge,df_gamm=df_gamm,data_fc=data_fc,param_ancova_pred=param$param_ancova_pred,var_exp)
              df_pred_ancova<-rbind(df_pred_ancova,data.frame(id_net=idx_net,data_pred_ancova$df_plot))
              plot_pred<-(data_pred_ancova$plot+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", model: ",idx_mod,"\nexpvar: ",var_exp,", sex: ",label_sex,", p value: p<",param$param_nbs$p_cdt_threshold,", #",as.character(idx_net),sep="")))
              ggsave(paste("atl-",atlas,"_var-",idx_var,"_mod-",idx_mod,"_trm-",idx_term,"_sex-",label_sex,"_pval-p_",param$param_nbs$p_cdt_threshold,"_idx-",as.character(idx_net),"_pred.png",sep=""),
                     plot=plot_pred,path=file.path(paths$output,"output","plot"),height=5,width=5,dpi=600)
              
            }
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
  fwrite(df_pred_ancova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_pred.csv",sep="")),row.names = F)
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
df_threshold_size<-df_size_net_p<-NULL
for (idx_mod in param$param_nbs$list_mod){
  for (idx_term in param$param_nbs$list_term){
    var_exp<-list_term[[idx_term]][["var_exp"]]
    for (idx_sex in list_sex){
      list_max_size<-df_max_size[df_max_size$model==idx_mod & df_max_size$term==var_exp
                                 & df_max_size$sex==idx_sex,"max_size"]
      list_max_size<-sort(list_max_size)
      df_size_net_subset<-df_size_net[df_size_net$model==idx_mod & df_size_net$term==var_exp & df_size_net$sex==idx_sex,]
      if(nrow(df_size_net_subset)>0){
        for (idx_row in seq(nrow(df_size_net_subset))){
          df_size_net_subset[idx_row,"p"]<-sum(list_max_size>df_size_net_subset[idx_row,"size"])/param$param_nbs$n_perm
        }
        df_size_net_p<-rbind(df_size_net_p,df_size_net_subset)
      }
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
fwrite(df_size_net_p,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_bfs_size.csv",sep="")),row.names = F)
