source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')

####

paths_=paths
list_atlas_=list_atlas
param=param_gam_fc_cs

####

print("Starting gam_fc_cs().")
nullobj<-func_createdirs(paths_,str_proc="gam_fc()",copy_log=T,list_param=param)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
#gam_fc_cs_core(paths_,atlas,param,list_sex=list(1,2),list_covar,
#               list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
#               calc_parallel=T,test_mod=F)


####

paths<-paths_
list_sex<-list(1,2)
list_mod<-param$list_mod_tanner
list_term<-param$list_term_tanner
idx_var<-idx_tanner
calc_parallel=T
test_mod=F

####

label_wave<-"c2m1"
wave_clin<-2
wave_mri<-1

####
print(paste("Atlas: ",atlas,", Measure: ",idx_var,", Wave: ",label_wave,sep=""))
wave_clin<-param$list_wave[[label_wave]]$clin
wave_mri<-param$list_wave[[label_wave]]$mri

# Prepare clinical data and demean
# QC subsetting condition must accord with MRI wave, but under the name of clinical wave
subset_subj<-param$subset_subj[wave_mri]
names(subset_subj)<-wave_clin
data_clin<-func_clinical_data_long(paths,wave_clin,subset_subj,list_covar,rem_na_clin=T,prefix=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src",sep=""),print_terminal=F)
df_clin<-func_demean_clin(data_clin$df_clin,separate_sex=T)$df_clin
fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_clin.csv",sep="")),row.names=F)
df_clin$wave<-wave_mri # Need to meet MRI wave for later joining

# Prepare FC data
print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc2(paths,atlas,param$key_group,list_wave=wave_mri,include_grp=T,abs_nfc=param$abs_nfc)
df_fc<-data_fc$df_fc; df_fc_grp<-data_fc$df_fc_grp
fwrite(df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_fc.csv",sep="")),row.names=F)
fwrite(df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_src_fc_grp.csv",sep="")),row.names=F)

# Calculate model
data_gamm<-func_calc_gamm(paths,df_clin,df_fc,df_fc_grp,data_fc,
                          calc_parallel,test_mod,
                          atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova_grp<-data_gamm$df_anova_grp
# Threshold and plot graph edges
data_plot<-func_threshold_gamm(paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                               atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
df_plot<-data_plot$df_plot; df_plot_grp<-data_plot$df_plot_grp
# Detect sub-network by breadth-first approach
data_bfs<-func_detect_subnset(paths,df_plot,df_gamm,data_fc,plot_result=F,
                              atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
#data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,data_bfs,data_fc,calc_parallel,plot_result=T,
#                               atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)

####

plot_result=T

####


print("Calculating permutation")

# Prepare parallelization cluster
test_mod<-F
if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
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
      for (idx_sex in list_sex){
        df_clin_perm_add<-df_clin[df_clin$sex==idx_sex,]
        df_clin_perm_add[,var_exp_perm]<-sample(df_clin_perm_add[,var_exp_perm])
        df_clin_perm<-rbind(df_clin_perm,df_clin_perm_add)
      }
      # Join FC and permuted clinical data
      df_join<-join_fc_clin(df_fc,df_clin_perm)
      
      # Calculate model
      data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
      df_gamm<-data_gamm$df_gamm
      df_anova<-data_gamm$df_anova
      for (idx_mod in param$param_nbs$list_mod){
        list_term_detect<-set_term$term_detect
        for (idx_term_detect in list_term_detect){
          var_exp_detect<-list_term[[idx_term_detect]][["var_exp"]]
          for (idx_sex in list_sex){
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
              for (p_cdt in param$param_nbs$p_cdt_threshold){
                df_sign<-df_gamm_subset[df_gamm_subset$p<p_cdt,]
                max_size<-func_bfs(df_sign)$max_size
                df_max_size<-rbind(df_max_size,data.frame(id_perm=idx_perm,model=idx_mod,term=var_exp_detect,sex=idx_sex,p_threshold=p_cdt,max_size=max_size))
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
for (idx_mod in param$param_nbs$list_mod){
  for (set_term in param$param_nbs$list_term){
    for (idx_term_detect in set_term$term_detect){
      var_exp_detect<-list_term[[idx_term_detect]][["var_exp"]]
      for (idx_sex in list_sex){
        for (p_cdt in param$param_nbs$p_cdt_threshold){
          list_max_size<-df_max_size[df_max_size$model==idx_mod & df_max_size$term==var_exp_detect
                                     & df_max_size$p_threshold==p_cdt & df_max_size$sex==idx_sex,"max_size"]
          if (length(list_max_size)>0){
            list_max_size<-sort(list_max_size)
            df_size_net_subset<-df_size_net[df_size_net$model==idx_mod & df_size_net$term==var_exp_detect
                                            & df_size_net$p_threshold==p_cdt & df_size_net$sex==idx_sex,]
            if(nrow(df_size_net_subset)>0){
              for (idx_row in seq(nrow(df_size_net_subset))){
                p_fwe<-sum(list_max_size>df_size_net_subset[idx_row,"size"])/param$param_nbs$n_perm
                df_size_net_subset[idx_row,"p_fwe"]<-p_fwe
                if (p_fwe<param$param_nbs$p_perm_threshold){
                  if (plot_result){
                    idx_net<-as.character(df_size_net_subset[idx_row,"id_net"])
                    df_net_subset<-df_net[df_net$model==idx_mod & df_net$term==var_exp_detect & df_net$p_threshold==p_cdt & df_net$sex==idx_sex & df_net$id_net==idx_net,]
                    df_node_subset<-df_node[df_node$model==idx_mod & df_node$term==var_exp_detect & df_node$p_threshold==p_cdt & df_node$sex==idx_sex & df_node$id_net==idx_net,
                                            c("node","degree","label_node","group_node")]
                    plot_subnet<-plot_net(df_edge=df_net_subset,df_node=df_node_subset,df_roi=data_fc$df_roi)
                    plot_subnet<-(plot_subnet+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,", expvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", #",as.character(idx_net),sep="")))
                    list_plot<-c(list_plot,list(list("plot"=plot_subnet,"height"=15,"width"=15,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                     "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_idx-",as.character(idx_net),"_subnet.png",sep=""))))
                    if(idx_term_detect %in% names(param$param_ancova_pred)){
                      df_pred_ancova_subset<-df_pred_ancova[df_pred_ancova$model==idx_mod & df_pred_ancova$term==var_exp_detect & df_pred_ancova$p_threshold==p_cdt & df_pred_ancova$sex==idx_sex & df_pred_ancova$id_net==idx_net,]
                      #plot_pred<-plot_pred_ancova(df_edge=network$df_edge,df_gamm=df_gamm,data_fc=data_fc,param_ancova_pred=param$param_ancova_pred,idx_term_detect,var_exp_detect)
                      plot_pred<-(plot_pred_ancova(df_pred_ancova_subset)+ggtitle(paste("atlas: ",atlas,", measure: ",idx_var,", wave: ",label_wave,", model: ",idx_mod,"\nexpvar: ",var_exp_detect,", sex: ",label_sex,", p value: p<",p_cdt,", #",as.character(idx_net),sep=""))+ xlab(list_term[[idx_term_detect]][["title"]]))
                      list_plot<-c(list_plot,list(list("plot"=plot_pred,"height"=5,"width"=5,"dpi"=600,"path"=file.path(paths$output,"output","plot"),
                                                       "filename"=paste("atl-",atlas,"_var-",idx_var,"_wav-",label_wave,"_mod-",idx_mod,"_trm-",idx_term_detect,"_sex-",label_sex,"_pval-p_",p_cdt,"_idx-",as.character(idx_net),"_pred.png",sep=""))))
                    }
                  }  
                }
              }
              df_fwep<-rbind(df_fwep,df_size_net_subset)
            }
            thr_size_nbs<-list_max_size[ceiling(length(list_max_size)*(1-param$param_nbs$p_perm_threshold))]
            df_threshold_size<-rbind(df_threshold_size,data.frame(model=idx_mod,term=var_exp_detect,sex=idx_sex,p_threshold=p_cdt,
                                                                  thr_size=thr_size_nbs))
            if (idx_sex==1){
              label_sex<-"m";title_sex<-"male";color_plt<-"steelblue2"
            }else{
              label_sex<-"f";title_sex<-"female";color_plt<-"lightcoral"
            }
            title_plot<-list_term[[idx_term_detect]][["title"]]
            plot_permutation(paths,list_max=list_max_size,thr_size_nbs,
                             atlas,var=idx_var,wave=label_wave,idx_mod,idx_term_detect,label_sex,title_plot,title_sex,p_cdt,color_plt)
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


