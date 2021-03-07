

paths_=paths
list_atlas_=list_atlas
param=param_gamm_fc

####

print("Starting gamm_fc().")
nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T,list_param=param)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
#data_fc<-prep_data_fc(paths_,atlas,param$key_group,abs_nfc=param$abs_nfc)
data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=c("1","2"),include_grp=T,abs_nfc=param$abs_nfc)
data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

####

paths<-paths_
list_sex<-list(1,2)
list_mod<-param$list_mod_tanner
list_term<-param$list_term_tanner
idx_var<-idx_tanner
calc_parallel=T
test_mod=F

####

df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                 prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
df_clin<-func_demean_clin(df_clin,thr_cont=6,separate_sex=T)$df_clin # thr_cont=4 to demean Tanner, =6 not to
fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_clin.csv",sep="")),row.names=F)

# Prepare FC data
df_fc<-data_fc$df_fc
df_fc_grp<-data_fc$df_fc_grp
fwrite(df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc.csv",sep="")),row.names=F)
fwrite(df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc_grp.csv",sep="")),row.names=F)

####

# Join clinical and FC data
df_join<-join_fc_clin(df_fc,df_clin)
df_join_grp<-join_fc_clin(df_fc_grp,df_clin)

# Calculate model
if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                              "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC"),
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
fwrite(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
fwrite(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
fwrite(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)
fwrite(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
fwrite(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
fwrite(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep="")),row.names = F)

####

# Threshold and plot graph edges
print("Thresholding and plotting results")
data_plot<-func_threshold_gamm(paths,df_gamm,df_gamm_grp,df_anova,df_anova_grp,data_fc,
                               atlas,param,list_sex,list_covar,list_mod,list_term,idx_var)
df_plot<-data_plot$df_plot; df_plot_grp<-data_plot$df_plot_grp

# Detect sub-network by breadth-first approach
print("Detecting subnetworks")
data_bfs<-func_detect_subnset(paths,df_plot,df_gamm,data_fc,
                              atlas,param,list_sex,list_covar,list_mod,list_term,idx_var)
df_net<-data_bfs$df_net; df_node<-data_bfs$df_node; df_size_net<-data_bfs$df_size_net; df_pred_ancova<-data_bfs$df_pred_ancova

# Permutation test
print("Calculating permutation")
t_start<-Sys.time()
data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,df_size_net,data_fc,calc_parallel,
                               atlas,param,list_sex,list_covar,list_mod,list_term,idx_var)
print(Sys.time()-t_start) # 8.54651 mins for 10 permutations without batch

####
