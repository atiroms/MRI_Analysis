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

for (label_wave in names(param$list_wave)){
  wave_clin<-param$list_wave[[label_wave]]$clin
  wave_mri<-param$list_wave[[label_wave]]$mri

  
}

####

label_wave<-"c2m1"
wave_clin<-2
wave_mri<-1

####

# Prepare clinical data and demean
# QC subsetting condition must accord with MRI wave, but under the name of clinical wave
subset_subj<-param$subset_subj[wave_mri]
names(subset_subj)<-wave_clin
data_clin<-func_clinical_data_long(paths,wave_clin,subset_subj,list_covar,rem_na_clin=T,prefix=paste("atl-",atlas,"_var-",idx_var,"_wave-",label_wave,sep=""),print_terminal=F)
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
data_bfs<-func_detect_subnset(paths,df_plot,df_gamm,data_fc,
                              atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
df_net<-data_bfs$df_net; df_node<-data_bfs$df_node; df_size_net<-data_bfs$df_size_net; df_pred_ancova<-data_bfs$df_pred_ancova
# Permutation test
data_nbs<-func_nbs_permutation(paths,df_fc,df_clin,df_size_net,data_fc,calc_parallel,
                               atlas,param,list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
