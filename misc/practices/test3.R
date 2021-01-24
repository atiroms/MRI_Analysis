paths_=paths
list_atlas_=list_atlas
key_group_='group_3'
subset_subj_=sex_diff_fc_cs_subset_subj
list_covar_=sex_diff_fc_cs_list_covar
list_mod_cs_=sex_diff_fc_cs_list_mod_cs
list_mod_diff_=sex_diff_fc_cs_list_mod_diff
list_plot_=sex_diff_fc_cs_list_plot
thr_p_cdt_=sex_diff_fc_cs_thr_p_cdt
thr_p_perm_=sex_diff_fc_cs_thr_p_perm
n_perm_=sex_diff_fc_cs_n_perm

####

print("Starting sex_diff_fc_cs()")
nullobj<-func_createdirs(paths_,str_proc="sex_diff_fc_cs()",copy_log=T)
# Increase memory limit
memory.limit(1000000)

####

atlas<-list_atlas[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,key_group_,include_diff=T)
data_clin<-func_clinical_data_long(paths_,list_wave=c("1","2"),subset_subj_,list_covar_,rem_na_clin=T,
                                   prefix=paste("atl-",atlas,sep=""),print_terminal=F)