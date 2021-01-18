paths_=paths
subset_subj_=sex_diff_fc_subset_subj
list_wave_=list_wave
list_atlas_=list_atlas
key_group_='group_3'
list_covar_=sex_diff_fc_list_covar
list_mod_=sex_diff_fc_list_mod
list_plot_=sex_diff_fc_list_plot
list_type_p_=sex_diff_fc_list_type_p
thr_p_=sex_diff_fc_thr_p

####

print("Starting sex_diff_fc_multi().")
nullobj<-func_createdirs(paths_,str_proc="sex_diff_fc_multi()",copy_log=T)
dict_roi <- func_dict_roi(paths_)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,key_group_)

# Prepare clinical data
data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=T,
                                   prefix=NULL,print_terminal=F)

# Calculate ROI-wise GAMM of FC
df_join<-join_fc_clin(data_fc$df_fc,data_clin$df_clin)
data_gamm<-iterate_gamm(df_join,data_fc$df_roi,list_mod_,calc_parallel=F,calc_identical=F,list_sex=list(c(1,2)))

# Calculate Group-wise GAMM of FC
df_join_grp<-join_fc_clin(data_fc$df_fc_grp,data_clin$df_clin)
data_gamm_grp<-iterate_gamm(df_join_grp,data_fc$df_grp,list_mod_,calc_parallel=F,calc_identical=T,list_sex=list(c(1,2)))

# Graphical output of ROI- and group-wise GAMM of FC
plot_gam_fc(paths_,df_gam=df_plot,df_gam_grp_sign=df_plot_grp,df_gam_grp_abs=NULL,atlas,
            list_mod,list_plot,list_type_p=list_type_p_,thr_p=thr_p_,waves=NULL,idx_var)