paths_=paths
subset_subj_=model_fp_subset_subj
list_atlas_=list_atlas
list_wave_=list_wave
skip_ancova=T
list_covar_tanner_=model_fp_list_covar_tanner
list_tanner_=model_fp_list_tanner
list_mod_tanner_=model_fp_list_mod_tanner
list_graph_tanner_=model_fp_list_graph_tanner
list_strat_tanner_=model_fp_list_strat_tanner
list_covar_hormone_=model_fp_list_covar_hormone
list_hormone_=model_fp_list_hormone
list_mod_hormone_=model_fp_list_mod_hormone
list_graph_hormone_=model_fp_list_graph_hormone



idx_tanner<-names(list_tanner_)[[1]]



list_covar_=list_covar
list_mod_=list_mod_tanner_
list_graph_=list_graph_tanner_
suffix=paste("var-",idx_tanner,sep="")



atlas<-"aal116"
measure<-"fc"
model<-names(list_mod_)[[1]]
idx_sex<-1

term<-list_term[1]
