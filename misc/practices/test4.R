#source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')
source('D:/atiro/GitHub/MRI_Analysis/analyze/connection.R')

####

paths_=paths
list_atlas_=list_atlas
param<-param_gamm_fc

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
#gamm_fc_core(paths_,data_fc,atlas,param,list(1,2),list_covar,
#             param$list_mod_tanner,param$list_term_tanner,idx_tanner,
#             calc_parallel=T,test_mod=F)

####

paths=paths_
list_sex<-list(1,2)
list_mod<-param$list_mod_tanner
list_term<-param$list_term_tanner
idx_var<-idx_tanner
calc_parallel<-T
#calc_parallel<-F
test_mod<-F

####

df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                 prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
df_clin<-func_demean_clin(df_clin,thr_cont=4,separate_sex=T)$df_clin # thr_cont=4 to demean Tanner, =5 not to
fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_clin_diff.csv",sep="")),row.names=F)
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

####
t_start<-Sys.time()
data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
print(Sys.time()-t_start)
# 12.81098 mins in Ubuntu OpenBLAS terminal-defined single threading 


data_gamm_grp<-iterate_gamm4(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
stopCluster(clust)

