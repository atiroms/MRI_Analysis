source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')

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
data_fc<-prep_data_fc(paths_,atlas,param$key_group,include_diff=T,abs_nfc=param$abs_nfc)
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
#calc_parallel=F
calc_parallel=T
test_mod=F
#test_mod=T

####

#print("Calculated result already exists.")
#df_gamm<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))
#df_gamm_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))
#df_anova<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))
#df_anova_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))
df_clin_diff<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_clin_diff.csv",sep=""))))
df_fc_diff<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc_diff.csv",sep=""))))
df_fc_grp_diff<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_fc_grp_diff.csv",sep=""))))

####

# Permutation for NBS
# Prepare parallelization cluster
if (calc_parallel){
  clust<-makeCluster(floor(detectCores()*3/4))
}else{
  clust<-makeCluster(1)
}
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                              "sort","gam","as.formula","summary.gam",
                              "anova.gam","as.numeric.factor"),
              envir=environment())

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
    #df_join_grp_diff<-join_fc_clin(df_fc_grp_diff,df_clin_perm)
    
    # Calculate model
    data_gamm<-iterate_gamm3(clust,df_join_diff,data_fc$df_edge,progressbar=F,test_mod=test_mod)
    #data_gamm_grp<-iterate_gamm3(clust,df_join_grp_diff,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
    
    for (idx_mod in param$param_nbs$list_mod){
      for (idx_sex in list_sex){
        # Subset GAMM result dataframe for plotting
        df_gamm<-data_gamm$df_gamm
        df_gamm_subset<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp & df_gamm$sex==idx_sex,]
        if (nrow(df_gamm_subset)==0){
          df_anova<-data_gamm$df_anova
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
}
stopCluster(clust)
write.csv(df_max_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_perm_max_size.csv",sep="")),row.names = F)

# Summarize permutation result
df_threshold_size<-NULL
for (idx_mod in param$param_nbs$list_mod){
  for (idx_term in param$param_nbs$list_term){
    var_exp<-list_term[[idx_term]][["var_exp"]]
    for (idx_sex in list_sex){
      list_max_size<-df_max_size[df_max_size$model==idx_mod & df_max_size$term==var_exp
                                  & df_max_size$sex==idx_sex,"max_size"]
      list_max_size<-sort(list_max_size)
      thr_size_nbs<-list_max_size[ceiling(length(list_max_size)*(1-param$param_nbs$p_perm_threshold))]
      df_threshold_size<-rbind(df_threshold_size,data.frame(model=idx_mod,term=var_exp,sex=idx_sex,
                                                            max_size=max_size,thr_size=thr_size_nbs))
    }
  }
}
write.csv(df_max_size,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_perm_thr_size.csv",sep="")),row.names = F)
