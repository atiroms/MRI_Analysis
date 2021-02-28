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
data_fc<-prep_data_fc(paths_,atlas,param$key_group,abs_nfc=param$abs_nfc)
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
df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                 prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
df_clin<-func_demean_clin(df_clin,thr_cont=6,separate_sex=T)$df_clin # thr_cont=4 to demean Tanner, =5 not to
df_join<-join_fc_clin(data_fc$df_fc,df_clin)
df_join_grp<-join_fc_clin(data_fc$df_fc_grp,df_clin)

####

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
data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
data_gamm_grp<-iterate_gamm4(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
stopCluster(clust)

# Add multiple comparison-corrected p values
df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))

df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))

# Save results
write.csv(df_gamm,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
write.csv(data_gamm$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
write.csv(df_anova,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep="")),row.names = F)
write.csv(df_gamm_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_aic,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
write.csv(df_anova_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep="")),row.names = F)


####
####

list_mod<-list("l"="value ~ age + tanner + (1|ID_pnTTC)")
# Prepare parallelization cluster
if (calc_parallel){
  clust<-makeCluster(floor(detectCores()*3/4))
}else{
  clust<-makeCluster(1)
}
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                              "as.formula","as.numeric.factor",
                              "lmer","summary","anova","AIC"),
              envir=environment())
# Calculate model
t_start<-Sys.time()
data_gamm<-iterate_glmm(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
print(Sys.time()-t_start) # Time difference of 12.04399 mins(MRO 4.0.2, single) Time difference of 3.16382 mins(R native 4.0.2, parallel)
stopCluster(clust)

####

list_mod<-list("l"="value ~ age + tanner + s(ID_pnTTC,bs='re')")
# Prepare parallelization cluster
if (calc_parallel){
  clust<-makeCluster(floor(detectCores()*3/4))
}else{
  clust<-makeCluster(1)
}
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod",
                              "as.formula","as.numeric.factor",
                              "gam","summary.gam","anova.gam"),
              envir=environment())
t_start<-Sys.time()
data_gamm<-iterate_gamm3(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
print(Sys.time()-t_start) # Time difference of 46.00904 mins(MRO 4.0.2, single) Time difference of 41.33795 mins (R native 4.0.2, parallel)
stopCluster(clust)

####

list_mod<-list("l"="value ~ age + tanner + (1|ID_pnTTC)")
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
t_start<-Sys.time()
data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
print(Sys.time()-t_start) # Time difference of 3.09792 mins (native R 4.0.2, parallel)
stopCluster(clust)

####

sessionInfo()
library(RhpcBLASctl)

get_num_cores()
get_num_procs()
blas_get_num_procs()
#blas_set_num_threads(threads)
omp_get_num_procs()
omp_get_max_threads()
#omp_set_num_threads(threads)

blas_set_num_threads(1)
omp_set_num_threads(1)

####

df_test<-df_join[df_join$from=="ho112_00001" & df_join$to=="ho112_00002",]
rownames(df_test)<-NULL
df_test_del<-df_test[-176,]

#library(lme4)
#library(languageR)
library(lmerTest)

t_start<-Sys.time()
mod1<-lmer(value~age+tanner+(1|ID_pnTTC),data=df_test)
x<-summary(mod1)
y<-anova(mod1)
print(Sys.time()-t_start)

t_start<-Sys.time()
mod3<-lmer(value~age+tanner+(1|ID_pnTTC),data=df_test_del)
x<-summary(mod3)
y<-anova(mod3)
print(Sys.time()-t_start)



t_start<-Sys.time()
mod2<-gam(value~age+tanner+s(ID_pnTTC,bs='re'),data=df_test)
x<-summary.gam(mod2)
y<-anova.gam(mod2)
print(Sys.time()-t_start)
summary(mod2)


####
