source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')
#source('D:/atiro/GitHub/MRI_Analysis/analyze/connection.R')

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

####
paths=paths_
key_group<-param$key_group
include_diff=T
abs_nfc=param$abs_nfc

#

data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

#gam_fc_diff_core(paths_,data_fc,atlas,param,list(1,2),
#                 list_covar,param$list_mod_tanner,param$list_term_tanner,idx_tanner,
#                 calc_parallel=T,test_mod=F)
#                 #calc_parallel=F,test_mod=F)

####

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
                              "sort","as.formula","as.numeric.factor",
                              "lm","summary","anova","AIC"),
              envir=environment())

#x<-iterate_ancova(clust,df_join=df_join_diff,df_edge=data_fc$df_edge,progressbar=T,test_mod=F)

# Calculate model
data_gamm<-iterate_ancova(clust,df_join_diff,data_fc$df_edge,progressbar=F,test_mod=test_mod)
data_gamm_grp<-iterate_ancova(clust,df_join_grp_diff,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
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

df_join_diff_subset<-df_join_diff[df_join_diff$from=="ho112_00001" & df_join_diff$to=="ho112_00002",]
formula_model<-list_mod$l

x<-aov(as.formula(formula_model),df_join_diff_subset)
View(x)

library(stats)
library(Rcpp)
library(RcppEigen)

#1 >> 6s
t_start<-Sys.time()
x<-summary(aov(as.formula(formula_model),df_join_diff_subset))[[1]]
print(Sys.time()-t_start) #0.0119
View(x)

#2 >> 6s
t_start<-Sys.time()
x1<-lm(as.formula(formula_model),df_join_diff_subset)
x2<-anova(x1)
x3<-summary(x1)
print(Sys.time()-t_start) #0.018
View(x1)
View(x2)
View(x3)



####

ancova_core<-function(df_src,list_mod_in=NULL,list_sex_in=NULL,
                      calc_parallel_in=NULL,test_mod_in=NULL){
  if(!is.null(list_mod_in)){list_mod<-list_mod_in}
  if(!is.null(list_sex_in)){list_sex<-list_sex_in}
  if(!is.null(calc_parallel_in)){calc_parallel<-calc_parallel_in}
  if(!is.null(test_mod_in)){test_mod<-test_mod_in}
  
  df_gamm<-df_anova<-df_aic<-data.frame()
  list_gamm_output<-NULL
  for (idx_mod in names(list_mod)){
    for (idx_sex in list_sex){
      label_sex<-paste(idx_sex,collapse="_")
      df_src_sex<-df_src[df_src$sex %in% idx_sex,]
      df_src_sex$value<-as.numeric(df_src_sex$value)
      mod<-try(lm(as.formula(list_mod[[idx_mod]]),data=df_src_sex), silent=F)
      
      #mod<-try(gam(as.formula(list_mod[[idx_mod]]),data=df_src_sex,method="REML",control=list(nthreads=1)), silent=F)
      if (class(mod)[1]!="try-error"){
        if (test_mod){
          list_gamm<-list("mod"=mod)
          names(list_gamm)<-paste("mod-",idx_mod,"_sex-",label_sex,sep="")
          list_gamm_output<-c(list_gamm_output,list_gamm)
        }
        p_table<-summary(mod)$coefficients
        df_gamm_add<-data.frame(term=rownames(p_table),estimate=p_table[,'Estimate'],
                                se=p_table[,'Std. Error'],F=NA,t=p_table[,'t value'],
                                p=p_table[,'Pr(>|t|)'])
        df_gamm<-rbind(df_gamm,
                       cbind(sex=label_sex,model=idx_mod,df_gamm_add))
        df_aic<-rbind(df_aic,
                      data.frame(sex=label_sex,model=idx_mod,aic=AIC(mod),aic_best=0))
        p_table_anova<-anova(mod)
        p_table_anova<-p_table_anova[rownames(p_table_anova)!="Residuals",c("Df","F value","Pr(>F)")]
        colnames(p_table_anova)<-c("df","F","p")
        df_anova<-rbind(df_anova,
                        cbind(sex=label_sex,model=idx_mod,term=rownames(p_table_anova),
                              p_table_anova))
      }
    } # Finished looping over sex
  }# Finished looping over model
  
  # Compare AICs of GAMM models
  df_aic_compare<-data.frame()
  for (idx_sex in list_sex){
    label_sex<-paste(idx_sex,collapse="_")
    df_aic_sex<-df_aic[df_aic$sex==label_sex,]
    df_aic_sex[which(df_aic_sex$aic==min(df_aic_sex$aic)),'aic_best']<-1
    df_aic_compare<-rbind(df_aic_compare,df_aic_sex)
  }
  
  # Prepare output dataframe
  if ("from" %in% colnames(df_src)){
    df_id<-df_src[1,c("from","to")]
    rownames(df_id)<-NULL
    df_gamm<-cbind(df_id,df_gamm)
    df_aic_compare<-cbind(df_id,df_aic_compare)
    df_anova<-cbind(df_id,df_anova)
  }
  
  return(list("df_gamm"=df_gamm,"df_aic"=df_aic_compare,"df_anova"=df_anova,"mod"=list_gamm_output))
}
