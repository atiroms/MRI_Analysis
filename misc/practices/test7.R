# Test ordered factor


path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full<-NULL
#path_exp_full<-"/media/atiroms/SSD_03/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"421_fc_aroma"
#dir_out<-"423.3_fc_gam_diff_aroma_test1" 
dir_out<-"424_fc_gamm_aroma_test17" 
#dir_out<-"424.1_fc_gamm_mix_aroma_test5" 
#dir_out<-"423.2_fc_gam_cs_aroma_test4" 
#dir_out<-"424_fc_gamm_aroma_test2"

#list_atlas<-c("aal116","gordon333","ho112","power264",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")
#list_atlas<-"aal116"
list_atlas<-"ho112"
#list_atlas<-c("ho112","power264")
#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100x7","schaefer200x7","schaefer400x7",
#              "schaefer100x17","schaefer200x17","schaefer400x17",
#              "shen268")


#**************************************************
# Libraries =======================================
#**************************************************
library(easypackages)
#libraries(ggplot2,GGally,igraph,ggrepel,colorRamps,tidyverse,parallel,mgcv,car,plyr,dplyr,data.table,pbapply,stringr)
libraries("ggplot2","colorRamps","tidyverse","parallel","mgcv","dplyr","data.table","pbapply","stringr","lmerTest")


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/gta_function.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)



####

paths_=paths
list_atlas_=list_atlas
param=param_gamm_fc

####

print("Starting gamm_fc().")
nullobj<-func_createdirs(paths_,str_proc="gamm_fc()",copy_log=T,list_param=param)
memory.limit(1000000)

atlas<-"ho112"

print(paste("Preparing FC data: ",atlas,sep=""))
#data_fc<-prep_data_fc(paths_,atlas,param$key_group,abs_nfc=param$abs_nfc)
data_fc<-prep_data_fc2(paths_,atlas,param$key_group,list_wave=c("1","2"),include_grp=T,
                       abs_nfc=param$abs_nfc,std_fc=param$std_fc,div_mean_fc=param$div_mean_fc)
data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-"gonadal"

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

####

list_mod=param$list_mod_tanner
list_term=param$list_term_tanner
idx_var=idx_tanner
calc_parallel=T
test_mod=F

####

# Prepare clinical data and demean
df_clin<-func_clinical_data_long(paths,param$list_wave,param$subset_subj,list_covar,rem_na_clin=T,
                                 prefix=paste("var-",idx_var,sep=""),print_terminal=F)$df_clin
# Select subjects with longitudinal data
if (param$force_long){
  list_id_subj<-df_clin[df_clin$wave==param$list_wave[1],'ID_pnTTC']
  for (wave in param$list_wave[-1]){
    list_id_subj<-sort(intersect(list_id_subj,df_clin[df_clin$wave==wave,'ID_pnTTC']))
  }
  df_clin<-df_clin[df_clin$ID_pnTTC %in% list_id_subj,]
}

####

df_clin<-func_std_clin(df_clin,separate_sex=T)$df_clin
fwrite(df_clin,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_clin.csv",sep="")),row.names=F)

# Prepare FC data
df_fc<-data_fc$df_fc; df_fc_grp<-data_fc$df_fc_grp
fwrite(df_fc,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_fc.csv",sep="")),row.names=F)
fwrite(df_fc_grp,file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_src_fc_grp.csv",sep="")),row.names=F)

label_wave<-"long"
# Calculate model

####

list_sex<-param$list_sex

####

print("Calculating GAMM/ANOVA")
# Join FC and clinical data
df_join<-join_fc_clin(df_fc,df_clin)
df_join_grp<-join_fc_clin(df_fc_grp,df_clin)

# Prepare parallelization cluster
if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                              "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC"),
              envir=environment())

# Calculate model
data_gamm<-iterate_gamm4(clust,df_join,data_fc$df_edge,progressbar=F,test_mod=test_mod)
data_gamm_grp<-iterate_gamm4(clust,df_join_grp,data_fc$df_edge_grp,progressbar=F,test_mod=test_mod)
stopCluster(clust)


# Add multiple comparison-corrected p values
df_gamm<-as.data.frame(add_mltcmp(data_gamm$df_gamm,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
df_anova<-as.data.frame(add_mltcmp(data_gamm$df_anova,data_fc$df_roi,list_mod,list_term,calc_seed_level=F))
df_gamm_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_gamm,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))
df_anova_grp<-as.data.frame(add_mltcmp(data_gamm_grp$df_anova,data_fc$df_grp,list_mod,list_term,calc_seed_level=F))



####
####
df_join_copy<-df_join
df_join$age<-(df_join$age-mean(df_join$age))/sd(df_join$age)

mod<-lmer(as.formula("value ~ age + tanner + (1|ID_pnTTC)"),data=df_join)
