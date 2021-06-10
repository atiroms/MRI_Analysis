# Test ordered factor


path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
path_exp_full<-NULL
#path_exp_full<-"/media/atiroms/SSD_03/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"421_fc_aroma"
#dir_out<-"423.3_fc_gam_diff_aroma_test1" 
dir_out<-"424_fc_gamm_aroma_test29" 
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
#**************************************************
# gamm_fc() =======================================
#**************************************************
param_gamm_fc<-list(
  "tfnbs"=T,
  "param_tfnbs"=list("e"=0.4,"h"=3.0,"n_thresh_h"=100),
  
  # Parameters for FC normalization
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=T, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  
  # Parameters for clinical data subsetting
  "force_long"=T, # use longitudinal data only
  #"omit_decreasing"="tanner", # omit subjects with longitudinally decreasing data of the variable
  "omit_decreasing"=NULL,
  #"group_tanner"=list("pre"=1,"early"=c(2,3),"late"=c(4,5)),
  "group_tanner"=NULL,
  
  "key_group"="group_3",
  "list_wave"=c(1,2),
  "list_sex"=list(1,2),
  "list_p"=list(list("type"="p","threshold"=0.001),
                list("type"="p","threshold"=0.005),
                list("type"="p","threshold"=0.01),
                list("type"="p_bh","threshold"=0.05)),
  "subset_subj"=list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                              list("key"="W1_rsfMRIexist","condition"="==1"),
                              list("key"="W1_Censor","condition"="<126")),
                     "2"=list(list("key"="W2_T1QC","condition"="==1"),
                              list("key"="W2_rsfMRIexist","condition"="==1"),
                              list("key"="W2_Censor","condition"="<126"))),
  "list_covar_tanner"=list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)","dtype"="factor"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="factor")),
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)","dtype"="ordered"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="ordered"),
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="ordered"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="ordered")),
  "list_mod_tanner"=list("l" = "value ~ age + tanner + (1|ID_pnTTC)"),
  "list_term_tanner"=list("a"=list("title"="Age effect","var_exp"="age"),
                          "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                          "t"=list("title"="Tanner effect","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)"),
                          "tl"=list("title"="Tanner effect","var_exp"="tanner.L"),
                          "tp-e"=list("title"="Tanner pre-early difference","var_exp"="tannerearly"),
                          "tp-l"=list("title"="Tanner pre-late difference","var_exp"="tannerlate")),
  "list_covar_hormone"=list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),

  "list_hormone"=NULL,
  #"list_mod_hormone"=list("l" = "value ~ age + hormone + (1|ID_pnTTC)",
  #                        "li"= "value ~ age * hormone + (1|ID_pnTTC)"),
  "list_mod_hormone"=list("l" = "value ~ age + hormone + (1|ID_pnTTC)"),

  "list_term_hormone"=list("a"=list("title"="Age effect","var_exp"="age"),
                           "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                           "h"=list("title"="Hormone effect","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)")),
  "param_nbs"=list(#"list_mod"=c("l","li"),
    "list_mod"="l",
    "list_term"=list(list("term_perm"="t","term_detect"=c("t","at","tl")),
                     list("term_perm"="h","term_detect"=c("h","ah"))),
    #"list_term"=list(list("term_perm"="t","term_detect"=c("t","tp-e","tp-l"))),
    #"p_cdt_threshold"=0.001,
    "p_cdt_threshold"=c(0.001,0.005,0.01),
    "p_perm_threshold"=0.05,
    #"n_perm"=1000),
    #"n_perm"=100),
  #"n_perm"=20),
  #"n_perm"=10),
    "n_perm"=3),
  "param_ancova_pred"=list("t"=data.frame(term=c("(Intercept)","tanner2","tanner3","tanner4","tanner5"),
                                          level=c(1,2,3,4,5)),
                           "at"=data.frame(term=c("age","age:tanner2","age:tanner3","age:tanner4","age:tanner5"),
                                           level=c(1,2,3,4,5)))
)
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
#test_mod=T

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
# Calculate model
data_gamm<-func_calc_gamm(paths,df_clin,df_fc,df_fc_grp,data_fc,calc_parallel,test_mod,
                          atlas,param,param$list_sex,list_covar,list_mod,list_term,idx_var,label_wave)
df_gamm<-data_gamm$df_gamm; df_anova<-data_gamm$df_anova; df_gamm_grp<-data_gamm$df_gamm_grp; df_anova_grp<-data_gamm$df_anova_grp

####

# Calculate threshold-free network based statistics
data_tfnbs<-func_iterate_tfnbs(paths,df_gamm,df_anova,data_fc,plot_result=T,return_nbs=T,
                               atlas,param,list_mod,list_term,idx_var,label_wave)
# Permutation test
func_tfnbs_permutation(paths,data_fc,df_clin,data_tfnbs$df_tfnbs,calc_parallel,plot_result=T,
                       atlas,param,list_mod,list_term,idx_var,label_wave)





####
####
estimate_base<-df_gamm_factor[df_gamm_factor$sex==1 & df_gamm_factor$term=="(Intercept)","estimate"]
list_estimate_factor<-estimate_base
for (tanner in 2:4){
  estimate<-df_gamm_factor[df_gamm_factor$sex==1 & df_gamm_factor$term==paste("tanner",as.character(tanner),sep=""),"estimate"]
  list_estimate_factor<-c(list_estimate_factor,estimate_base+estimate)
}

intercept<-df_gamm_ordered[df_gamm_ordered$sex==1 & df_gamm_ordered$term=="(Intercept)","estimate"]
term1<-df_gamm_ordered[df_gamm_ordered$sex==1 & df_gamm_ordered$term=="tanner.L","estimate"]
term2<-df_gamm_ordered[df_gamm_ordered$sex==1 & df_gamm_ordered$term=="tanner.Q","estimate"]
term3<-df_gamm_ordered[df_gamm_ordered$sex==1 & df_gamm_ordered$term=="tanner.C","estimate"]
list_estimate_ordered<-NULL
for (tanner in 1:4){
  contr_tanner<-contr.poly(4)[tanner,]
  estimate<-intercept+term1*contr_tanner[".L"]+term2*contr_tanner[".Q"]+term3*contr_tanner[".C"]
  list_estimate_ordered<-c(list_estimate_ordered,estimate)
}


####
####
df_join_copy<-df_join
df_join$age<-(df_join$age-mean(df_join$age))/sd(df_join$age)

mod<-lmer(as.formula("value ~ age + tanner + (1|ID_pnTTC)"),data=df_join)
