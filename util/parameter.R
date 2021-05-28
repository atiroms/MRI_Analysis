#**************************************************
# Description =====================================
#**************************************************
# R script for commonly used parameters in MRI data analysis


#**************************************************
# gamm_fc() =======================================
#**************************************************
param_gamm_fc<-list(
  # Parameters for FC normalization
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=T, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  
  # Parameters for clinical data subsetting
  "force_long"=F, # use longitudinal data only
  #"omit_decreasing"="tanner", # omit subjects with longitudinally decreasing data of the variable
  "omit_decreasing"=NULL,
  
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
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)")),
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)","dtype"="factor"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)","dtype"="ordered"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="ordered"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="ordered"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="ordered")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor")),
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)","dtype"="factor"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor")),
  #"list_mod_tanner"=list("l" = "value ~ age + tanner + (1|ID_pnTTC)",
  #                       "li"= "value ~ age * tanner + (1|ID_pnTTC)"),
  #                       "a"= "value ~ age + s(tanner,k=3) + s(ID_pnTTC,bs='re')"),
  #"list_mod_tanner"=list("l" = "value ~ age + tanner + (1|ID_pnTTC)",
  #                       "li"= "value ~ age * tanner + (1|ID_pnTTC)"),
  "list_mod_tanner"=list("l" = "value ~ age + tanner + (1|ID_pnTTC)"),
  "list_term_tanner"=list("a"=list("title"="Age effect","var_exp"="age"),
                          "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                          "t"=list("title"="Tanner effect","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)"),
                          "tl"=list("title"="Tanner effect","var_exp"="tanner.L")),
  "list_covar_hormone"=list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  #"list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
  #                    "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
  #                    "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
  #                    "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_hormone"=NULL,
  #"list_mod_hormone"=list("l" = "value ~ age + hormone + (1|ID_pnTTC)",
  #                        "li"= "value ~ age * hormone + (1|ID_pnTTC)"),
  "list_mod_hormone"=list("l" = "value ~ age + hormone + (1|ID_pnTTC)"),
  #"a"= "value ~ age + s(hormone,k=3) + s(ID_pnTTC,bs='re')"),
  #"a"= "value ~ s(age,k=3) + s(hormone,k=3) + s(ID_pnTTC,bs='re')",
  #"q"="value ~ poly(age,2) + poly(hormone,2) + s(ID_pnTTC,bs='re')")
  "list_term_hormone"=list("a"=list("title"="Age effect","var_exp"="age"),
                           "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                           "h"=list("title"="Hormone effect","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)")),
  "param_nbs"=list(#"list_mod"=c("l","li"),
                   "list_mod"="l",
                   "list_term"=list(list("term_perm"="t","term_detect"=c("t","at","tl")),
                                    list("term_perm"="h","term_detect"=c("h","ah"))),
                   #"p_cdt_threshold"=0.001,
                   "p_cdt_threshold"=c(0.001,0.005,0.01),
                   "p_perm_threshold"=0.05,
                   #"n_perm"=1000),
                   #"n_perm"=100),
                   "n_perm"=20),
                   #"n_perm"=10),
                   #"n_perm"=3),
  "param_ancova_pred"=list("t"=data.frame(term=c("(Intercept)","tanner2","tanner3","tanner4","tanner5"),
                                          level=c(1,2,3,4,5)),
                           "at"=data.frame(term=c("age","age:tanner2","age:tanner3","age:tanner4","age:tanner5"),
                                           level=c(1,2,3,4,5)))
)


#**************************************************
# gamm_fc() mixing both sex =======================
#**************************************************
param_gamm_fc_mix<-list(
  # Parameters for FC normalization
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=F, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  
  # Parameters for clinical data subsetting
  "force_long"=F, # use longitudinal data only
  "omit_decreasing"=NULL,
  #"omit_decreasing"=c("tanner_m","tanner_f"), # omit subjects with longitudinally decreasing data of the variable 
  "fill_na_tanner"=T,
  
  "key_group"="group_3",
  "list_wave"=c(1,2),
  "list_sex"=list(c(1,2)),
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
  "list_covar_tanner"=list("tanner_m"=NULL,
                           "tanner_f"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("gonadal"=list("male"=list("1"="W1_Tanner_Male_Genitals","2"="W2_Tanner_Male_Genitals",
                                                "label"="Tanner stage (male gonadal)","dtype"="factor"),
                                    "female"=list("1"="W1_Tanner_Female_Breast","2"="W2_Tanner_Female_Breast",
                                                  "label"="Tanner stage (female gonadal)","dtype"="factor")),
                     "adrenal"=list("male"=list("1"="W1_Tanner_Male_Pubic_Hair","2"="W2_Tanner_Male_Pubic_Hair",
                                                "label"="Tanner stage (male gonadal)","dtype"="factor"),
                                    "female"=list("1"="W1_Tanner_Female_Pubic_Hair","2"="W2_Tanner_Female_Pubic_Hair",
                                                  "label"="Tanner stage (female gonadal)","dtype"="factor"))),
  #"list_mod_tanner"=list("l" = "value ~ age + tanner_m + tanner_f + (1|ID_pnTTC)",
  #                       "li" = "value ~ age + tanner_m + age:tanner_m + tanner_f + age:tanner_f + (1|ID_pnTTC)"),
  "list_mod_tanner"=list("l" = "value ~ age + tanner_m + tanner_f + (1|ID_pnTTC)"),
  "list_term_tanner"=list("a"=list("title"="Age effect","var_exp"="age"),
                          "t_m"=list("title"="Tanner effect (Male)","var_exp"="tanner_m"),
                          "t_f"=list("title"="Tanner effect (Female)","var_exp"="tanner_f"),
                          "at_m"=list("title"="Age by Tanner (Male) interaction","var_exp"="age:tanner_m"),
                          "at_f"=list("title"="Age by Tanner (Female) interaction","var_exp"="age:tanner_f")),
  "list_hormone"=NULL,
  "param_nbs"=list(#"list_mod"=c("l","li"),
                   "list_mod"="l",
                   "list_term"=list(list("term_perm"="t_m","term_detect"=c("t_m","at_m")),
                                    list("term_perm"="t_f","term_detect"=c("t_f","at_f")),
                                    list("term_perm"="h","term_detect"=c("h","ah"))),
                   #"p_cdt_threshold"=0.001,
                   "p_cdt_threshold"=c(0.001,0.005,0.01),
                   "p_perm_threshold"=0.05,
                   #"n_perm"=1000),
                   #"n_perm"=100),
                   #"n_perm"=10),
                   "n_perm"=3),
  "param_ancova_pred"=list("t_m"=data.frame(term=c("(Intercept)","tanner_m2","tanner_m3","tanner_m4","tanner_m5"),
                                          level=c(1,2,3,4,5)),
                           "at_m"=data.frame(term=c("age","age:tanner_m2","age:tanner_m3","age:tanner_m4","age:tanner_m5"),
                                           level=c(1,2,3,4,5)),
                           "t_f"=data.frame(term=c("(Intercept)","tanner_f2","tanner_f3","tanner_f4","tanner_f5"),
                                            level=c(1,2,3,4,5)),
                           "at_f"=data.frame(term=c("age","age:tanner_f2","age:tanner_f3","age:tanner_f4","age:tanner_f5"),
                                             level=c(1,2,3,4,5)))
)



#**************************************************
# ca_str_cs() =====================================
#**************************************************
param_ca_str_cs<-list(
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=T, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  "list_dim_ca"=c(10,20,40),
  "key_group"="group_3",
  "subset_subj" = list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
                                  list("key"="W1_rsfMRIexist","condition"="==1"),
                                  list("key"="W1_Censor","condition"="<126")),
                       "2"  =list(list("key"="W2_T1QC","condition"="==1"),
                                  list("key"="W2_rsfMRIexist","condition"="==1"),
                                  list("key"="W2_Censor","condition"="<126")),
                       "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
                                  list("key"="W1_rsfMRIexist","condition"="==1"),
                                  list("key"="W1_Censor","condition"="<126"),
                                  list("key"="W2_T1QC","condition"="==1"),
                                  list("key"="W2_rsfMRIexist","condition"="==1"),
                                  list("key"="W2_Censor","condition"="<126"))),
  "list_sex" = list("male"="==1","female"="==2","all"=" %in% c(1,2)"),
  "list_wave_mri"=c("1","2","2-1"),
  "list_wave_clin"=c("1","2"),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)","dtype"="factor"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_mod_tanner"=list("l" = "value ~ age + tanner","li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a"=list("title"="Age","var_exp"="age"),
                          "s(a)"=list("title"="Age","var_exp"="s(age)"),
                          "t"=list("title"="Tanner","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=NULL,
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                      "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                      "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                      "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_mod_hormone"=list("l" = "value ~ age + hormone","li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a"=list("title"="Age","var_exp"="age"),
                           "s(a)"=list("title"="Age","var_exp"="s(age)"),
                           "h"=list("title"="Hormone","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone","var_exp"="s(hormone)"))
)


#**************************************************
# gam_str() ====================================== 
#**************************************************
param_gam_str<-list(
  "key_group"="group_3",
  "key_global_covar"="eTIV",
  "list_type_measure"=list(list("measure"="volume","global"=T),
                           list("measure"="thickness","global"=F),
                           list("measure"="area","global"=F)),
  #"list_type_measure"=list(list("measure"="volume","global"=F),
  #                         list("measure"="thickness","global"=F),
  #                         list("measure"="area","global"=F)),
  "subset_subj" = list("1"  =list(list("key"="W1_T1QC","condition"="==1")),
                       "2"  =list(list("key"="W2_T1QC","condition"="==1")),
                       "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
                                  list("key"="W2_T1QC","condition"="==1"))),
  "list_sex" = list("male"="==1","female"="==2","all"=" %in% c(1,2)"),
  "list_wave"=list("c1m1"=list("clin"=1,"mri"=1),
                   "c1m2"=list("clin"=1,"mri"=2),
                   "c2m1"=list("clin"=2,"mri"=1),
                   "c2m2"=list("clin"=2,"mri"=2)),
  #"list_wave_mri"=c("1","2","2-1"),
  #"list_wave_clin"=c("1","2"),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)"),
                     "max_fct"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)","dtype"="factor"),
                     "full_fct"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
                     "gonadal_fct"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor"),
                     "adrenal_fct"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_mod_tanner"=list("l" = "value ~ age + tanner","li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a"=list("title"="Age","var_exp"="age"),
                          "s(a)"=list("title"="Age","var_exp"="s(age)"),
                          "t"=list("title"="Tanner","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=NULL,
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                      "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                      "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                      "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_mod_hormone"=list("l" = "value ~ age + hormone","li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a"=list("title"="Age","var_exp"="age"),
                           "s(a)"=list("title"="Age","var_exp"="s(age)"),
                           "h"=list("title"="Hormone","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone","var_exp"="s(hormone)"))
)

#**************************************************
# ca_fc_cs() ====================================== 
#**************************************************
param_ca_fc_cs<-list(
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=T, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  "list_dim_ca"=c(10,20,40),
  "key_group"="group_3",
  "subset_subj" = list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
                                 list("key"="W1_rsfMRIexist","condition"="==1"),
                                 list("key"="W1_Censor","condition"="<126")),
                      "2"  =list(list("key"="W2_T1QC","condition"="==1"),
                                 list("key"="W2_rsfMRIexist","condition"="==1"),
                                 list("key"="W2_Censor","condition"="<126")),
                      "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
                                 list("key"="W1_rsfMRIexist","condition"="==1"),
                                 list("key"="W1_Censor","condition"="<126"),
                                 list("key"="W2_T1QC","condition"="==1"),
                                 list("key"="W2_rsfMRIexist","condition"="==1"),
                                 list("key"="W2_Censor","condition"="<126"))),
  "list_sex" = list("male"="==1","female"="==2","all"=" %in% c(1,2)"),
  "list_wave_mri"=c("1","2","2-1"),
  "list_wave_clin"=c("1","2"),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)","dtype"="factor"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_mod_tanner"=list("l" = "value ~ age + tanner","li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a"=list("title"="Age","var_exp"="age"),
                          "s(a)"=list("title"="Age","var_exp"="s(age)"),
                          "t"=list("title"="Tanner","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=NULL,
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                      "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                      "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                      "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_mod_hormone"=list("l" = "value ~ age + hormone","li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a"=list("title"="Age","var_exp"="age"),
                           "s(a)"=list("title"="Age","var_exp"="s(age)"),
                           "h"=list("title"="Hormone","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone","var_exp"="s(hormone)"))
)


#**************************************************
# variance_fp() ===================================
#**************************************************
param_variance_fp<-list(
  "list_wave"=list("long"=list("clin"=c(1,2),"mri"=c(1,2)),
                   "c1m1"=list("clin"=1,"mri"=1),
                   "c1m2"=list("clin"=1,"mri"=2),
                   "c2m1"=list("clin"=2,"mri"=1),
                   "c2m2"=list("clin"=2,"mri"=2)),
  #"list_wave"=c(1,2),
  "list_covar"=list("tanner"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
                    #"age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                    "sex"=list("1"="Sex","2"="Sex","label"="Sex")),
  "subset_subj"=list("1"=list(list("key"="W1_T1QC","condition"="==1"),list("key"="W1_rsfMRIexist","condition"="==1"),list("key"="W1_Censor","condition"="<126")),
                     "2"=list(list("key"="W2_T1QC","condition"="==1"),list("key"="W2_rsfMRIexist","condition"="==1"),list("key"="W2_Censor","condition"="<126")))
  
)


#**************************************************
# gam_fc_cs() =====================================
#**************************************************
param_gam_fc_cs<-list(
  "abs_nfc"=T, # absolute value for negative functional connectivity
  "key_group"="group_3",
  #"list_wave"=list("c1m1"=list("clin"=1,"mri"=1),
  #                 "c1m2"=list("clin"=1,"mri"=2),
  #                 "c2m1"=list("clin"=2,"mri"=1),
  #                 "c2m2"=list("clin"=2,"mri"=2)),
  #"list_wave"=list("c1m2"=list("clin"=1,"mri"=2)),
  "list_wave"=list("c2m1"=list("clin"=2,"mri"=1)),
  "list_p"=list(list("type"="p","threshold"=0.001),
                list("type"="p","threshold"=0.005),
                list("type"="p","threshold"=0.01),
                list("type"="p_bh","threshold"=0.05)),
  "subset_subj"=list("1"=list(list("key"="W1_T1QC","condition"="==1"),list("key"="W1_rsfMRIexist","condition"="==1"),list("key"="W1_Censor","condition"="<126")),
                     "2"=list(list("key"="W2_T1QC","condition"="==1"),list("key"="W2_rsfMRIexist","condition"="==1"),list("key"="W2_Censor","condition"="<126"))),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)","dtype"="factor"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor")),
  #"list_tanner"=NULL,
  "list_mod_tanner"=list("l" = "value ~ age + tanner","li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a"=list("title"="Age","var_exp"="age"),
                          "s(a)"=list("title"="Age","var_exp"="s(age)"),
                          "t"=list("title"="Tanner","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=NULL,
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  #"list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
  #                    "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
  #                    "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
  #                    "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  #"list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
  #                    "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
  #                    "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  #"list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone")),
  "list_hormone"=NULL,
  "list_mod_hormone"=list("l" = "value ~ age + hormone","li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a"=list("title"="Age","var_exp"="age"),
                           "s(a)"=list("title"="Age","var_exp"="s(age)"),
                           "h"=list("title"="Hormone","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone","var_exp"="s(hormone)")),
  "param_nbs"=list("list_mod"=c("l","li"),
                   "list_term"=list(list("term_perm"="t","term_detect"=c("t","at")),
                                    list("term_perm"="h","term_detect"=c("h","ah"))),
                   #"p_cdt_threshold"=0.001,
                   "p_cdt_threshold"=c(0.001,0.005,0.01),
                   "p_perm_threshold"=0.05,
                   "n_perm"=1000),
                   #"n_perm"=100),
                   #"n_perm"=10),
                   #"n_perm"=3),
  "param_ancova_pred"=list("t"=data.frame(term=c("(Intercept)","tanner2","tanner3","tanner4","tanner5"),level=c(1,2,3,4,5)),
                           "at"=data.frame(term=c("age","age:tanner2","age:tanner3","age:tanner4","age:tanner5"),level=c(1,2,3,4,5)))
)

#**************************************************
# gam_fc_diff() ===================================
#**************************************************
param_gam_fc_diff<-list(
  # Parameters for FC normalization
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "std_fc"=T, # standardize z values with demeaning and division with sd
  "div_mean_fc"=F, # normalize z values with division with mean
  
  # Parameters for clinical data subsetting
  "omit_decreasing"="tanner", # omit subjects with longitudinally decreasing data of the variable
  #"omit_decreasing"=NULL,
  
  "key_group"="group_3",
  "list_wave"=c(1,2),
  "list_p"=list(list("type"="p","threshold"=0.001),
                list("type"="p","threshold"=0.005),
                list("type"="p","threshold"=0.01),
                list("type"="p_bh","threshold"=0.05)),
  "subset_subj"=list("1"=list(list("key"="W1_T1QC","condition"="==1"),list("key"="W1_rsfMRIexist","condition"="==1"),list("key"="W1_Censor","condition"="<126")),
                     "2"=list(list("key"="W2_T1QC","condition"="==1"),list("key"="W2_rsfMRIexist","condition"="==1"),list("key"="W2_Censor","condition"="<126"))),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  "list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max","label"="Tanner stage (max)","dtype"="factor"),
                     "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)","dtype"="factor"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)","dtype"="factor")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),"2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  #"list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),"2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor")),
  "list_mod_tanner"=list("cs" = "value ~ diff_age + ses1_tanner + ses2_tanner",
                         "diff" = "value ~ diff_age + diff_tanner"),
                         #"l" = "value ~ ses1_age + ses1_tanner + ses2_age + ses2_tanner",
                         #"li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a1"=list("title"="1st wave Age","var_exp"="ses1_age"),
                          "t1"=list("title"="1st wave Tanner","var_exp"="ses1_tanner"),
                          "a2"=list("title"="2nd wave Age","var_exp"="ses2_age"),
                          "t2"=list("title"="2nd wave Tanner","var_exp"="ses2_tanner"),
                          "ad"=list("title"="Age difference","var_exp"="diff_age"),
                          "td"=list("title"="Tanner difference","var_exp"="diff_tanner"),
                          "am"=list("title"="Age mean","var_exp"="mean_age"),
                          "tm"=list("title"="Tanner mean","var_exp"="mean_tanner")),
  "list_covar_hormone"=list("hormone"=NULL,
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                      "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                      "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                      "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  #"list_hormone"=NULL,
  "list_mod_hormone"=list("cs" = "value ~ diff_age + ses1_hormone + ses2_hormone",
                          "diff" = "value ~ diff_age + diff_hormone"),
                          #"l" = "value ~ ses1_age + ses1_hormone + ses2_age + ses2_hormone",
                          #"li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a1"=list("title"="1st wave Age","var_exp"="ses1_age"),
                           "h1"=list("title"="1st wave Hormone","var_exp"="ses1_hormone"),
                           "a2"=list("title"="2nd wave Age","var_exp"="ses2_age"),
                           "h2"=list("title"="2nd wave Hormone","var_exp"="ses2_hormone"),
                           "ad"=list("title"="Age difference","var_exp"="diff_age"),
                           "hd"=list("title"="Hormone difference","var_exp"="diff_hormone"),
                           "am"=list("title"="Age mean","var_exp"="mean_age"),
                           "hm"=list("title"="Hormone mean","var_exp"="mean_hormone")),
  "param_nbs"=list("list_mod"=c("cs","diff"),
                   "list_term"=list(list("term_perm"="t1","term_detect"="t1"),
                                    list("term_perm"="t2","term_detect"="t2"),
                                    list("term_perm"="td","term_detect"="td"),
                                    list("term_perm"="tm","term_detect"="tm"),
                                    list("term_perm"="h1","term_detect"="h1"),
                                    list("term_perm"="h2","term_detect"="h2"),
                                    list("term_perm"="hd","term_detect"="hd"),
                                    list("term_perm"="hm","term_detect"="hm")),
                   #"p_cdt_threshold"=0.001,
                   "p_cdt_threshold"=c(0.001,0.005,0.01),
                   "p_perm_threshold"=0.05,
                   #"n_perm"=1000),
                   "n_perm"=20),
                   #"n_perm"=10),
                   #"n_perm"=3),
  "param_ancova_pred"=list("t1"=data.frame(term=c("(Intercept)","ses1_tanner2","ses1_tanner3","ses1_tanner4","ses1_tanner5"),
                                           level=c(1,2,3,4,5)),
                           "t2"=data.frame(term=c("(Intercept)","ses2_tanner2","ses2_tanner3","ses2_tanner4","ses2_tanner5"),
                                           level=c(1,2,3,4,5)),
                           "td"=data.frame(term=c("(Intercept)","diff_tanner1","diff_tanner2","diff_tanner3","diff_tanner4"),
                                           level=c(0,1,2,3,4)))
)


#**************************************************
# std_clin() ======================================
#**************************************************
param_std_clin<-list(
  "list_wave"=c(1,2),
  "subset_subj"=list("1"=list(),"2"=list()),
  "list_tanner"=list("max" =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                     "full"=list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                     "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                    "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)"),
                     "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                    "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                    "label"="Tanner stage (adrenal)")),
  "list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                      "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                      "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                      "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_covar"=list("age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                    "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_mod"=list("l"="value ~ age + s(ID_pnTTC,bs='re')",
                  "q"="value ~ poly(age,2) + s(ID_pnTTC,bs='re')",
                  "c"="value ~ poly(age,3) + s(ID_pnTTC,bs='re')",
                  "a"="value ~ s(age,k=3) + s(ID_pnTTC,bs='re')"),
                  #"a"="value ~ s(age) + s(ID_pnTTC,bs='re')"),
  "spec_graph"=list("x_axis"="age",
                    "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                  "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                    "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                 "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))
)


#**************************************************
# sex_diff_fc() ===================================
#**************************************************

sex_diff_fc_subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                                         list("key"="W1_rsfMRIexist","condition"="==1"),
                                         list("key"="W1_Censor","condition"="<126")),
                                "2"=list(list("key"="W2_T1QC","condition"="==1"),
                                         list("key"="W2_rsfMRIexist","condition"="==1"),
                                         list("key"="W2_Censor","condition"="<126")),
                                "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
                                           list("key"="W1_rsfMRIexist","condition"="==1"),
                                           list("key"="W1_Censor","condition"="<126"),
                                           list("key"="W2_T1QC","condition"="==1"),
                                           list("key"="W2_rsfMRIexist","condition"="==1"),
                                           list("key"="W2_Censor","condition"="<126")))
sex_diff_fc_list_covar<-list("age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                             "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex"))
sex_diff_fc_list_mod_diff <- list("int"  = "value ~ sex*mean_age")
sex_diff_fc_list_mod_long <- list("lin"  = "value ~ sex + age + s(ID_pnTTC,bs='re')",
                                  "int"  = "value ~ sex*age + s(ID_pnTTC,bs='re')")
sex_diff_fc_list_mod_cs   <- list("lin"  = "value ~ sex + age",
                                  "int"  = "value ~ sex*age")
sex_diff_fc_list_plot <- list("s"     =list("title"="Sex effect","var_exp"="sex2"),
                              "sxm(a)"=list("title"="Sex by mean(age) interaction","var_exp"="sex2:mean_age"),
                              "sxa"   =list("title"="Sex by age interaction","var_exp"="sex2:age"))
sex_diff_fc_thr_p_cdt <- 0.001
sex_diff_fc_thr_p_perm <- 0.05
sex_diff_fc_n_perm <- 1000
#sex_diff_fc_n_perm <- 3


#**************************************************
# model_fp_multi() ================================
#**************************************************
model_fp_subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                                      list("key"="W1_rsfMRIexist","condition"="==1"),
                                      list("key"="W1_Censor","condition"="<126")),
                             "2"=list(list("key"="W2_T1QC","condition"="==1"),
                                      list("key"="W2_rsfMRIexist","condition"="==1"),
                                      list("key"="W2_Censor","condition"="<126")))
model_fp_list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                                 "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                                 "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
model_fp_list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                           "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                           "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                          "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                          "label"="Tanner stage (gonadal)"),
                           "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                          "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                          "label"="Tanner stage (adrenal)"))
model_fp_list_mod_tanner <- list("ld" ="value ~ diff_age + diff_tanner",
                                 "ldm"="value ~ diff_age + diff_tanner + mean_tanner",
                                 "ad" ="value ~ s(diff_age,k=3) + s(diff_tanner,k=3)",
                                 "adm"="value ~ s(diff_age,k=3) + s(mean_tanner,k=3) + s(diff_tanner,k=3)")
model_fp_list_graph_tanner <-list("d(a)"=list("title"="Age diff effect","x_axis"="diff_age",
                                              "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                            "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                              "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                           "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                                  "d(t)"=list("title"="Tanner diff effect","x_axis"="diff_tanner",
                                              "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                            "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                              "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                           "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))),
                                  "m(t)"=list("title"="Tanner mean effect","x_axis"="mean_tanner",
                                              "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                            "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                              "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                           "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))))
model_fp_list_strat_tanner <-list("5by5"=list("1"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                              "2"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)),
                                  "3by3"=list("1"=list("12"=c(1,2),"3"=3,"45"=c(4,5)),
                                              "2"=list("12"=c(1,2),"3"=3,"45"=c(4,5))),
                                  "2by2"=list("1"=list("12"=c(1,2),"345"=c(3,4,5)),
                                              "2"=list("123"=c(1,2,3),"45"=c(4,5))))
model_fp_list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                                  "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                                  "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
model_fp_list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                            "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                            "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                            "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
model_fp_list_mod_hormone <- list("ld" ="value ~ diff_age + diff_hormone",
                                  "ldm"="value ~ diff_age + diff_hormone+ mean_hormone",
                                  "ad" ="value ~ s(diff_age,k=3) + s(diff_hormone,k=3)",
                                  "adm"="value ~ s(diff_age,k=3) + s(mean_hormone,k=3) + s(diff_hormone,k=3)")
model_fp_list_graph_hormone <-list("d(a)"=list("title"="Age diff effect","x_axis"="diff_age",
                                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                                   "d(h)"=list("title"="Hormone diff effect","x_axis"="diff_hormone",
                                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                                   "m(h)"=list("title"="Hormone mean effect","x_axis"="mean_hormone",
                                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))


# OBSOLETE ====

##**************************************************
## gam_fc_cs_multi() ===============================
##**************************************************
#gam_fc_subset_subj <- list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
#                                      list("key"="W1_rsfMRIexist","condition"="==1"),
#                                      list("key"="W1_Censor","condition"="<126")),
#                           "2"  =list(list("key"="W2_T1QC","condition"="==1"),
#                                      list("key"="W2_rsfMRIexist","condition"="==1"),
#                                      list("key"="W2_Censor","condition"="<126")),
#                           "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
#                                      list("key"="W1_rsfMRIexist","condition"="==1"),
#                                      list("key"="W1_Censor","condition"="<126"),
#                                      list("key"="W2_T1QC","condition"="==1"),
#                                      list("key"="W2_rsfMRIexist","condition"="==1"),
#                                      list("key"="W2_Censor","condition"="<126")))
#gam_fc_list_waves<-list("c1m1"  =list("wave_clin"="1","wave_mri"="1"),
#                        "c1m2"  =list("wave_clin"="1","wave_mri"="2"),
#                        "c1m2-1"=list("wave_clin"="1","wave_mri"="2-1"),
#                        "c2m1"  =list("wave_clin"="2","wave_mri"="1"),
#                        "c2m2"  =list("wave_clin"="2","wave_mri"="2"),
#                        "c2m2-1"=list("wave_clin"="2","wave_mri"="2-1"))
#gam_fc_list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
#                               "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
#                               "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
#gam_fc_list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
#                         "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
#                         "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
#                                        "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
#                                        "label"="Tanner stage (gonadal)"),
#                         "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
#                                        "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
#                                        "label"="Tanner stage (adrenal)"))
#gam_fc_list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
#                                "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                                "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
#gam_fc_list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
#                          "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
#                          "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
#                          "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
#gam_fc_list_mod_tanner <- list("l"= "value ~ age + tanner",
#                               "a"= "value ~ s(age,k=3) + s(tanner,k=3)")
#                               #"q"="value ~ poly(age,2) + poly(tanner,2)")
#gam_fc_list_mod_hormone <- list("l"= "value ~ age + hormone",
#                                "a"= "value ~ s(age,k=3) + s(hormone,k=3)")
#                                #"q"="value ~ poly(age,2) + poly(hormone,2)")
#gam_fc_list_plot_tanner <-list("t"=list("title"="Tanner effect","var_exp"="tanner"),
#                               "s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)"))
#gam_fc_list_plot_hormone <-list("h"=list("title"="Hormone effect","var_exp"="hormone"),
#                                "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)"))


#**************************************************
# plot_clin(), plot_pair(), plot_long() ===========
#**************************************************

#list_wave <- c(1,2)
#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))
#subset_subj <- list("1"=list(),
#                    "2"=list())
#list_mod <- list("l"="value ~ age + s(ID_pnTTC,bs='re')",
#                 "a"="value ~ s(age,k=3) + s(ID_pnTTC,bs='re')",
#                 "q"="value ~ poly(age,2) + s(ID_pnTTC,bs='re')")
#
#list_tanner<-list("max" =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
#                  "full"=list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
#                  "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
#                                 "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
#                                 "label"="Tanner stage (gonadal)"),
#                  "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
#                                 "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
#                                 "label"="Tanner stage (adrenal)"))
#list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
#                        "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                        "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
#spec_graph_tanner<-list("title"="Tanner vs Age","x_axis"="age",
#                        "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                      "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                        "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                     "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))
#
#list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
#                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
#                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
#                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
#list_covar_hormone<-list("hormone"=list("1"="W1_Hormone",   "2"="W2_Hormone",   "label"="Hormone"),
#                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
#spec_graph_hormone<-list("title"="Hormone vs Age","x_axis"="age",
#                         "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                       "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                         "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                      "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))
#
#list_pair<-list(c("gonadal","testo"),c("gonadal","dheas"),c("adrenal","testo"),c("adrenal","dheas"),
#                c("max","testo"),c("max","dheas"))
#list_mod_pair <- list("l"="value ~ hormone + s(ID_pnTTC,bs='re')",
#                      "a"="value ~ s(hormone,k=3) + s(ID_pnTTC,bs='re')",
#                      "q"="value ~ poly(hormone,2) + s(ID_pnTTC,bs='re')")
#list_covar_pair<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
#                      "hormone"=list("1"="W1_Hormone",   "2"="W2_Hormone",   "label"="Hormone"),
#                      "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                      "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
#spec_graph_pair<-list("title"="Tanner vs Hormone","x_axis"="hormone",
#                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))


#**************************************************
# model_fp() ======================================
#**************************************************
#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))
#list_covar<-list("sdq_td"=list("1"="W1_SDQ_tdJ",     "2"="W2_SDQ_tdJ",      "label"="SDQ_td"),
#                 "age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                 "sex"=list("1"="Sex","2"="Sex","label"="Sex"))
#list_mod <- list("ld"= "value ~ diff_age + diff_sdq_td",
#                 "ldm"= "value ~ diff_age + diff_sdq_td + mean_sdq_td")
#list_graph <-list("adiff"=list("title"="Age diff effect","x_axis"="diff_age",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
#                  "sdiff"=list("title"="SDQ diff effect","x_axis"="diff_sdq_td",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
#                  "smean"=list("title"="SDQ mean effect","x_axis"="mean_sdq_td",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
#list_strat_tanner <-NULL
#list_covar<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
#                 "age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                 "sex"=list("1"="Sex","2"="Sex","label"="Sex"))
#list_mod <- list("lin_diff"="value ~ diff_age + diff_tanner",
#                 "lin_diff_mean"="value ~ diff_age + diff_tanner + mean_tanner",
#                 "add_diff"="value ~ s(diff_age,k=3) + s(diff_tanner,k=3)",
#                 "add_diff_mean"="value ~ s(diff_age,k=3) + s(mean_tanner,k=3) + s(diff_tanner,k=3)")
#list_graph <-list("adiff"=list("title"="Age diff effect","x_axis"="diff_age",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
#                  "tdiff"=list("title"="Tanner diff effect","x_axis"="diff_tanner",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
#                  "tmean"=list("title"="Tanner mean effect","x_axis"="mean_tanner",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
#list_strat_tanner <-list("5by5"=list("1"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
#                                     "2"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)),
#                         "3by3"=list("1"=list("12"=c(1,2),"3"=3,"45"=c(4,5)),
#                                     "2"=list("12"=c(1,2),"3"=3,"45"=c(4,5))),
#                         "2by2"=list("1"=list("12"=c(1,2),"345"=c(3,4,5)),
#                                     "2"=list("123"=c(1,2,3),"45"=c(4,5))))

#list_covar<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
#                 "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
#                 "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
#                 "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"),
#                 "age"  =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
#                 "sex"  =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
#list_mod <- list("lin_diff"="value ~ diff_age + diff_testo",
#                 "lin_diff_mean"="value ~ diff_age + diff_testo + mean_testo",
#                 "add_diff"="value ~ s(diff_age,k=3) + s(diff_testo,k=3)",
#                 "add_diff_mean"="value ~ s(diff_age,k=3) + s(mean_testo,k=3) + s(diff_testo,k=3)")
#list_graph <-list("diff"=list("title"="Testosterone diff effect",
#                              "x_axis"="diff_testo",
#                              "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                            "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                              "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                           "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
#                  "mean"=list("title"="Testosterone mean effect",
#                               "x_axis"="mean_testo",
#                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
#                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
#                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
#                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
#list_strat_tanner <-NULL

#**************************************************
# gamm_fc() =======================================
#**************************************************

#list_covar<-list("sdq_td"=list("1"="W1_SDQ_tdJ",     "2"="W2_SDQ_tdJ",      "label"="SDQ_td"),
#                 "age"  =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
#                 "sex"  =list("1"="Sex",            "2"="Sex",            "label"="Sex"))

#list_mod <- list("l"= "value ~ age + testo + s(ID_pnTTC,bs='re')")
#"a"= "value ~ s(age,k=3) + s(testo,k=3) + s(ID_pnTTC,bs='re')",
#"q"="value ~ poly(age,2) + poly(testo,2) + s(ID_pnTTC,bs='re')")

#list_mod <- list("l"= "value ~ age + sdq_td + s(ID_pnTTC,bs='re')")

#list_plot <-list(#"a"=list("title"="Age effect","var_exp"="age"),
#"sa"=list("title"="Age effect","var_exp"="s(age)"),
#"pa1"=list("title"="Age effect","var_exp"="poly(age, 2)1"),
#"pa2"=list("title"="Age effect","var_exp"="poly(age, 2)2"),
#"t"=list("title"="Testosterone effect","var_exp"="testo")
#"st"=list("title"="Testosterone effect","var_exp"="s(testo)")
#"st"=list("title"="Testosterone effect","var_exp"="s(testo)"),
#"pt1"=list("title"="Testosterone effect","var_exp"="poly(testo, 2)1"),
#"pt2"=list("title"="Testosterone effect","var_exp"="poly(testo, 2)2")
#)

#list_plot <-list("sdq"=list("title"="SDQ effect","var_exp"="sdq_td"))


#**************************************************
# obsolete parameters in connection.R =============
#**************************************************
#list_type_p=c("p","p_bh","seed_p_bh")
#list_type_p="p_bh"
#thr_p <- 0.05
#
#list_cost<-seq(0.15,0.40,0.01)
#absolute<-T
#threshold<-NA
#
#list_dim_ca<-c(5,10,20,40)
