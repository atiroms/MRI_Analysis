#**************************************************
# Description =====================================
#**************************************************
# R script for commonly used parameters in MRI data analysis


#**************************************************
# gam_fc_diff() ===================================
#**************************************************
param_gam_fc_diff<-list(
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "key_group"="group_3",
  "list_wave"=c(1,2),
  "list_p"=list(list("type"="p","threshold"=0.001),
                list("type"="p_bh","threshold"=0.05)),
  "subset_subj"=list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                              list("key"="W1_rsfMRIexist","condition"="==1"),
                              list("key"="W1_Censor","condition"="<126")),
                     "2"=list(list("key"="W2_T1QC","condition"="==1"),
                              list("key"="W2_rsfMRIexist","condition"="==1"),
                              list("key"="W2_Censor","condition"="<126"))),
  "list_covar_tanner"=list("tanner"=NULL,
                           "age"   =list("1"="W1_Age_at_MRI", "2"="W2_Age_at_MRI", "label"="Age"),
                           "sex"   =list("1"="Sex",           "2"="Sex",           "label"="Sex")),
  #"list_tanner"=list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max",
  #                                  "label"="Tanner stage (max)","dtype"="factor"),
  #                   "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full",
  #                                  "label"="Tanner stage (full)","dtype"="factor"),
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
  #                                  "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)","dtype"="factor"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
  #                                  "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)","dtype"="factor")),
  "list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                    "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor")),
  "list_mod_tanner"=list("l" = "value ~ ses1_age + ses1_tanner + ses2_age + ses2_tanner"),
                         #"li"= "value ~ age * tanner"),
  "list_term_tanner"=list("a1"=list("title"="1st wave Age effect","var_exp"="ses1_age"),
                          "t1"=list("title"="1st wave Tanner effect","var_exp"="ses1_tanner"),
                          "a2"=list("title"="2nd wave Age effect","var_exp"="ses2_age"),
                          "t2"=list("title"="2nd wave Tanner effect","var_exp"="ses2_tanner")),
                          #"a"=list("title"="Age effect","var_exp"="age"),
                          #"s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                          #"t"=list("title"="Tanner effect","var_exp"="tanner"),
                          #"at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          #"s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  "list_hormone"=NULL,
  "list_mod_hormone"=list("l" = "value ~ ses1_age + ses1_hormone",
                          "li"= "value ~ age * hormone"),
  "list_term_hormone"=list("a"=list("title"="Age effect","var_exp"="age"),
                           "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                           "h"=list("title"="Hormone effect","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)")),
  "param_nbs"=list("list_mod"="l",
                   "list_term"=c("t1","t2"),
                   #"list_term"="t2",
                   "p_cdt_threshold"=0.001,
                   "p_perm_threshold"=0.05,
                   "n_perm"=60)
)


#**************************************************
# gamm_fc() =======================================
#**************************************************
param_gamm_fc<-list(
  "abs_nfc"=F, # absolute value for negative functional connectivity
  "key_group"="group_3",
  "list_wave"=c(1,2),
  "list_p"=list(list("type"="p","threshold"=0.001),
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
  #                   "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
  #                                  "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
  #                                  "label"="Tanner stage (gonadal)"),
  #                   "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
  #                                  "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
  #                                  "label"="Tanner stage (adrenal)")),
  "list_tanner"=list("gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                    "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                    "label"="Tanner stage (gonadal)","dtype"="factor")),
  "list_mod_tanner"=list("l" = "value ~ age + tanner + s(ID_pnTTC,bs='re')",
                         "li"= "value ~ age * tanner + s(ID_pnTTC,bs='re')"),
  #"a"= "value ~ age + s(tanner,k=3) + s(ID_pnTTC,bs='re')"),
  #"a"= "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')",
  #"q"="value ~ poly(age,2) + poly(tanner,2) + s(ID_pnTTC,bs='re')")
  "list_plot_tanner"=list("a"=list("title"="Age effect","var_exp"="age"),
                          "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                          "t"=list("title"="Tanner effect","var_exp"="tanner"),
                          "at"=list("title"="Age by Tanner interaction","var_exp"="age:tanner"),
                          "s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)")),
  "list_covar_hormone"=list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                            "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                            "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex")),
  #"list_hormone"=list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
  #                    "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
  #                    "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
  #                    "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S")),
  "list_hormone"=NULL,
  "list_mod_hormone"=list("l" = "value ~ age + hormone + s(ID_pnTTC,bs='re')",
                          "li"= "value ~ age * hormone + s(ID_pnTTC,bs='re')"),
  #"a"= "value ~ age + s(hormone,k=3) + s(ID_pnTTC,bs='re')"),
  #"a"= "value ~ s(age,k=3) + s(hormone,k=3) + s(ID_pnTTC,bs='re')",
  #"q"="value ~ poly(age,2) + poly(hormone,2) + s(ID_pnTTC,bs='re')")
  "list_plot_hormone"=list("a"=list("title"="Age effect","var_exp"="age"),
                           "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                           "h"=list("title"="Hormone effect","var_exp"="hormone"),
                           "ah"=list("title"="Age by Hormone interaction","var_exp"="age:hormone"),
                           "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)"))
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
# ca_fc_cs_multi() ================================ 
#**************************************************
ca_fc_subset_subj <- list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
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
                                     list("key"="W2_Censor","condition"="<126")))
ca_fc_list_sex <- list("male"="==1","female"="==2","all"=" %in% c(1,2)")
ca_fc_list_wave_mri<-list("m1"="1","m2"="2","m2-1"="2-1")
ca_fc_list_wave_clin<-list("c1"="1","c2"="2")
ca_fc_list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                              "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                              "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
ca_fc_list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                        "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                        "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                       "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                       "label"="Tanner stage (gonadal)"),
                        "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                       "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                       "label"="Tanner stage (adrenal)"))
ca_fc_list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                               "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                               "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
ca_fc_list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                         "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                         "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                         "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))


#**************************************************
# gam_fc_cs_multi() ===============================
#**************************************************
gam_fc_subset_subj <- list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
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
                                      list("key"="W2_Censor","condition"="<126")))
gam_fc_list_waves<-list("c1m1"  =list("wave_clin"="1","wave_mri"="1"),
                        "c1m2"  =list("wave_clin"="1","wave_mri"="2"),
                        "c1m2-1"=list("wave_clin"="1","wave_mri"="2-1"),
                        "c2m1"  =list("wave_clin"="2","wave_mri"="1"),
                        "c2m2"  =list("wave_clin"="2","wave_mri"="2"),
                        "c2m2-1"=list("wave_clin"="2","wave_mri"="2-1"))
gam_fc_list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                               "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                               "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
gam_fc_list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                         "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                         "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                        "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                        "label"="Tanner stage (gonadal)"),
                         "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                        "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                        "label"="Tanner stage (adrenal)"))
gam_fc_list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                                "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                                "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
gam_fc_list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                          "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                          "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                          "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
gam_fc_list_mod_tanner <- list("l"= "value ~ age + tanner",
                               "a"= "value ~ s(age,k=3) + s(tanner,k=3)")
                               #"q"="value ~ poly(age,2) + poly(tanner,2)")
gam_fc_list_mod_hormone <- list("l"= "value ~ age + hormone",
                                "a"= "value ~ s(age,k=3) + s(hormone,k=3)")
                                #"q"="value ~ poly(age,2) + poly(hormone,2)")
gam_fc_list_plot_tanner <-list("t"=list("title"="Tanner effect","var_exp"="tanner"),
                               "s(t)"=list("title"="Tanner effect","var_exp"="s(tanner)"))
gam_fc_list_plot_hormone <-list("h"=list("title"="Hormone effect","var_exp"="hormone"),
                                "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)"))


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

list_dim_ca<-c(5,10,20,40)
