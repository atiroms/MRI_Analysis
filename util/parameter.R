#**************************************************
# Description =====================================
#**************************************************
# R script for commonly used parameters in MRI data analysis


#**************************************************
# gamm_fc_multi() =================================
#**************************************************

subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                             list("key"="W1_rsfMRIexist","condition"="==1"),
                             list("key"="W1_Censor","condition"="<126")),
                    "2"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))

list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                        "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                        "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                  "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                  "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                 "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                 "label"="Tanner stage (gonadal)"),
                  "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                 "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                 "label"="Tanner stage (adrenal)"))
list_mod_tanner <- list("l"= "value ~ age + tanner + s(ID_pnTTC,bs='re')",
                        "a"= "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')",
                        "q"="value ~ poly(age,2) + poly(tanner,2) + s(ID_pnTTC,bs='re')")
list_plot_tanner <- list("a"=list("title"="Age effect","var_exp"="age"),
                         "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                         "t"=list("title"="Tanner stage effect","var_exp"="tanner"),
                         "s(t)"=list("title"="Tanner stage effect","var_exp"="s(tanner)"))
list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
list_mod_hormone <- list("l"= "value ~ age + hormone + s(ID_pnTTC,bs='re')",
                         "a"= "value ~ s(age,k=3) + s(hormone,k=3) + s(ID_pnTTC,bs='re')",
                         "q"="value ~ poly(age,2) + poly(hormone,2) + s(ID_pnTTC,bs='re')")
list_plot_hormone <- list("a"=list("title"="Age effect","var_exp"="age"),
                          "s(a)"=list("title"="Age effect","var_exp"="s(age)"),
                          "h"=list("title"="Hormone effect","var_exp"="hormone"),
                          "s(h)"=list("title"="Hormone effect","var_exp"="s(hormone)"))


#**************************************************
# gamm_fc() =================================
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