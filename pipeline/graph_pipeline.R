#**************************************************
# Description =====================================
#**************************************************
# R script to sequentially analize graph data


#**************************************************
# Create path list ================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
dir_in<-dir_out<-NULL

func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
                    path_exp_=path_exp,
                    dir_in_=dir_in,
                    dir_out_=dir_out){
  path_root<-NA
  for(p in list_path_root){
    if(file.exists(p)){
      path_root<-p
    }
  }
  if(is.na(path_root)){
    print("Error: root path could not be found.")
  }
  path_script <- file.path(path_root,"GitHub/MRI_Analysis")
  path_common <- file.path(path_root,"DropBox/MRI_img/pnTTC/puberty/common")
  path_in     <- file.path(path_root,path_exp_,dir_in_)
  path_out    <- file.path(path_root,path_exp_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))
source(file.path(paths$script,"util/plot.R"))
source(file.path(paths$script,"analyze/long_graph.R"))


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/graph_GV"
list_id_dir_graph<-list("acompcor"=203,"acompcor_gsr"=233)
list_wave<-c(1,2)

list_metric_local=c("efficiency_local_bin")
list_metric_global=c("efficiency_bin","efficiency_local_bin")

subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                             list("key"="W1_rsfMRIexist","condition"="==1"),
                             list("key"="W1_Censor","condition"="<126")),
                    "2"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))

list_tanner<-list("max" =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                  "full"=list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                  "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                 "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                 "label"="Tanner stage (gonadal)"),
                  "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                 "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                 "label"="Tanner stage (adrenal)"))
list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                        "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                        "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
list_mod_tanner_diff <- list("lin_diff"="value ~ diff_age + diff_tanner",
                        "lin_diff_mean"="value ~ diff_age + diff_tanner + mean_tanner",
                        "add_diff"="value ~ s(diff_age,k=3) + s(diff_tanner,k=3)",
                        "add_diff_mean"="value ~ s(diff_age,k=3) + s(mean_tanner,k=3) + s(diff_tanner,k=3)")
list_graph_tanner_diff <-list("diff"=list("title"="Tanner diff effect",
                                     "x_axis"="diff_tanner",
                                     "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                   "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                     "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                  "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))),
                         "mean"=list("title"="Tanner mean effect",
                                     "x_axis"="mean_tanner",
                                     "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                   "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                     "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                  "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))))
list_mod_tanner <- list("lin"  ="value ~ age + tanner + s(ID_pnTTC,bs='re')",
                         "add"  ="value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')",
                         "quad" ="value ~ poly(age,2) + poly(tanner,2) + s(ID_pnTTC,bs='re')")
list_graph_tanner <-list("t"=list("title"="Tanner effect","x_axis"="tanner",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))

list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
list_covar_hormone<-list("hormone"=list("1"="W1_Hormone",   "2"="W2_Hormone",   "label"="Hormone"),
                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
list_mod_hormone_diff <- list("lin_diff"="value ~ diff_age + diff_hormone",
                         "lin_diff_mean"="value ~ diff_age + diff_hormone+ mean_hormone",
                         "add_diff"="value ~ s(diff_age,k=3) + s(diff_hormone,k=3)",
                         "add_diff_mean"="value ~ s(diff_age,k=3) + s(mean_hormone,k=3) + s(diff_hormone,k=3)")
list_graph_hormone_diff <-list("diff"=list("title"="Hormone diff effect",
                                      "x_axis"="diff_hormone",
                                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                          "mean"=list("title"="Hormone mean effect",
                                      "x_axis"="mean_hormone",
                                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
list_mod_hormone <- list("lin"  ="value ~ age + hormone + s(ID_pnTTC,bs='re')",
                         "add"  ="value ~ s(age,k=3) + s(hormone,k=3) + s(ID_pnTTC,bs='re')",
                         "quad" ="value ~ poly(age,2) + poly(hormone,2) + s(ID_pnTTC,bs='re')")
list_graph_hormone <-list("h"=list("title"="Hormone effect","x_axis"="hormone",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))


#**************************************************
# Iterate gamm_gta() over Tanner and homone data===
#**************************************************
gamm_gta_multi<-function(list_id_dir_=list_id_dir_graph,list_wave_=list_wave,subset_subj_=subset_subj,
                         list_metric_global_=list_metric_global,list_metric_local_=list_metric_local,
                         list_tanner_=list_tanner,list_covar_tanner_=list_covar_tanner,
                         list_mod_tanner_=list_mod_tanner,list_graph_tanner_=list_graph_tanner,
                         list_mod_tanner_diff_=list_mod_tanner_diff,list_graph_tanner_diff_=list_graph_tanner_diff,
                         list_hormone_=list_hormone,list_covar_hormone_=list_covar_hormone,
                         list_mod_hormone_=list_mod_hormone,list_graph_hormone_=list_graph_hormone,
                         list_mod_hormone_diff_=list_mod_hormone_diff,list_graph_hormone_diff_=list_graph_hormone_diff
                         ){
  
  print("Starting gamm_gta_multi()")
  for (suffix_dir in names(list_id_dir_)){
    id_dir_graph<-list_id_dir_[[suffix_dir]]
    dir_in<-paste(as.character(id_dir_graph),"graph",suffix_dir,sep='_')
    id_dir_gamm<-id_dir_graph+1
    for (idx_tanner in names(list_tanner_)){
      id_dir_gamm<-id_dir_gamm+0.1
      print(paste("Preproc: ",suffix_dir,", Tanner: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
      list_covar_tanner_[["tanner"]]<-list_tanner_[[idx_tanner]]
      dir_out<-paste(as.character(id_dir_gamm),"gamm_graph",suffix_dir,idx_tanner,sep='_')
      paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
      nullobj<-gamm_gta(paths_=paths,subset_subj_=subset_subj_,list_covar=list_covar_tanner_,
                        list_wave_=list_wave_,list_metric_global_=list_metric_global_,
                        list_metric_local_=list_metric_local_,list_mod_=list_mod_tanner_,
                        list_graph_=list_graph_tanner_,list_mod_diff_=list_mod_tanner_diff_,
                        list_graph_diff_=list_graph_tanner_diff)
    }
    for (idx_hormone in names(list_hormone_)){
      id_dir_gamm<-id_dir_gamm+0.1
      print(paste("Preproc: ",suffix_dir,", Hormone: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
      list_covar_hormone_[["hormone"]]<-list_hormone_[[idx_hormone]]
      dir_out<-paste(as.character(id_dir_gamm),"gamm_graph",suffix_dir,idx_hormone,sep='_')
      paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
      nullobj<-gamm_gta(paths_=paths,subset_subj_=subset_subj_,list_covar=list_covar_hormone_,
                        list_wave_=list_wave_,list_metric_global_=list_metric_global_,
                        list_metric_local_=list_metric_local_,list_mod_=list_mod_hormone_,
                        list_graph_=list_graph_hormone_,list_mod_diff_=list_mod_hormone_diff_,
                        list_graph_diff_=list_graph_hormone_diff)
    }
  }
  print("Finished gamm_gta_multi()")
}