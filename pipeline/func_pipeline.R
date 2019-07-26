#**************************************************
# Description =====================================
#**************************************************
# R script to sequentially analize functional data.


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
source(file.path(paths$script,"analyze_img/timeseries.R"))
source(file.path(paths$script,"analyze_img/connection.R"))
source(file.path(paths$script,"analyze_img/fingerprint.R"))


#**************************************************
# Parameters ======================================
#**************************************************

id_dir_ts<-101
suffix_dir<-"acompcor"
list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"
n_permutation<-1000
#n_permutation<-100

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))
list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))
list_mod <- list("lin_diff_a"=
                   "value ~ diff_age + sex",
                 "lin_diff_t"=
                   "value ~ sex + sex:diff_tanner",
                 "lin_diff_at"=
                   "value ~ diff_age + sex + sex:diff_tanner",
                 "add_diff_a"=
                   "value ~ s(diff_age,k=3) + sex",
                 "add_diff_t"=
                   "value ~ sex + s(diff_tanner,k=3,by=sex)",
                 "add_diff_at"=
                   "value ~ s(diff_age,k=3) + sex + s(diff_tanner,k=3,by=sex)"
                 )
list_graph <-list("a"=list("title"="Effect of age difference",
                           "x_axis"="diff_age",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),
                                                    "color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),
                                                      "color"="lightcoral","alpha"=1))),
                  "st"=list("title"="Effect of Tanner stage difference",
                            "x_axis"="diff_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "sat"=list("title"="Age difference-Tanner stage difference interaction",
                             "x_axis"="diff_age",
                             "smooth"=list("Male delta TS = -1"=list("fix"=list("sex"=1,"diff_tanner"=-1),
                                                                     "color"="Steelblue2","alpha"=0.4,"ribbon"=F),
                                           "Male delta TS = 1"=list("fix"=list("sex"=1,"diff_tanner"=1),
                                                                    "color"="steelblue2","alpha"=0.7,"ribbon"=F),
                                           "Male delta TS = 3"=list("fix"=list("sex"=1,"diff_tanner"=3),
                                                                    "color"="steelblue2","alpha"=1,"ribbon"=F),
                                           "Female delta TS = -1"=list("fix"=list("sex"=2,"diff_tanner"=-1),
                                                                       "color"="lightcoral","alpha"=0.4,"ribbon"=F),
                                           "Female delta TS = 1"=list("fix"=list("sex"=2,"diff_tanner"=1),
                                                                      "color"="lightcoral","alpha"=0.7,"ribbon"=F),
                                           "Female delta TS = 3"=list("fix"=list("sex"=2,"diff_tanner"=3),
                                                                      "color"="lightcoral","alpha"=1,"ribbon"=F)),
                             "point"=list("Male"=list("subset"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1),
                                          "Female"=list("subset"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1))))


#**************************************************
# Libraries =======================================
#**************************************************



#**************************************************
# Timeseries to GAMM of Fingerprints ==============
#**************************************************
pipe_func<-function(paths_=paths,
                    id_dir_ts_=id_dir_ts,
                    suffix_dir_=suffix_dir,
                    list_atlas_=list_atlas,
                    list_wave_=list_wave,
                    list_covar_=list_covar,
                    list_mod_=list_mod,list_graph_=list_graph,
                    subset_subj_=subset_subj,
                    n_permutation_=n_permutation){
  
  print('Starting pipe_func().')
  
  # Timeseries to functional connectivity
  id_dir_ts<-id_dir_ts_
  dir_in<-paste(as.character(id_dir_ts),"ts",suffix_dir_,sep='_')
  id_dir_fc<-id_dir_ts+1
  dir_out<-paste(as.character(id_dir_fc),"fc",suffix_dir_,sep='_')
  paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
  nullobj<-fc(paths_=paths,list_atlas_=list_atlas_)
  
  # Functional connetivity to fingerprint
  dir_in<-paste(as.character(id_dir_fc),"fc",suffix_dir_,sep='_')
  id_dir_fp<-id_dir_fc+1
  dir_out<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
  paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
  nullobj<-fp(paths_=paths,list_atlas_=list_atlas_,subset_subj_=subset_subj_)
  
  # Fingerprint to identification of fingerprints
  dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
  id_dir_idfp<-id_dir_fp+1
  dir_out<-paste(as.character(id_dir_idfp),"idfp",suffix_dir_,sep='_')
  paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
  nullobj<-identify_fp(paths_=paths,list_atlas_=list_atlas_,
                       list_wave_=list_wave_,
                       subset_subj_=subset_subj_,
                       n_permutation_=n_permutation_)
  
  # Fingerprint to GLM / ANCOVA of fingerprint difference
  dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
  id_dir_gamfp<-id_dir_fp+2
  dir_out<-paste(as.character(id_dir_gamfp),"fp_model",suffix_dir_,sep='_')
  paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
  nullobj<-glm_ancova_fp(paths_=paths,list_atlas_=list_atlas_,
                         list_wave_=list_wave_,list_covar_=list_covar_,
                         list_mod_=list_mod_,list_graph_=list_graph_,
                         subset_subj_=subset_subj_)
  
  print('Finished pipe_func().')
  
}
