#**************************************************
# Description =====================================
#**************************************************
# R script to sequentially analize structural data.


#**************************************************
# Create path list ================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/str_FS"
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
source(file.path(paths$script,"analyze/long_structure.R"))


#**************************************************
# Parameters ======================================
#**************************************************

dir_in <-"01_extract"
id_dir_start<-02

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/str_FS"

file_input<-"fs_measure.csv"

list_wave <- c(1,2)

list_measure <-c("volume","thickness","area")
#list_measure <-"volume"

list_str_group<-c("cortex","subcortex","white matter","global","misc")
#list_str_group<-"subcortex"
#list_str_group<-c("global","misc")
#list_str_group<-c("cortex","subcortex","global")

list_hormone<-list("testo"=list("1"="W1_Testosterone",
                                "2"="W2_Testosterone",
                                "label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",
                                "2"="W2_Cortisol",
                                "label"="Cortisol"),
                   "dhea"=list("1"="W1_DHEA",
                               "2"="W2_DHEA",
                               "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",
                                "2"="W2_DHEAS",
                                "label"="DHEA-S"))

list_covar<-list("age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild","value"=1)))

list_mod <- list("lin"=
                   "value ~ age + hormone + s(ID_pnTTC,bs='re')",
                 "add"=
                   "value ~ s(age,k=3) + s(hormone,k=3) + s(ID_pnTTC,bs='re')",
                 "quad"=
                   "value ~ poly(age,2) + poly(hormone,2) + s(ID_pnTTC,bs='re')")

list_graph <-list("a"=list("title"="Age effect",
                           "x_axis"="age",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),
                                                    "color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),
                                                      "color"="lightcoral","alpha"=1))),
                  "h"=list("title"="Hormone effect",
                           "x_axis"="hormone",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),
                                                    "color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),
                                                      "color"="lightcoral","alpha"=1))))


#**************************************************
# Iterate over hormones ===========================
#**************************************************

pipe_str<-function(dir_in_=dir_in,id_dir_start_=id_dir_start,list_str_group_=list_str_group,
                   file_input_=file_input,
                   list_measure_=list_measure,list_hormone_=list_hormone,
                   list_wave_=list_wave,list_covar_=list_covar,
                   list_mod_=list_mod,list_graph_=list_graph,
                   subset_subj_=subset_subj
                   ){
  print("Starting pipe_str()")
  id_dir_cnt<-id_dir_start_
  
  for (hormone in names(list_hormone_)){
    print(paste("Calculating hormone:",list_hormone_[[hormone]][["label"]],sep=" "))
    id_dir_cnt<-id_dir_cnt+0.1
    dir_out<-paste(as.character(id_dir_cnt),"gamm",hormone,sep='_')
    paths<-func_path(dir_in_=dir_in_,dir_out_=dir_out)
    list_covar<-c(list_hormone_[hormone],list_covar_)
    list_mod<-list_mod_
    for (idx_mod in names(list_mod_)){
      list_mod[idx_mod]<-gsub("hormone",hormone,list_mod[idx_mod],fixed=T)
    }
    list_graph<-list_graph_
    list_graph[["h"]][["title"]]<-paste(list_hormone_[[hormone]][["label"]],"effect",sep=' ')
    list_graph[["h"]][["x_axis"]]<-hormone

    nullobj<-gamm_str(paths_=paths,subset_subj_=subset_subj_,list_covar_=list_covar,file_input_=file_input_,
                      list_wave_=list_wave_,list_measure_=list_measure_,list_str_group_=list_str_group_,
                      list_mod_=list_mod,list_graph_=list_graph,key_group_='group_3')
  }
  
  print("Finished pipe_str()")

}
