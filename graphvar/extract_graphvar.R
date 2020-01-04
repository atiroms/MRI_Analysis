#**************************************************
# Description =====================================
#**************************************************
# R script to extract XGraphVar-calculated graph metrices


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/graph_GV"

list_id_dir<-list("acompcor"=201,
                  "aroma"=211,
                  "acompcor_gsr"=231,
                  "aroma_gsr"=241,
                  "acompcor"=301,
                  "aroma"=311,
                  "acompcor_gsr"=331,
                  "aroma_gsr"=341)

dir_in<-"201_fc_mat_acompcor"
dir_out<-"202_graph_acompcor"

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100","schaefer200","schaefer400","shen268")
list_atlas<-"power264"


#**************************************************
# Libraries =======================================
#**************************************************
library(R.matlab)
library(tidyr)
library(parallel)


#**************************************************
# Create path list ================================
#**************************************************
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


#**************************************************
# Data extraction =================================
#**************************************************
extract_gv<-function(paths_=paths,list_atlas_=list_atlas){
  
}