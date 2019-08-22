#**************************************************
# Description =====================================
#**************************************************
# R script to extract XCP-preprocessed functional correlation data for use in GraphVar


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"102_fc_acompcor"
dir_out<-"107_corrmat"

list_atlas<-c("aal116","glasser360","gordon333","power264",
              "schaefer100","schaefer200","schaefer400","shen268")


#**************************************************
# Libraries =======================================
#**************************************************
library(R.matlab)
library(dplyr)


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
# Data extraction =================================
#**************************************************

prep_graphvar<-function(paths_=paths,
                        list_atlas_=list_atlas
                        ){
  
  for (atlas in list_atlas_){
    df_fc<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    list_node<-sort(unique(c(unique(df_fc$from),unique(df_fc$to))))
    
  }
}
