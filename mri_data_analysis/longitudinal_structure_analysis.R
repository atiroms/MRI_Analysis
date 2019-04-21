#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal structural properties from FreeSurfer data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS"
dir_in <-"01_extract"
dir_out <-"03_gamm"
file_input<-"fs_measure.csv"

list_wave <- c(1,2)
list_measure <-c("volume","thickness","area")
list_covar<-list("tanner_max"=c("W1_Tanner_Max","W2_Tanner_Max"),
                 "age"=c("W1_Age_at_MRI","W2_Age_at_MRI"))

#key_global_covar<-"BrainSegVolNotVent"
key_global_covar<-"eTIV"

subset_subj <- list("1"=list(list("key"="W1_T1QC_new_mild","value"=1),list("key"="Sex","value"=1)),
                    "2"=list(list("key"="W2_T1QC_new_mild","value"=1),list("key"="Sex","value"=1)))
#subset_subj <- list(list("column"="W1_T1QC","value"=1),list("column"="Sex","value"=2))
#subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_5sub","value"=1),list("column"="Sex","value"=1))


#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(tidyverse)
library(tidyr)
library(mgcv)
library(dplyr)
library(ggplot2)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms"),
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
  path_common <- file.path(path_root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
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
source(file.path(paths$script,"functionality/function.R"))
#source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
#source(file.path(script_dir,"Functionalities/LI_Functions.R"))


#**************************************************
# GAMM of structural measures =====================
#**************************************************
gamm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                   list_wave_=list_wave,list_measure_=list_measure,
                   key_global_covar_=key_global_covar
                   ){
  print("Starting gamm_str().")
  data_clinical<-func_clinical_data(paths_,wave,subset_subj_)
  
  nullobj<-func_createdirs(paths_,copy_log=T)
  dict_roi<-func_dict_roi(paths_)
  df_str<-read.csv(file.path(paths_$input,"output",file_input_))
  df_str$value[which(is.nan(df_str$value))]<-0
  for (wave in list_wave){
    
  }
  
  
}


