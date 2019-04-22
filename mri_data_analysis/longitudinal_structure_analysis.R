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
#list_measure <-c("volume","thickness","area")
list_measure <-"volume"
list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI"))
#list_covar<-list("tanner"=list("1"="W1_Tanner_Full",
#                               "2"="W2_Tanner_Full"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI"))

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild","value"=1),
                             list("key"="Sex","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild","value"=1),
                             list("key"="Sex","value"=1)))

#list_str_group<-c("cortex","subcortex","white matter","global","misc")
list_str_group<-"subcortex"

#key_global_covar<-"BrainSegVolNotVent"
key_global_covar<-"eTIV"



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

paths_=paths
subset_subj_=subset_subj
list_covar_=list_covar
file_input_=file_input
list_wave_=list_wave
list_measure_=list_measure
list_str_group_=list_str_group
key_global_covar_=key_global_covar


gamm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                   list_wave_=list_wave,list_measure_=list_measure,list_str_group_=list_str_group,
                   key_global_covar_=key_global_covar
                   ){
  print("Starting gamm_str().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  df_clin<-func_clinical_data_long(paths_,list_wave_)
  data_subset_clin<-func_subset_clin(df_clin,
                                     list_wave_,list_measure_,subset_subj_,
                                     list_covar,
                                     rem_na_clin=T)
  df_clin_subset<-data_subset_clin$df_clin
  
  # Load structural data
  print('Loading structural data.')
  df_str<-read.csv(file.path(paths_$input,"output",file_input_))
  df_str$value[which(is.nan(df_str$value))]<-0
  dict_roi<-func_dict_roi(paths_)
  df_str_subset<-func_subset_str(df_str,
                                 list_measure_,
                                 list_str_group_,
                                 dict_roi)
  
  # Join clinical and structural data frames
  print('Joining clinical and structural data.')
  df_join<-inner_join(df_str_subset,df_clin_subset,by=c('ID_pnTTC','wave'))
  df_join$ID_pnTTC<-as.factor(df_join$ID_pnTTC)
  df_join$wave<-as.factor(df_join$wave)
  df_join$measure<-as.factor(df_join$measure)
  
  # Calculate GAMM
  print('Calculating GAMM.')
  df_out<-data.frame(matrix(nrow=0,ncol=8))
  colnames(df_out)<-c("measure","roi","label_roi","term_smooth","F","p")
  for (measure in list_measure){
    df_join_measure<-df_join[df_join$measure==measure,]
    list_roi<-unique(df_join_measure$roi)
    list_roi<-list_roi[order(list_roi)]
    for (roi in list_roi){
      label_roi<-dict_roi[dict_roi$id==id,'label']
      df_join_measure_roi<-df_join_measure[df_join_measure$roi==roi,]
      mod_gamm<-gam(value ~ s(age) + s(tanner,k=3) + s(ID_pnTTC,bs='re'),data=df_join_measure_roi)
      s_table<-summary.gam(mod_gamm)$s.table
      df_out_add<-data.frame(measure=measure,roi=roi,label_roi=label_roi,
                             term_smooth=rownames(s_table),F=s_table[,'F'],p=s_table[,'p-value'])
      df_out<-rbind(df_out,df_out_add)
      plot<-plot_gamm(mod_gamm,'tanner')
      plot<-(plot
             + ggtitle('GAMM model')
             + xlab('Tanner stage')
             + ylab(measure))
      ggsave(paste("fc_heatmap_",atlas,"_",sprintf("%05d", id_subj),".eps",sep=""),plot=fig_fc_heatmap,device=cairo_ps,
             path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    }
  }
  
  df_join$ID_pnTTC<-as.factor(df_join$ID_pnTTC)
  df_join$wave<-as.factor(df_join$wave)
  df_join$measure<-as.factor(df_join$measure)
  
}
