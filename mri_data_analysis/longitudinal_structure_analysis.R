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
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"))
#list_covar<-list("tanner"=list("1"="W1_Tanner_Full",
#                               "2"="W2_Tanner_Full",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild","value"=1)))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="Sex","value"=2)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="Sex","value"=2)))

subset_subj <- list("1"=list(list("key"="W1_T1QC_new_mild","value"=1),
                             list("key"="Sex","value"=2)),
                    "2"=list(list("key"="W2_T1QC_new_mild","value"=1),
                             list("key"="Sex","value"=2)))

#str_mod <- "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ s(age) + s(tanner,k=3) + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ age + tanner + s(ID_pnTTC,bs='re')"
str_mod <- "value ~ age*tanner + s(ID_pnTTC,bs='re')"

#list_str_group<-c("cortex","subcortex","white matter","global","misc")
list_str_group<-"subcortex"

#color<-"black"
#color<-"steelblue2"
color<-"lightcoral"

#key_global_covar<-"BrainSegVolNotVent"
#key_global_covar<-"eTIV"


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
                   list_wave_=list_wave,list_measure_=list_measure,list_str_group_=list_str_group,
                   str_mod_=str_mod,
                   color_=color
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
  gam.control(maxit=1000)
  df_out_term<-data.frame(matrix(nrow=0,ncol=9))
  colnames(df_out_term)<-c("measure","roi","label_roi","term","F","t","p")
  for (measure in list_measure_){
    print(paste('Calculating measurements of ',measure,sep=''))
    df_join_measure<-df_join[df_join$measure==measure,]
    list_roi<-as.character(unique(df_join_measure$roi))
    list_roi<-list_roi[order(list_roi)]
    for (roi in list_roi){
      label_roi<-as.character(dict_roi[dict_roi$id==roi,'label'])
      print(paste('Calculating ',roi,' = ',label_roi,sep=''))
      df_join_measure_roi<-df_join_measure[df_join_measure$roi==roi,]
      #mod_gamm<-gam(value ~ s(age) + s(tanner,k=3) + s(ID_pnTTC,bs='re'),data=df_join_measure_roi)
      formula_gamm<-as.formula(str_mod_)
      mod_gamm<-gam(formula_gamm,data=df_join_measure_roi)
      s_table<-summary.gam(mod_gamm)$s.table
      p_table<-summary.gam(mod_gamm)$p.table
      df_out_term_add<-rbind(data.frame(measure=measure,roi=roi,label_roi=label_roi,
                                        term=rownames(s_table),F=s_table[,'F'],t=NA,p=s_table[,'p-value']),
                             data.frame(measure=measure,roi=roi,label_roi=label_roi,
                                        term=rownames(p_table),F=NA,t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']))
      df_out_term<-rbind(df_out_term,df_out_term_add)
      for (covar in names(list_covar_)){
        plot<-plot_gamm(mod_gamm,covar,color_)
        label_covar<-list_covar_[[covar]][["label"]]
        plot<-(plot
               + ggtitle(paste('GAMM ',label_roi,sep=''))
               + xlab(label_covar)
               + ylab(capitalize(measure))
               + theme(legend.position = "none"))
        ggsave(paste("gamm_",measure,"_",roi,"_",covar,".eps",sep=""),plot=plot,device=cairo_ps,
               path=file.path(paths$output,"output"),dpi=300,height=5,width=5,limitsize=F)
      }
    }
  }
  rownames(df_out_term)<-NULL
  print('Saving results.')
  write.csv(df_out_term, file.path(paths_$output,"output","gamm.csv"),row.names = F)
  print('Finished gamm_str().')
}
