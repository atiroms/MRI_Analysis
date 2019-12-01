#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal structural properties from FreeSurfer data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/Puberty/Stats/str_FS"
dir_in <-"01_extract"
dir_out <-"02_glm_test"

file_input<-"fs_measure.csv"

list_wave <- c(1,2)

#list_measure <-c("volume","thickness","area")
list_measure <-"volume"

#list_str_group<-c("cortex","subcortex","white matter","global","misc")
list_str_group<-"subcortex"
#list_str_group<-c("global","misc")
#list_str_group<-c("cortex","subcortex","global")

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage (max)"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))


subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild","value"=1)))

#list_mod <- list("a+s+st"=
#                   "value ~ age + sex + sex:tanner + s(ID_pnTTC,bs='re')",
#                 "a+s+st+sat"=
#                   "value ~ age + sex + sex:tanner + sex:age:tanner + s(ID_pnTTC,bs='re')",
#                 "a+s+st2"=
#                   "value ~ age + sex + sex:poly(tanner,2) + s(ID_pnTTC,bs='re')",
#                 "a+s+st2+sat2"=
#                   "value ~ age + sex + sex:poly(tanner,2) + sex:age:poly(tanner,2) + s(ID_pnTTC,bs='re')",
#                 "a2+s+st"=
#                   "value ~ poly(age,2) + sex + sex:tanner + s(ID_pnTTC,bs='re')",
#                 "a2+s+st+sa2t"=
#                   "value ~ poly(age,2) + sex + sex:tanner + sex:poly(age,2):tanner + s(ID_pnTTC,bs='re')",
#                 "a2+s+st2"=
#                   "value ~ poly(age,2) + sex + sex:poly(tanner,2) + s(ID_pnTTC,bs='re')",
#                 "a2+s+st2+sa2t2"=
#                   "value ~ poly(age,2) + sex + sex:poly(tanner,2) + sex:poly(age,2):poly(tanner,2) + s(ID_pnTTC,bs='re')")

#list_mod <- list("gam_ast"=
#                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')",
#                 "glm_ast"=
#                   "value ~ age + sex + sex:tanner + s(ID_pnTTC,bs='re')",
#                 "gam_at"=
#                   "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')",
#                 "glm_at"=
#                   "value ~ age + tanner + s(ID_pnTTC,bs='re')")
list_mod <- list("gam_at"=
                   "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')",
                 "glm_at"=
                   "value ~ age + tanner + s(ID_pnTTC,bs='re')")

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
                  "st"=list("title"="Tanner stage effect",
                            "x_axis"="tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "sat"=list("title"="Age-Tanner stage interaction",
                             "x_axis"="age",
                             "smooth"=list("Male TS = 1"=list("fix"=list("sex"=1,"tanner"=1),
                                                              "color"="Steelblue2","alpha"=0.4,"ribbon"=F),
                                           "Male TS = 3"=list("fix"=list("sex"=1,"tanner"=3),
                                                              "color"="steelblue2","alpha"=0.7,"ribbon"=F),
                                           "Male TS = 5"=list("fix"=list("sex"=1,"tanner"=5),
                                                              "color"="steelblue2","alpha"=1,"ribbon"=F),
                                           "Female TS = 1"=list("fix"=list("sex"=2,"tanner"=1),
                                                                "color"="lightcoral","alpha"=0.4,"ribbon"=F),
                                           "Female TS = 3"=list("fix"=list("sex"=2,"tanner"=3),
                                                                "color"="lightcoral","alpha"=0.7,"ribbon"=F),
                                           "Female TS = 5"=list("fix"=list("sex"=2,"tanner"=5),
                                                                "color"="lightcoral","alpha"=1,"ribbon"=F)),
                             "point"=list("Male"=list("subset"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1),
                                          "Female"=list("subset"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1))))


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
#source(file.path(paths$script,"util/glm_function.R"))
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# GAMM of structural measures =====================
#**************************************************

glm_core<-function(df_src,roi,label_roi,group,measure,list_mod_,list_graph_,list_covar_,paths_){
  print(paste("GLM/GAM Group: ",group,", ROI: ",roi,"=",label_roi,", Measure: ",measure,".",  sep=""))
  df_out_aic_add<-df_out_lm_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    list_plot<-list()
    list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
    for (idx_sex in list_sex){
      df_src_sex<-df_src[df_src$sex==idx_sex,]
      mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex)
      p_table<-summary.gam(mod)$p.table
      if (is.null(summary.gam(mod)$s.table)){
        df_out_lm_add_add<-data.frame(roi=roi,label_roi=label_roi,group=group,measure=measure,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                      estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                      t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
        
      }else{
        s_table<-summary.gam(mod)$s.table
        df_out_lm_add_add<-rbind(data.frame(roi=roi,label_roi=label_roi,group=group,measure=measure,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                            estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                            t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                                 data.frame(roi=roi,label_roi=label_roi,group=group,measure=measure,sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                            estimate=NA,se=NA,F=s_table[,'F'],
                                            t=NA,p=s_table[,'p-value']))
      }
      df_out_lm_add<-rbind(df_out_lm_add,df_out_lm_add_add)
      df_out_aic_add<-rbind(df_out_aic_add,
                            data.frame(roi=roi,label_roi=label_roi,group=group,measure=measure,sex=idx_sex,
                                       model=idx_mod,aic=mod$aic,aic_best_among_models=0))
      
      # Graphical output of GLM/GAM results
      for (idx_graph in names(list_graph_)){
        if (list_graph_[[idx_graph]][["x_axis"]] %in% colnames(mod$model)){
          # Add sex-wise lines/plots to existent plot, initialize if absent
          plot<-plot_gamm(plot_in=list_plot[[idx_graph]],mod_gamm=mod,
                          df_join_measure_roi=df_src_sex,
                          spec_graph=list_graph_[[idx_graph]])
          list_plot[[idx_graph]]<-plot
          
          # Output
          if (idx_sex==list_sex[length(list_sex)]){
            plot<-(plot
                   + ggtitle(paste(list_graph_[[idx_graph]][["title"]],label_roi,measure,idx_mod,sep=" "))
                   + xlab(capitalize(list_graph_[[idx_graph]][["x_axis"]]))
                   + ylab(capitalize(measure))
                   + theme(legend.position = "none"))
            filename_plot<-paste("roi-",roi,"_msr-",measure,"_mod-",idx_mod,
                                 "_plt-",idx_graph,"_glm.eps",sep="")
            ggsave(filename_plot,plot=plot,device=cairo_ps,
                   path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
          }
        }
      }
    }
  }
  
  # Compare AICs of GLM models
  df_out_aic_add_sex_rbind<-data.frame()
  for (idx_sex in list_sex){
    df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
    df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                       'aic_best_among_models']<-1
    df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
  }
  
  return(list("df_out_lm_add"=df_out_lm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind))
}


glm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                  list_wave_=list_wave,list_measure_=list_measure,list_str_group_=list_str_group,
                  list_mod_=list_mod,list_graph_=list_graph,key_group_='group_3'
                  ){
  print("Starting gamm_str().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
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
  df_join<-inner_join(df_str_subset,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex','measure')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  
  write.csv(df_join,file.path(paths_$output,"output","src.csv"),row.names=F)
  
  # Calculate GAMM
  print('Calculating GAMM.')
  df_out_lm<-df_out_aic<-NULL
  
  for (measure in list_measure_){
    for (str_group in list_str_group_){
      df_join_measure_group<-df_join[df_join$measure==measure & df_join$group==str_group,]
      list_roi<-sort(unique(as.character(df_join_measure_group$roi)))
      for (roi in list_roi){
        label_roi=as.character(dict_roi[dict_roi$id==roi,"label"])
        df_src=df_join_measure_group[df_join_measure_group$roi==roi,]
        if (dim(df_src)[2]>0){
          out_glm<-glm_core(df_src,roi,label_roi,group=str_group,measure,
                            list_mod_,list_graph_,list_covar_,paths_)
          df_out_lm<-rbind(df_out_lm,out_glm$df_out_lm_add)
          df_out_aic<-rbind(df_out_aic,out_glm$df_out_aic_add)
        }
      }
    }
  }
  
  # Data saving
  rownames(df_out_lm)<-rownames(df_out_aic)<-NULL
  write.csv(df_out_lm, file.path(paths_$output,"output","glm.csv"),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output","glm_aic.csv"),row.names = F)
  
  print('Finished gamm_str().')
}
