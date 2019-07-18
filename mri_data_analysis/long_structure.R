#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal structural properties from FreeSurfer data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS"
dir_in <-"01_extract"
#dir_out <-"03_gamm"
dir_out <-"06_gamm_smooth"
#dir_out <-"08_gamm_subcortex_male"
#dir_out <-"10_gamm_subcortex_female_misc"
file_input<-"fs_measure.csv"

list_wave <- c(1,2)

#list_measure <-c("volume","thickness","area")
list_measure <-"volume"

#list_str_group<-c("cortex","subcortex","white matter","global","misc")
#list_str_group<-"subcortex"
#list_str_group<-c("global","misc")
list_str_group<-c("cortex","subcortex","global")

#list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
#                               "2"="W2_Tanner_Max",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"))
#list_covar<-list("tanner"=list("1"="W1_Tanner_Full",
#                               "2"="W2_Tanner_Full",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"))
#list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
#                               "2"="W2_Tanner_Max",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"),
#                 "sex"=list("1"="Sex",
#                            "2"="Sex",
#                            "label"="Sex"))

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=2)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=2)))

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

#list_mod <- list("a+s+st"=
#                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')",
#                 "a+s+st+sat"=
#                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + ti(age,tanner,k=2,by=sex) + s(ID_pnTTC,bs='re')")

list_mod <- list("a+s+st"=
                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')")

#list_mod <- list("a+t"=
#                   "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')")

#list_mod <- list("a+t"=
#                   "value ~ s(age,k=3) + s(tanner,k=3)")

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

#list_graph <-list("a"=list("title"="Age effect",
#                           "x_axis"="age",
#                           "smooth"=list("Male"=list("fix"=list("sex"=1),
#                                                     "color"="steelblue2","alpha"=1,"ribbon"=T)),
#                           "point"=list("Male"=list("subset"=NULL,
#                                                    "color"="steelblue2","alpha"=1))),
#                  "t"=list("title"="Tanner stage effect",
#                           "x_axis"="tanner",
#                           "smooth"=list("Male"=list("fix"=list("sex"=1),
#                                                     "color"="steelblue2","alpha"=1,"ribbon"=T)),
#                           "point"=list("Male"=list("subset"=NULL,
#                                                    "color"="steelblue2","alpha"=1))))

#list_graph <-list("a"=list("title"="Age effect",
#                           "x_axis"="age",
#                           "smooth"=list("Female"=list("fix"=list("sex"=2),
#                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
#                           "point"=list("Female"=list("subset"=NULL,
#                                                      "color"="lightcoral","alpha"=1))),
#                  "t"=list("title"="Tanner stage effect",
#                           "x_axis"="tanner",
#                           "smooth"=list("Female"=list("fix"=list("sex"=2),
#                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
#                           "point"=list("Female"=list("subset"=NULL,
#                                                       "color"="lightcoral","alpha"=1))))



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
source(file.path(paths$script,"util/function.R"))
#source(file.path(paths$script,"util/glm_function.R"))
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# GAMM of structural measures =====================
#**************************************************

#paths_=paths
#subset_subj_=subset_subj
#list_covar_=list_covar
#file_input_=file_input
#list_wave_=list_wave
#list_measure_=list_measure
#list_str_group_=list_str_group
#list_mod_=list_mod
#list_graph_=list_graph


gamm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                   list_wave_=list_wave,list_measure_=list_measure,list_str_group_=list_str_group,
                   list_mod_=list_mod,list_graph_=list_graph
                   ){
  print("Starting gamm_str().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  df_clin<-func_clinical_data_long(paths_,list_wave_)
  data_subset_clin<-func_subset_clin(df_clin,
                                     list_wave_,list_measure_,subset_subj_,
                                     list_covar_,
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
  for (key in c('ID_pnTTC','wave','sex','measure')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  
  # Calculate GAMM
  print('Calculating GAMM.')
  df_out_term<-data.frame(matrix(nrow=0,ncol=9))
  colnames(df_out_term)<-c("measure","roi","label_roi","group_roi","model","term","F","t","p")
  df_out_model<-data.frame(matrix(nrow=0,ncol=7))
  colnames(df_out_model)<-c("measure","roi","label_roi","group_roi","model","aic","aic_best_among_models")
  for (measure in list_measure_){
    print(paste('Calculating measurements of ',measure,sep=''))
    df_join_measure<-df_join[df_join$measure==measure,]
    list_roi<-as.character(unique(df_join_measure$roi))
    list_roi<-list_roi[order(list_roi)]
    for (roi in list_roi){
      label_roi<-as.character(dict_roi[dict_roi$id==roi,'label'])
      group_roi<-as.character(dict_roi[dict_roi$id==roi,'group'])
      print(paste('Calculating ',roi,' = ',label_roi,sep=''))
      df_join_measure_roi<-df_join_measure[df_join_measure$roi==roi,]
      list_mod_gamm<-list()
      df_out_model_add<-data.frame()
      for (mod in names(list_mod_)){
        list_mod_gamm[[mod]]<-gam(as.formula(list_mod_[[mod]]),data=df_join_measure_roi)
        p_table<-summary.gam(list_mod_gamm[[mod]])$p.table
        s_table<-summary.gam(list_mod_gamm[[mod]])$s.table
        df_out_term_add<-rbind(data.frame(measure=measure,roi=roi,label_roi=label_roi,group_roi=group_roi,model=mod,
                                          term=rownames(p_table),F=NA,t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                               data.frame(measure=measure,roi=roi,label_roi=label_roi,group_roi=group_roi,model=mod,
                                          term=rownames(s_table),F=s_table[,'F'],t=NA,p=s_table[,'p-value']))
        df_out_term<-rbind(df_out_term,df_out_term_add)
        df_out_model_add<-rbind(df_out_model_add,
                                data.frame(measure=measure,roi=roi,label_roi=label_roi,group_roi=group_roi,model=mod,
                                           aic=list_mod_gamm[[mod]]$aic,aic_best_among_models=0))
        
        for (idx_graph in names(list_graph_)){
          plot<-plot_gamm(mod_gamm=list_mod_gamm[[mod]],
                          df_join_measure_roi,
                          spec_graph=list_graph_[[idx_graph]])
          axis_x<-list_graph_[[idx_graph]][["x_axis"]]
          label_x<-list_covar_[[axis_x]][["label"]]
          plot<-(plot
                 + ggtitle(paste(list_graph_[[idx_graph]][["title"]],label_roi,sep=' '))
                 + xlab(label_x)
                 + ylab(capitalize(measure))
                 + theme(legend.position = "none"))
          filename_plot<-paste("gamm_",measure,"_",roi,"_",mod,"_",idx_graph,".eps",sep="")
          ggsave(filename_plot,plot=plot,device=cairo_ps,
                 path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
          
        }
      }
      # compare AICs of models
      df_out_model_add[which(df_out_model_add$aic==min(df_out_model_add$aic)),'aic_best_among_models']<-1
      df_out_model<-rbind(df_out_model,df_out_model_add)
      rownames(df_out_term)<-rownames(df_out_model)<-NULL
      write.csv(df_out_term, file.path(paths_$output,"output","gamm.csv"),row.names = F)
      write.csv(df_out_model,file.path(paths_$output,"output","aic.csv"),row.names = F)
      # compare models
      #if (length(list_mod_)==2){
      #  anova_mod<-anova.gam(list_mod_gamm[[1]],list_mod_gamm[[2]],test="F")
      #  df_out_model[dim(df_out_model)[1]+1,]<-c(measure,roi,label_roi,group_roi,anova_mod[2,"F"],anova_mod[2,"Pr(>F)"])
      #}
    }
  }
  #print('Saving results.')
  #rownames(df_out_term)<-rownames(df_out_model)<-NULL
  #write.csv(df_out_term, file.path(paths_$output,"output","gamm.csv"),row.names = F)
  #write.csv(df_out_model,file.path(paths_$output,"output","aic.csv"),row.names = F)
  print('Finished gamm_str().')
}
