#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal graph measures.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/graph_GV"
dir_in <-"203_graph_acompcor"
dir_out <-"204_gamm_graph_acompcor"

list_wave<-c(1,2)

#list_metric_local=c("degrees_und","efficiency_local_bin")
#list_metric_global=c("efficiency_bin","charpath_B_radius","charpath_B_diameter")
list_metric_local=c("efficiency_local_bin")
list_metric_global=c("efficiency_bin","efficiency_local_bin")

#list_covar<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
#                 "age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                 "sex"=list("1"="Sex","2"="Sex","label"="Sex"))
list_covar<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                 "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                 "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                 "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"),
                 "age"  =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                 "sex"  =list("1"="Sex",            "2"="Sex",            "label"="Sex"))
subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                             list("key"="W1_rsfMRIexist","condition"="==1"),
                             list("key"="W1_Censor","condition"="<126")),
                    "2"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))
list_mod <- list("lin"  ="value ~ age + testo + s(ID_pnTTC,bs='re')",
                 "add"  ="value ~ s(age,k=3) + s(testo,k=3) + s(ID_pnTTC,bs='re')",
                 "quad" ="value ~ poly(age,2) + poly(testo,2) + s(ID_pnTTC,bs='re')")
list_graph <-list("t"=list("title"="Testosterone effect","x_axis"="testo",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
list_mod_diff <- list("lin_diff"="value ~ diff_age + diff_testo",
                      "lin_diff_mean"="value ~ diff_age + diff_testo + mean_testo",
                      "add_diff"="value ~ s(diff_age,k=3) + s(diff_testo,k=3)",
                      "add_diff_mean"="value ~ s(diff_age,k=3) + s(mean_testo,k=3) + s(diff_testo,k=3)")
list_graph_diff <-list("diff"=list("title"="Testosterone diff effect",
                                   "x_axis"="diff_testo",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                       "mean"=list("title"="Testosterone mean effect",
                                   "x_axis"="mean_testo",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))

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

gamm_core<-function(df_src,roi,label_roi,group,measure,list_mod_,list_graph_,list_covar_,paths_){
  print(paste("GLM/GAM ",measure,' of ',roi," (",label_roi,")",sep=""))
  df_out_aic_add<-df_out_lm_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    list_plot<-list()
    list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
    for (idx_sex in list_sex){
      df_src_sex<-df_src[df_src$sex==idx_sex,]
      #mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex)
      mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
      if (class(mod)[1]=="try-error"){
        print(paste("Error fiting ",idx_mod, ", sex= ",idx_sex,".",sep=''))
      }else{
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
                            df_in=df_src_sex,
                            spec_graph=list_graph_[[idx_graph]])
            list_plot[[idx_graph]]<-plot
            
            # Output
            if (idx_sex==list_sex[length(list_sex)]){
              idx_axis_x<-list_graph_[[idx_graph]][["x_axis"]]
              if (substring(idx_axis_x,1,5)=="diff_"){
                label_axis_x<-list_covar_[[substring(idx_axis_x,6)]][["label"]]
                label_axis_x<-paste(label_axis_x,"Difference",sep=" ")
              }else if(substring(idx_axis_x,1,5)=="mean_"){
                label_axis_x<-list_covar_[[substring(idx_axis_x,6)]][["label"]]
                label_axis_x<-paste(label_axis_x,"Mean",sep=" ")
              }else{
                label_axis_x<-list_covar_[[idx_axis_x]][["label"]]
              }
              plot<-(plot
                     + ggtitle(paste(list_graph_[[idx_graph]][["title"]],' on ',label_roi,' ',measure,'\n',idx_mod,sep=''))
                     + xlab(label_axis_x)
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
  }
  
  # Compare AICs of GAMM models
  df_out_aic_add_sex_rbind<-data.frame()
  for (idx_sex in list_sex){
    df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
    df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                       'aic_best_among_models']<-1
    df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
  }
  
  return(list("df_out_lm_add"=df_out_lm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind))
}

paths_=paths
subset_subj_=subset_subj
list_covar_=list_covar
list_wave_=list_wave
list_metric_global_=list_metric_global
list_metric_local_=list_metric_local
list_mod_=list_mod
list_graph_=list_graph
list_mod_diff_=list_mod_diff
list_graph_diff_=list_graph_diff
gamm_gta<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,list_wave_=list_wave,
                   list_metric_global_=list_metric_global,list_metric_local_=list_metric_local,
                   list_mod_=list_mod,list_graph_=list_graph,list_mod_diff_=list_mod_diff,
                   list_graph_diff_=list_graph_diff
                   ){
  print("Starting gamm_gta().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  dict_roi<-func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
  # Load global graph data
  print('Loading global graph data.')
  df_gta_global<-read.csv(file.path(paths_$input,"output",paste("atl-power264_graph_global.csv")))
  colnames(df_gta_global)[colnames(df_gta_global)=="ses"]<-"wave"
  
  # Join clinical and global graph data frames
  print('Joining clinical and global graph data.')
  df_join<-inner_join(df_gta_global,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  write.csv(df_join,file.path(paths_$output,"output","src_global.csv"),row.names=F)
  
  # Calculate GAMM
  print('Calculating GAMM of global data.')
  df_out_global<-df_out_global_aic<-NULL
  for(metric in list_metric_global_){
    df_join_metric<-df_join[df_join$metric==metric,]

    out_glm<-gamm_core(df_src=df_join_metric,roi="global",label_roi="Global Metric",
                       group="whole",measure=metric,
                       list_mod_,list_graph_,list_covar_,paths_)
    df_out_global<-rbind(df_out_global,out_glm$df_out_lm_add)
    df_out_global_aic<-rbind(df_out_global_aic,out_glm$df_out_aic_add)
  }
  
  # Data saving
  rownames(df_out_global)<-rownames(df_out_global_aic)<-NULL
  write.csv(df_out_global, file.path(paths_$output,"output","gamm_global.csv"),row.names = F)
  write.csv(df_out_global_aic,file.path(paths_$output,"output","gamm_global_aic.csv"),row.names = F)
  
  # Load local graph data
  print('Loading local graph data.')
  df_gta_local<-read.csv(file.path(paths_$input,"output",paste("atl-power264_graph_local.csv")))
  colnames(df_gta_local)[colnames(df_gta_local)=="ses"]<-"wave"
  
  # Join clinical and local graph data frames
  print('Joining clinical and local graph data.')
  df_join<-inner_join(df_gta_local,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  write.csv(df_join,file.path(paths_$output,"output","src_local.csv"),row.names=F)
  
  # Calculate GAMM
  print('Calculating GAMM of local data.')
  df_out_local<-df_out_local_aic<-NULL
  for(metric in list_metric_local_){
    df_join_metric<-df_join[df_join$metric==metric,]
    list_roi<-sort(unique(as.character(df_join_metric$roi)))
    for (roi in list_roi){
      label_roi=as.character(dict_roi[dict_roi$id==roi,"label"])
      df_src=df_join_metric[df_join_metric$roi==roi,]
      if (dim(df_src)[2]>0){
        out_glm<-gamm_core(df_src,roi,label_roi,group="whole",measure=metric,
                           list_mod_,list_graph_,list_covar_,paths_)
        df_out_local<-rbind(df_out_local,out_glm$df_out_lm_add)
        df_out_local_aic<-rbind(df_out_local_aic,out_glm$df_out_aic_add)
      }
    }
  }
  
  # Data saving
  rownames(df_out_local)<-rownames(df_out_local_aic)<-NULL
  write.csv(df_out_local, file.path(paths_$output,"output","gamm_local.csv"),row.names = F)
  write.csv(df_out_local_aic,file.path(paths_$output,"output","gamm_local_aic.csv"),row.names = F)
  
  #************************************************
  # Analyses using parameter difference and average 
  #************************************************
  # Prepare clinical data difference and average dataframe
  list_subj<-list()
  for (wave in list_wave_){
    list_subj<-c(list_subj,
                 list(sort(unique(as.numeric.factor(df_join[df_join$wave==wave,"ID_pnTTC"])))))
  }
  list_id_subj_twice<-intersect(list_subj[[1]],list_subj[[2]])
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  df_clin_diff<-func_clinical_data_join(df_src=df_clin,
                                        list_id_subj=list_id_subj_twice,
                                        list_covar=list_covar_)
  colnames(df_clin_diff)[colnames(df_clin_diff)=="ses"]<-"wave"
  for (key in c('ID_pnTTC','sex')){
    if (key %in% colnames(df_clin_diff)){
      df_clin_diff[,key]<-as.factor(df_clin_diff[,key])
    }
  }
  
  # Calculate GAMM
  print('Calculating GAMM of global data change.')
  df_out_global<-df_out_global_aic<-NULL
  for(metric in list_metric_global_){
    df_gta_global_metric<-df_gta_global[df_gta_global$metric==metric,]
    df_join_metric<-df_clin_diff
    for (idx_row in seq(nrow(df_clin_diff))){
      id_subj<-as.numeric.factor(df_clin_diff[idx_row,"ID_pnTTC"])
      value_ses1<-df_gta_global_metric[df_gta_global_metric$ID_pnTTC==id_subj
                                       & df_gta_global_metric$wave==1,"value"]
      value_ses2<-df_gta_global_metric[df_gta_global_metric$ID_pnTTC==id_subj
                                       & df_gta_global_metric$wave==2,"value"]
      df_join_metric[idx_row,"value"]<-value_ses2-value_ses1
    }
    
    out_glm<-gamm_core(df_src=df_join_metric,roi="global",label_roi="Global Metric",
                       group="whole",measure=metric,
                       list_mod_diff_,list_graph_diff_,list_covar_,paths_)
    df_out_global<-rbind(df_out_global,out_glm$df_out_lm_add)
    df_out_global_aic<-rbind(df_out_global_aic,out_glm$df_out_aic_add)
  }
  
  # Data saving
  rownames(df_out_global)<-rownames(df_out_global_aic)<-NULL
  write.csv(df_out_global, file.path(paths_$output,"output","gamm_global_diff.csv"),row.names = F)
  write.csv(df_out_global_aic,file.path(paths_$output,"output","gamm_global_diff_aic.csv"),row.names = F)
  
  print('Finished gamm_gta().')
}
