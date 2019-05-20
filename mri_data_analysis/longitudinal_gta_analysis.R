#**************************************************
# Description =====================================
#**************************************************
# R script for longitudinal analysis of graph theoretical measures.
# Inputs are gta measures computed with gta_bin()


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for gamm_gta()
path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_src<-"57_gta_weight"
dir_dst<-"58_gamm_gta_weight"

list_wave <- c(1,2)

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"
list_atlas<-"glasser360"

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

list_mod <- list("additive"=
                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex)",
                 "additive_mixed"=
                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')",
                 "linear"=
                   "value ~ age + sex + sex:tanner",
                 "linear_mixed"=
                   "value ~ age + sex + sex:tanner + s(ID_pnTTC,bs='re')")

list_plot<-list("a"=list("title"="Age effect",
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
                                                     "color"="lightcoral","alpha"=1))))
                #"sat"=list("title"="Age-Tanner stage interaction",
                #           "x_axis"="age",
                #           "smooth"=list("Male TS = 1"=list("fix"=list("sex"=1,"tanner"=1),
                #                                            "color"="Steelblue2","alpha"=0.4,"ribbon"=F),
                #                         "Male TS = 3"=list("fix"=list("sex"=1,"tanner"=3),
                #                                            "color"="steelblue2","alpha"=0.7,"ribbon"=F),
                #                         "Male TS = 5"=list("fix"=list("sex"=1,"tanner"=5),
                #                                            "color"="steelblue2","alpha"=1,"ribbon"=F),
                #                         "Female TS = 1"=list("fix"=list("sex"=2,"tanner"=1),
                #                                              "color"="lightcoral","alpha"=0.4,"ribbon"=F),
                #                         "Female TS = 3"=list("fix"=list("sex"=2,"tanner"=3),
                #                                              "color"="lightcoral","alpha"=0.7,"ribbon"=F),
                #                         "Female TS = 5"=list("fix"=list("sex"=2,"tanner"=5),
                #                                              "color"="lightcoral","alpha"=1,"ribbon"=F)),
                #           "point"=list("Male"=list("subset"=list("sex"=1),
                #                                    "color"="steelblue2","alpha"=1),
                #                        "Female"=list("subset"=list("sex"=2),
                #                                      "color"="lightcoral","alpha"=1)))
                #)


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
                    dir_src_=dir_src,
                    dir_dst_=dir_dst){
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
  path_src     <- file.path(path_root,path_exp_,dir_src_)
  path_dst    <- file.path(path_root,path_exp_,dir_dst_)
  output <- list("script"=path_script,"input"=path_src,"output"=path_dst,
                 "common"=path_common,"dir_in"=dir_src_,"dir_out"=dir_dst_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))
#source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
#source(file.path(paths$script,"functionality/gta_function.R"))


#**************************************************
# GAMM function ===================================
#**************************************************
paths_=paths
subset_subj_=subset_subj
list_covar_=list_covar
list_wave_=list_wave
list_mod_=list_mod
list_plot_=list_plot
list_atlas_=list_atlas
atlas=list_atlas_[1]

gamm_gta<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                   list_wave_=list_wave,
                   list_mod_=list_mod,list_plot_=list_plot,list_atlas_=list_atlas
                   ){
  print("Starting gamm_gta().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  df_clin<-func_clinical_data_long(paths_,list_wave_)
  data_subset_clin<-func_subset_clin(df_clin,
                                     list_wave_,subset_subj_,
                                     list_covar_,
                                     rem_na_clin=T)
  df_clin_subset<-data_subset_clin$df_clin
  
  dict_roi<-func_dict_roi(paths_)
  
  for (atlas in list_atlas_){
    print(paste('Calculating atlas: ',atlas,sep=''))
    # Load GTA data
    print('Loading GTA data.')
    file_src<-paste("atl-",atlas,"_gta_weight.csv",sep="")
    df_gta<-read.csv(file.path(paths_$input,"output",file_src))
    colnames(df_gta)[colnames(df_gta)=="ses"]<-"wave"
    df_gta$value[which(is.nan(df_gta$value))]<-0

    # Join clinical and structural data frames
    print('Joining clinical and GTA data.')
    df_join<-inner_join(df_gta,df_clin_subset,by=c('ID_pnTTC','wave'))
    for (key in c('ID_pnTTC','wave','sex','measure')){
      if (key %in% colnames(df_join)){
        df_join[,key]<-as.factor(df_join[,key])
      }
    }
    
    # Calculate GAMM
    print('Calculating GAMM.')
    df_out_term<-data.frame(matrix(nrow=0,ncol=10))
    colnames(df_out_term)<-c("cost","node","label_node","group_node","metric","model","term","F","t","p")
    df_out_model<-data.frame(matrix(nrow=0,ncol=8))
    colnames(df_out_model)<-c("cost","node","label_node","group_node","metric","model","aic","aic_best_among_models")
    
    list_node<-as.character(unique(df_join$node))
    list_node<-list_node[order(list_node)]
    list_node<-c("graph",list_node[list_node!="graph"])
    
    list_metric_node<-as.character(unique(df_join[df_join$node==list_node[2],"metric"]))
    list_metric_node<-list_metric_node[order(list_metric_node)]
    list_metric_graph<-as.character(unique(df_join[df_join$node=="graph","metric"]))
    list_metric_graph<-list_metric_graph[order(list_metric_graph)]
    
    #list_cost<-as.character(unique(df_join$cost))
    #list_cost<-list_cost[order(list_cost)]
    #list_cost<-c("average",list_cost[list_cost!="average"])
    list_cost<-"average"
    for (cost in list_cost){
      if ("cost" %in% colnames(df_join)){
        df_join_cost<-df_join[df_join$cost==cost,]
      }else{
        df_join_cost<-df_join
      }
      print(paste('Calculating cost ',cost,sep=''))
      for (node in list_node){
        if (node=="graph"){
          label_node<-"Graph"
          group_node<-"global"
          list_metric<-list_metric_graph
        }else{
          label_node<-as.character(dict_roi[dict_roi$id==node,'label'])
          group_node<-as.character(dict_roi[dict_roi$id==node,'group'])
          list_metric<-list_metric_node
        }
        print(paste('Calculating ',node,' (',label_node,')',sep=''))
        df_join_cost_node<-df_join_cost[df_join_cost$node==node,]
        for (metric in list_metric){
          df_join_cost_node_metric<-df_join_cost_node[df_join_cost_node$metric==metric,]
          list_mod_gamm<-list()
          df_out_model_add<-data.frame()
          for (mod in names(list_mod_)){
            list_mod_gamm[[mod]]<-try(gam(as.formula(list_mod_[[mod]]),data=df_join_cost_node_metric),silent=T)
            if(class(list_mod_gamm[[mod]])[[1]]!="try-error"){
              p_table<-summary.gam(list_mod_gamm[[mod]])$p.table
              s_table<-summary.gam(list_mod_gamm[[mod]])$s.table
              if (!is.null(s_table)){
                df_out_term_add<-rbind(data.frame(cost=cost,node=node,label_node=label_node,group_node=group_node,metric=metric,model=mod,
                                                  term=rownames(p_table),F=NA,t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                                       data.frame(cost=cost,node=node,label_node=label_node,group_node=group_node,metric=metric,model=mod,
                                                  term=rownames(s_table),F=s_table[,'F'],t=NA,p=s_table[,'p-value']))
              }
              df_out_term<-rbind(df_out_term,df_out_term_add)
              df_out_model_add<-rbind(df_out_model_add,
                                      data.frame(cost=cost,node=node,label_node=label_node,group_node=group_node,metric=metric,model=mod,
                                                 aic=list_mod_gamm[[mod]]$aic,aic_best_among_models=0))
              
              for (idx_plot in names(list_plot_)){
                plot<-plot_gamm(mod_gamm=list_mod_gamm[[mod]],
                                df_join_measure_roi=df_join_cost_node_metric,
                                spec_graph=list_plot_[[idx_plot]])
                axis_x<-list_plot_[[idx_plot]][["x_axis"]]
                label_x<-list_covar_[[axis_x]][["label"]]
                plot<-(plot
                       + ggtitle(paste(list_plot_[[idx_plot]][["title"]],":",label_node,sep=' '))
                       + xlab(label_x)
                       + ylab(capitalize(metric))
                       + theme(legend.position = "none"))
                filename_plot<-paste("atl-",atlas,"_cost-",cost,"_node-",node,"_metric-",metric,"_mod-",mod,"_plt-",idx_plot,"_gamm.eps",sep="")
                ggsave(filename_plot,plot=plot,device=cairo_ps,
                       path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
                
              }
            }
          }
          # compare AICs of models
          if (length(df_out_model_add)>0){
            df_out_model_add[which(df_out_model_add$aic==min(df_out_model_add$aic)),'aic_best_among_models']<-1
            df_out_model<-rbind(df_out_model,df_out_model_add)
          }
          rownames(df_out_term)<-rownames(df_out_model)<-NULL
          fileprefix_out<-paste("atl-",atlas,"_gamm.csv",sep="")
          write.csv(df_out_term, file.path(paths_$output,"output",paste("atl-",atlas,"_gamm.csv",sep="")),row.names = F)
          write.csv(df_out_model,file.path(paths_$output,"output",paste("atl-",atlas,"_aic.csv",sep="")),row.names = F)
          # compare models
          #if (length(list_mod_)==2){
          #  anova_mod<-anova.gam(list_mod_gamm[[1]],list_mod_gamm[[2]],test="F")
          #  df_out_model[dim(df_out_model)[1]+1,]<-c(measure,roi,label_node,group_node,anova_mod[2,"F"],anova_mod[2,"Pr(>F)"])
          #}
        }
      }
    }
    print(paste("Finished Calculating atlas: ",atlas, sep=""))
  }
  print("Finished gamm_gta()")
}
