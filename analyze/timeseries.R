#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for fc()
#path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in <-"300_ts_acompcor"
#dir_out <-"301_fc_acompcor"
#dir_in <-"310_ts_aroma"
#dir_out <-"311_fc_aroma"
#dir_in <-"330_ts_acompcor_gsr"
#dir_out <-"331_fc_acompcor_gsr"
#dir_in <-"340_ts_aroma_gsr"
#dir_out <-"341_fc_aroma_gsr"
#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100","schaefer200","schaefer400","shen268")
#list_atlas<-c("aal116","gordon333","power264","shen268")
# list_atlas<-"aal116"

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_CONN"
dir_in <-"56.1_ts_conn"
dir_out <-"56.3_fc_conn"
#list_atlas<-"cnn"
list_atlas<-"hoa"

#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)
library(ggplot2)
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
  path_io     <- file.path(path_root,path_exp_)
  path_in     <- file.path(path_io,dir_in_)
  path_out    <- file.path(path_io,dir_out_)
  output <- list("script"=path_script,"io"=path_io,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# Functional data loading =========================
#**************************************************
#func_data_timeseries<-function(paths,data_clinical,subset_roi,atlas_=atlas){
func_data_timeseries<-function(paths__,atlas_=atlas,key_group="group_3"){
  #print("Starting to load timeseries data.")
  df_timeseries <- read.csv(file.path(paths__$input,"output",paste("atl-",atlas_,"_ts.csv",sep="")))
  #df_timeseries <- df_timeseries[is.element(df_timeseries$ID_pnTTC,
  #                                          data_clinical_$list_id_subj),]
  list_id_roi <- colnames(df_timeseries)[c(-1,-2,-3)]
  dict_roi <- func_dict_roi(paths__)
  dict_roi <- dict_roi[is.element(dict_roi$id,list_id_roi),]
  #dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- sort(as.character(dict_roi$id))
  list_label_roi<-list_group_roi<-NULL
  for (id_roi in list_id_roi){
    list_label_roi<-c(list_label_roi,as.character(dict_roi[dict_roi$id==id_roi,"label"]))
    list_group_roi<-c(list_group_roi,as.character(dict_roi[dict_roi$id==id_roi,"group_3"]))
  }
  n_roi <- length(list_id_roi)
  
  # Create dataframe of ROI-wise BOLD timeseries
  df_timeseries <- cbind(df_timeseries[,c(1,2,3)],
                         df_timeseries[,is.element(colnames(df_timeseries),
                                                   list_id_roi)])
  
  # Create dataframe of ROI group-wise BOLD timeseries by averaging over ROIs within each group
  list_group<-sort(unique(list_group_roi))
  n_group<-length(list_group)
  df_timeseries_group<-df_timeseries[,c(1,2,3)]
  for (group in list_group){
    list_id_roi_group<-list_id_roi[list_group_roi==group]
    df_timeseries_group_add<-df_timeseries[,is.element(colnames(df_timeseries),list_id_roi_group)]
    df_timeseries_group<-cbind(df_timeseries_group,col_mean=rowMeans(df_timeseries_group_add))
    colnames(df_timeseries_group)[colnames(df_timeseries_group)=="col_mean"]<-group
  }
  
  # Create list of IDs and sessions that exist in the data
  list_ses_exist <- sort(unique(df_timeseries$ses))
  list_id_subj_exist<-list()
  n_subj_exist<-NULL
  for (ses in list_ses_exist){
    df_timeseries_ses<-df_timeseries[df_timeseries$ses==ses,]
    list_id_subj_ses<-sort(unique(df_timeseries_ses$ID_pnTTC))
    for (id_subj in list_id_subj_ses){
      list_id_subj_exist<-c(list_id_subj_exist,list(c(ses,id_subj)))
    }
    #list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_timeseries_ses$ID_pnTTC))
    n_subj_exist<-c(n_subj_exist,length(list_id_subj_ses))
  }
  #n_subj_exist <- length(list_id_subj_exist)
  list_ses_exist <- sort(unique(df_timeseries$ses))
  output <- list("df_timeseries"=df_timeseries,"df_timeseries_group"=df_timeseries_group,
                 "list_id_roi"=list_id_roi,"list_label_roi"=list_label_roi,"list_group_roi"=list_group_roi,
                 "list_group"=list_group,
                 "dict_roi"=dict_roi,"n_roi"=n_roi,"n_group"=n_group,
                 "list_id_subj_exist"=list_id_subj_exist,
                 "n_subj_exist"=n_subj_exist,
                 "list_ses_exist"=list_ses_exist)
  #print("Finished loading timeseries data.")
  return(output)
}


#**************************************************
# Functional correlation of each subs =============
#**************************************************

# Core function of functional correlation to be parallelized by fc()
fc_core<-function(data_ts){
  id_subj<-data_ts$subj
  ses<-data_ts$ses
  df_ts<-data_ts$df_ts
  data_fc<-func_cor(input=df_ts)
  
  df_fc_flat<-data_fc$cor_flat
  df_fc_flat<-cbind(ses=ses,ID_pnTTC=id_subj,df_fc_flat)
  colnames(df_fc_flat)<-c("ses","ID_pnTTC","from","to","r","p","z_r")
  file_tmp<-paste("TMP_atl-",atlas,"_sub-",sprintf("%05d", id_subj),"_ses-",sprintf("%02d",ses),suffix_file,".csv",sep="")
  path_file_tmp<-file.path(paths_$output,"output",file_tmp)
  #list_path_tmp<-c(list_path_tmp,path_file_tmp)
  write.csv(df_fc_flat,path_file_tmp,row.names=F)
  #df_fc_stack<-rbind(df_fc_stack,cbind(ses=ses,ID_pnTTC=id_subj,df_fc_flat))
  
  # Convert 'id' to 'label' for heatmap plotting
  df_fc_roilabel<-data_fc$cor
  colnames(df_fc_roilabel)<-rownames(df_fc_roilabel)<-list_label_roi
  
  # Heatmap plot of FC correlation matrix
  plot_fc_heatmap<-plot_cor_heatmap(input=df_fc_roilabel)
  plot_fc_heatmap<-(plot_fc_heatmap
                    + ggtitle(paste("FC, wave",as.character(ses),sprintf("%05d", id_subj),atlas,sep=" "))
                    + theme(plot.title = element_text(hjust = 0.5),
                            axis.title=element_blank()))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_sub-",sprintf("%05d", id_subj),"_ses-",sprintf("%02d",ses),suffix_file,".eps",sep=""),
         plot=plot_fc_heatmap,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
  
  #print(paste("Finished Wave: ",as.character(ses),", Subject: ",as.character(id_subj),sep=""))
  
  return(path_file_tmp)
}


# Main function of functional correlation
fc<-function(paths_=paths,
             list_atlas_=list_atlas){
  print("Starting fc().")
  nullobj<-func_createdirs(paths_,str_proc="fc()")
  
  for (atlas in list_atlas_){
    print(paste("Atlas: ",atlas,", loading data.",sep=""))
    data_timeseries<-func_data_timeseries(paths__=paths_,atlas_=atlas)
    
    # Split timeseries data to each session/subject as input for parallel processing
    list_data_ts<-list_data_ts_group<-list()
    for (ses_subj in data_timeseries$list_id_subj_exist){
      ses<-ses_subj[[1]]
      id_subj<-ses_subj[[2]]
      df_timeseries_ses_subj<-data_timeseries$df_timeseries[which(data_timeseries$df_timeseries$ID_pnTTC==id_subj),]
      df_timeseries_ses_subj<-df_timeseries_ses_subj[which(df_timeseries_ses_subj$ses==ses),c(-1,-2,-3)]
      list_data_ts<-c(list_data_ts,list(list("ses"=ses,"subj"=id_subj,
                                             "df_ts"=df_timeseries_ses_subj)))
      
      df_timeseries_group_ses_subj<-data_timeseries$df_timeseries_group[which(data_timeseries$df_timeseries_group$ID_pnTTC==id_subj),]
      df_timeseries_group_ses_subj<-df_timeseries_group_ses_subj[which(df_timeseries_group_ses_subj$ses==ses),c(-1,-2,-3)]
      list_data_ts_group<-c(list_data_ts_group,list(list("ses"=ses,"subj"=id_subj,
                                                         "df_ts"=df_timeseries_group_ses_subj)))
    }
    
    # Parallel computing of ROI-wise FC for each subject/session
    list_label_roi<-data_timeseries$list_label_roi
    suffix_file="_fc"
    print(paste("Atlas: ",atlas,", calculating ROI-wise FC in parallel.",sep=""))
    clust<-makeCluster(floor(detectCores()*3/4))
    clusterExport(clust,
                  varlist=c("paths_","atlas","list_label_roi","func_cor","FisherZ",
                            "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave","suffix_file"),
                  envir=environment())
    list_path_tmp<-parSapply(clust,list_data_ts,fc_core)
    stopCluster(clust)
    
    # Bind results of ROI-wise FC from temporary files
    print(paste("Atlas: ",atlas,", binding ROI-wise FC results.",sep=""))
    df_fc_stack<-data.frame()
    for (path_tmp in list_path_tmp){
      df_tmp<-read.csv(path_tmp)
      df_fc_stack<-rbind(df_fc_stack,df_tmp)
      file.remove(path_tmp)
    }
    colnames(df_fc_stack)<-c("ses","ID_pnTTC","from","to","r","p","z_r")
    write.csv(df_fc_stack, file.path(paths_$output,"output",paste("atl-",atlas,"_fc.csv",sep="")),row.names = F)
    df_fc_stack<-NULL
    
    # Parallel computing of group-wise FC for each subject/session
    list_label_roi<-data_timeseries$list_group
    suffix_file="_fc_grp"
    print(paste("Atlas: ",atlas,", calculating group-wise FC in parallel.",sep=""))
    clust<-makeCluster(floor(detectCores()*3/4))
    clusterExport(clust,
                  varlist=c("paths_","atlas","list_label_roi","func_cor","FisherZ",
                            "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave","suffix_file"),
                  envir=environment())
    list_path_tmp<-parSapply(clust,list_data_ts_group,fc_core)
    stopCluster(clust)
    
    # Bind results of group-wise FC from temporary files
    print(paste("Atlas: ",atlas,", binding group-wise FC results.",sep=""))
    df_fc_stack<-data.frame()
    for (path_tmp in list_path_tmp){
      df_tmp<-read.csv(path_tmp)
      df_fc_stack<-rbind(df_fc_stack,df_tmp)
      file.remove(path_tmp)
    }
    colnames(df_fc_stack)<-c("ses","ID_pnTTC","from","to","r","p","z_r")
    write.csv(df_fc_stack, file.path(paths_$output,"output",paste("atl-",atlas,"_fc_grp.csv",sep="")),row.names = F)
    df_fc_stack<-NULL

  }
  print("Finished fc().")
}


#**************************************************
# Functional correlation of all subs ==============
#**************************************************
fc_all<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_timeseries<-func_data_timeseries(paths_,data_clinical,subset_roi_)
  data_clinical$list_id_subj<-data_timeseries$list_id_subj_exist
  data_clinical$n_subj<-data_timeseries$n_subj_exist
  nullobj<-func_createdirs(paths_)
  fc<-func_corr(input=data_timeseries$df_timeseries[,c(-1,-2)],
                dict_roi=data_timeseries$dict_roi,
                paths_,prefix_outputfile="FC",
                plot=T,save=T)
  graph<-Corr2Graph(fc)
  fig_circ<-CircularPlot(graph,
                         pvalue_type="p_Benjamini_Hochberg",
                         input_title = "Functional Correlation for All Subjects")
  output<-list("fc"=fc$corr,"fc_flat"=corr$corr_flat,"fig_corrmat"=graph,"fig_circ"=fig_circular)
  return(output)
}