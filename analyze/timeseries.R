#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for fc()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in <-"53_ts_acompcor"
#dir_out <-"54_fc_acompcor"
dir_out <-"59_fc_test"


#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
list_atlas<-"aal116"


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
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# Functional data loading =========================
#**************************************************
#func_data_timeseries<-function(paths,data_clinical,subset_roi,atlas_=atlas){
func_data_timeseries<-function(paths__,atlas_=atlas){
  #print("Starting to load timeseries data.")
  df_timeseries <- read.csv(file.path(paths__$input,"output",paste("atl-",atlas_,"_ts.csv",sep="")))
  #df_timeseries <- df_timeseries[is.element(df_timeseries$ID_pnTTC,
  #                                          data_clinical_$list_id_subj),]
  list_id_roi <- colnames(df_timeseries)[c(-1,-2,-3)]
  dict_roi <- func_dict_roi(paths__)
  dict_roi <- dict_roi[is.element(dict_roi$id,list_id_roi),]
  #dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- sort(as.character(dict_roi$id))
  list_label_roi<-NULL
  for (id_roi in list_id_roi){
    list_label_roi<-c(list_label_roi,as.character(dict_roi[dict_roi$id==id_roi,"label"]))
  }
  
  n_roi <- length(list_id_roi)
  df_timeseries <- cbind(df_timeseries[,c(1,2,3)],
                         df_timeseries[,is.element(colnames(df_timeseries),
                                                   list_id_roi)])
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
  output <- list("df_timeseries"=df_timeseries,"list_id_roi"=list_id_roi,"list_label_roi"=list_label_roi,
                 "dict_roi"=dict_roi,"n_roi"=n_roi,
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
  file_tmp<-paste("TMP_atl-",atlas,"_sub-",sprintf("%05d", id_subj),"_ses-",sprintf("%02d",ses),"_fc.csv",sep="")
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
                    + ggtitle(paste("Functional connectivity, wave",as.character(ses),sprintf("%05d", id_subj),atlas,sep=" "))
                    + theme(plot.title = element_text(hjust = 0.5),
                            axis.title=element_blank()))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_sub-",sprintf("%05d", id_subj),"_ses-",sprintf("%02d",ses),"_fc.eps",sep=""),
         plot=plot_fc_heatmap,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
  
  #print(paste("Finished Wave: ",as.character(ses),", Subject: ",as.character(id_subj),sep=""))
  
  return(path_file_tmp)
}


# Main function of functional correlation
fc<-function(paths_=paths,
             list_atlas_=list_atlas){
  print("Starting fc().")
  nullobj<-func_createdirs(paths_)
  
  for (atlas in list_atlas_){
    print(paste("Atlas: ",atlas,", loading data.",sep=""))
    data_timeseries<-func_data_timeseries(paths__=paths_,atlas_=atlas)
    
    # Split timeseries data to each session/subject as input for parallel processing
    list_data_ts<-list()
    for (ses_subj in data_timeseries$list_id_subj_exist){
      ses<-ses_subj[[1]]
      id_subj<-ses_subj[[2]]
      df_timeseries_ses_subj<-data_timeseries$df_timeseries[which(data_timeseries$df_timeseries$ID_pnTTC==id_subj),]
      df_timeseries_ses_subj<-df_timeseries_ses_subj[which(df_timeseries_ses_subj$ses==ses),c(-1,-2,-3)]
      list_data_ts<-c(list_data_ts,list(list("ses"=ses,"subj"=id_subj,"df_ts"=df_timeseries_ses_subj)))
    }
    list_label_roi<-data_timeseries$list_label_roi
    
    # Parallel computing of fc for each subject/session
    print(paste("Atlas: ",atlas,", calculating FC in parallel.",sep=""))
    clust<-makeCluster(floor(detectCores()*3/4))
    clusterExport(clust,
                  varlist=c("paths_","atlas","list_label_roi","func_cor",
                            "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave"),
                  envir=environment())
    list_path_tmp<-parSapply(clust,list_data_ts,fc_core)
    stopCluster(clust)
    
    # Bind results in temporary files
    print(paste("Atlas: ",atlas,", binding results."))
    df_fc_stack<-data.frame()
    for (path_tmp in list_path_tmp){
      df_tmp<-read.csv(path_tmp)
      df_fc_stack<-rbind(df_fc_stack,df_tmp)
      file.remove(path_tmp)
      #print(paste("Finished binding: ",path_tmp,sep=""))
    }
    colnames(df_fc_stack)<-c("ses","ID_pnTTC","from","to","r","p","z_r")
    write.csv(df_fc_stack, file.path(paths_$output,"output",paste("atl-",atlas,"_fc.csv",sep="")),row.names = F)
    df_fc_stack<-NULL
    #print("Finished saving all subject results.")
    #print(paste("Finished calculating for atlas: ",atlas, sep=""))
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