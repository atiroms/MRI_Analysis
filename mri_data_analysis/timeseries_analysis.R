#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for fc()
path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
#dir_in   <- "05_ts_temp"
#dir_out  <- "13_fc_temp"
#dir_in   <- "06_ts_t1w"
#dir_out  <- "14_fc_t1w"
dir_in   <- "07_ts_temponly"
dir_out  <- "15_fc_temponly"
#dir_in   <- "08_ts_36p_1mm"
#dir_out  <- "16_fc_36p_1mm"
#dir_in   <- "09_ts_36p_2mm"
#dir_out  <- "17_fc_36p_2mm"
#dir_in   <- "10_ts_36p_native"
#dir_out  <- "18_fc_36p_native"
#dir_in   <- "11_ts_aroma_2mm"
#dir_out  <- "19_fc_aroma_2mm"
#dir_in   <- "12_ts_acompcor_2mm"
#dir_out  <- "20_fc_acompcor_2mm"
subset_subj <- list(list("column"="W1_5sub","value"=1))
subset_roi  <- c("Uncertain","Default mode","Sensory/somatomotor Hand",
                 "Sensory/somatomotor Mouth","Fronto-parietal Task Control",
                 "Cingulo-opercular Task Control","Subcortical","Salience",
                 "Auditory","Visual","Dorsal attention","Ventral attention",
                 "Memory retrieval?","Cerebellar")


#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)


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
                 "common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))
source(file.path(paths$script,"functionality/graph.R"))


#**************************************************
# Functional data loading =========================
#**************************************************
func_data_timeseries<-function(paths,df_clinical,subset_roi){
  df_timeseries <- read.csv(file.path(paths$input,"output","timeseries.csv"))
  df_timeseries <- df_timeseries[is.element(df_timeseries$ID_pnTTC,
                                            df_clinical$list_id_subj),]
  list_id_roi <- colnames(df_timeseries)[c(-1,-2)]
  dict_roi <- func_dict_roi(paths)
  dict_roi <- dict_roi[is.element(dict_roi$ID_long,list_id_roi),]
  dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- as.character(dict_roi$ID_long)
  n_roi <- length(list_id_roi)
  df_timeseries <- cbind(df_timeseries[,c(1,2)],
                         df_timeseries[,is.element(colnames(df_timeseries),
                                                   list_id_roi)])
  list_id_subj_exists <- sort(unique(df_timeseries$ID_pnTTC))
  n_subj_exists <- length(list_id_subj_exists)
  output <- list("df_timeseries"=df_timeseries,"list_id_roi"=list_id_roi,
                 "dict_roi"=dict_roi,"n_roi"=n_roi,
                 "list_id_subj_exists"=list_id_subj_exists,
                 "n_subj_exists"=n_subj_exists)
  return(output)
}


#**************************************************
# Functional correlation of each subs =============
#**************************************************
fc<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  print("Starting to calculate FC.")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_timeseries<-func_data_timeseries(paths_,data_clinical,subset_roi_)
  data_clinical$list_id_subj<-df_timeseries$list_id_subj_exist
  data_clinical$n_subj_exists<-df_timeseries$n_subj_exists
  nullobj<-func_createdirs(paths_)
  fc_stack<-data.frame(matrix(ncol=5,nrow=0))
  for (id in df_clinical$list_id_subj){
    fc<-func_corr(input=data_timeseries$df_timeseries[which(data_timeseries$df_timeseries$ID_pnTTC==id),c(-1,-2)],
                  dict_roi=data_timeseries$dict_roi,
                  paths_,
                  prefix_outputfile=paste("fc",sprintf("%05d", id),sep="_"),
                  plot=T,save=T)$corr_flat
    fc_stack<-rbind(fc_stack,cbind(ID_pnTTC=rep(id,nrow(fc)),fc))
    print(paste("Finished calculating FC for subject",as.character(id),sep=" "))
  }
  colnames(fc_stack)<-c("ID_pnTTC","from","to","r","p")
  write.csv(fc_stack, file.path(paths_$output,"output","fc.csv"),row.names = F)
  print("Finished calculating all FCs.")
  return(fc_stack)
}


#**************************************************
# Functional correlation of all subs ==============
#**************************************************
fc_all<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_timeseries<-func_data_timeseries(paths_,data_clinical,subset_roi_)
  data_clinical$list_id_subj<-data_timeseries$list_id_subj_exists
  data_clinical$n_subj<-data_timeseries$n_subj_exists
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