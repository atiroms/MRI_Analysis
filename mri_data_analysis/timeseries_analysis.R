#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for fc()
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
#dir_in   <- "05_ts_temp"
#dir_out  <- "13_fc_temp"
#dir_in   <- "06_ts_t1w"
#dir_out  <- "14_fc_t1w"
#dir_in   <- "07_ts_temponly"
#dir_out  <- "15_fc_temponly"
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
#dir_in   <- "24_ts_acompcor_2mm"
#dir_out  <- "25_fc_acompcor_2mm"

dir_in   <- "27_ts_acompcor"
dir_out  <- "28_fc_acompcor"

#dir_in <-"03_3_ts_acompcor"
#dir_out <-"07_fc_acompcor"

subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_T1QC_new_mild_rsfMRIexist","value"=1))

list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")

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
  path_common <- file.path(path_root,"Dropbox/MRI/pnTTC/Puberty/Stats/CommonData")
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
#func_data_timeseries<-function(paths,data_clinical,subset_roi,atlas_=atlas){
func_data_timeseries<-function(paths__,data_clinical_,atlas_=atlas){
  print("    Starting to load timeseries data.")
  df_timeseries <- read.csv(file.path(paths__$input,"output",paste("timeseries_",atlas_,".csv",sep="")))
  df_timeseries <- df_timeseries[is.element(df_timeseries$ID_pnTTC,
                                            data_clinical_$list_id_subj),]
  list_id_roi <- colnames(df_timeseries)[c(-1,-2)]
  dict_roi <- func_dict_roi(paths__)
  dict_roi <- dict_roi[is.element(dict_roi$id,list_id_roi),]
  #dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- as.character(dict_roi$id)
  n_roi <- length(list_id_roi)
  df_timeseries <- cbind(df_timeseries[,c(1,2)],
                         df_timeseries[,is.element(colnames(df_timeseries),
                                                   list_id_roi)])
  list_id_subj_exist <- sort(unique(df_timeseries$ID_pnTTC))
  n_subj_exist <- length(list_id_subj_exist)
  output <- list("df_timeseries"=df_timeseries,"list_id_roi"=list_id_roi,
                 "dict_roi"=dict_roi,"n_roi"=n_roi,
                 "list_id_subj_exist"=list_id_subj_exist,
                 "n_subj_exist"=n_subj_exist)
  print("    Finished loading timeseries data.")
  return(output)
}


#**************************************************
# Functional correlation of each subs =============
#**************************************************
fc<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi,list_atlas_=list_atlas){
  print("Starting to calculate FC.")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_)
  
  for (atlas in list_atlas_){
    print(paste("  Starting to calculate for atlas: ",atlas,sep=""))
    data_timeseries<-func_data_timeseries(paths__=paths_,data_clinical_=data_clinical,atlas_=atlas)
    data_clinical$list_id_subj_exist<-data_timeseries$list_id_subj_exist
    data_clinical$n_subj_exist<-data_timeseries$n_subj_exist
    
    df_fc_stack<-data.frame()
    for (id_subj in data_clinical$list_id_subj_exist){
      data_fc<-func_cor(input=data_timeseries$df_timeseries[which(data_timeseries$df_timeseries$ID_pnTTC==id_subj),c(-1,-2)])
      
      df_fc_flat<-data_fc$cor_flat
      write.csv(df_fc_flat, file.path(paths_$output,"output",paste(atlas,sprintf("%05d", id_subj),"fc.csv",sep="_")),row.names=F)
      df_fc_stack<-rbind(df_fc_stack,cbind(ID_pnTTC=rep(id_subj,nrow(df_fc_flat)),df_fc_flat))
      
      # Convert ID_long to label_proper for heatmap plotting
      df_fc_roilabel<-data.frame(data_fc$cor$r)
      dict_roi<-data_timeseries$dict_roi
      for(i in seq(ncol(df_fc_roilabel))){
        colnames(df_fc_roilabel)[i]<-as.character(dict_roi[which(dict_roi$id==colnames(df_fc_roilabel)[i]),"label"])
      }
      df_fc_roilabel<-rownames_to_column(df_fc_roilabel, "row")
      for(i in seq(nrow(df_fc_roilabel))){
        df_fc_roilabel$row[i]<-as.character(dict_roi[which(dict_roi$id==df_fc_roilabel$row[i]),"label"])
      }
      
      # Heatmap plot of FC correlation matrix
      fig_fc_heatmap<-cor_heatmap(input=df_fc_roilabel)
      fig_fc_heatmap<-fig_fc_heatmap + ggtitle(paste(sprintf("%05d", id_subj),"Functional Connectivity",sep=" "))+ theme(plot.title = element_text(hjust = 0.5))
      
      # Save heatmap plot
      ggsave(paste(atlas,sprintf("%05d", id_subj),"fc_heatmap.eps",sep="_"),plot=fig_fc_heatmap,device=cairo_ps,
             path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
      
      print(paste("    Finished calculating FC for subject",as.character(id_subj),sep=" "))
    }
    colnames(df_fc_stack)<-c("ID_pnTTC","from","to","r","p")
    print("    Starting to save all subject results.")
    write.csv(df_fc_stack, file.path(paths_$output,"output",paste(atlas,"_fc.csv",sep="")),row.names = F)
    print("    Finished saving all subject results.")
    print(paste("  Finished calculating for atlas: ",atlas, sep=""))
  }
  print("Finished calculating all FCs.")
  return(df_fc_stack)
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