#**************************************************
# Description =====================================
#**************************************************
# R script to extract CONN-preprocessed timeseries data


#**************************************************
# Parameters ======================================
#***************************F***********************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/fs"

dir_in<-"D:/MRI_img/pnTTC/c1c2_func/56_conn"
dir_out<-"D:/MRI_img/pnTTC/c1c2_func/56.1_ts_conn"
list_atlas<-c("hoa","cnn","power264")
#list_atlas<-"cnn"


#**************************************************
# Libraries =======================================
#**************************************************
library(R.matlab)
library(rlist)


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
  #path_in     <- file.path(path_root,path_exp_,dir_in_)
  #path_out    <- file.path(path_root,path_exp_,dir_out_)
  path_in<-dir_in_
  path_out<-dir_out_
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()

#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))


#**************************************************
# Extract CONN results ============================
#**************************************************

extract_conn<-function(paths_=paths,
                       prefix_in="ROI_Subject",
                       suffix_in="_Condition001.mat",
                       folder_conn="conn_project01/results/preprocessing",
                       list_atlas_=list_atlas){
  print("Starting extract_conn().")
  nullobj<-func_createdirs(paths_)
  dict_roi <- func_dict_roi(paths_)
  
  df_clin<-read.csv(file.path(paths_$input,"output","df_clin.csv"))
  
  # prepare list of output dataframes
  list_df_ts<-data.frame()
  for (atlas in list_atlas_){
    list_df_ts<-c(list_df_ts,list(data.frame()))
    names(list_df_ts)[length(list_df_ts)]<-atlas
  }
  
  # Iterate over CONN image entries
  for (id_img in seq(dim(df_clin)[1])){
    print(paste("Extracting ",as.character(id_img)," out of ",as.character(dim(df_clin)[1])," image.",sep=""))
    id_subj<-as.numeric(df_clin[id_img,"ID_pnTTC"])
    ses<-as.numeric(df_clin[id_img,"ses"])
    
    # Read CONN result in Matlab file
    file_conn<-paste(prefix_in,sprintf("%03d",id_img),suffix_in,sep="")
    data_conn<-readMat(file.path(paths_$input,"output",folder_conn,file_conn))
    
    # Extract ROI names from Matlab file
    list_roi_conn<-NULL
    for (id_roi_conn in seq(length(data_conn$names))){
      list_roi_conn<-c(list_roi_conn,data_conn$names[[id_roi_conn]][[1]][[1]])
    }
    
    # Prepare atlas-wise list of timeseries for later concatenation
    list_list_ts<-list()
    for (atlas in list_atlas_){
      list_list_ts<-c(list_list_ts,list(list()))
      names(list_list_ts)[length(list_list_ts)]<-atlas
    }
    
    # Extract Timeseries data from Matlab file from individual images
    for (atlas in list_atlas_){
      list_label_conn_atlas<-dict_roi[dict_roi$atlas==atlas,"label_conn"]
      list_id_roi_conn_atlas<-which(is.element(list_roi_conn,list_label_conn_atlas))
      for (id_roi_conn in list_id_roi_conn_atlas){
        label_roi_conn<-list_roi_conn[id_roi_conn]
        id_roi<-as.character(dict_roi[which(dict_roi$label_conn==label_roi_conn),"id"])
        ts_roi<- data_conn$data[[id_roi_conn]][[1]][,1]
        list_list_ts[[atlas]]<-c(list_list_ts[[atlas]],list(ts_roi))
        names(list_list_ts[[atlas]])[length(list_list_ts[[atlas]])]<-id_roi
      }
    }
    
    # Concatenate atlas-wise list of timeseries
    for (atlas in list_atlas_){
      df_ts_atlas<-data.frame(list.cbind(list_list_ts[[atlas]]))
      df_ts_atlas<-data.frame(ses=ses,ID_pnTTC=id_subj,timeframe=seq(dim(df_ts_atlas)[1]),df_ts_atlas)
      list_df_ts[[atlas]]<-rbind(list_df_ts[[atlas]],df_ts_atlas)
    }
  
  }
  
  print("Saving results.")
  for (atlas in list_atlas_){
    write.csv(list_df_ts[[atlas]],file.path(paths_$output,"output",paste("atl-",atlas,"_ts.csv",sep="")),row.names=F)
  }
  
  print("Finished extract_conn().")
}
  



#**************************************************
# OBSOLETE ========================================
#**************************************************
# #### Description ####
# 
# # R script to extract ROI average BOLD signal from CONN data
# # execute ExtractData() for data extraction.
# # execute InspectROI() for inspecting ROI labels.
# 
# 
# #### Parameters ####
# 
# parent_dir <- "D:/atiroms"
# #parent_dir <- "C:/Users/atiro"
# 
# script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
# input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_HO")
# #input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_DK")
# output_dir <- file.path(input_dir,"Functional_data")
# 
# input_fileprefix<-"ROI_Subject"
# input_filesuffix<-"_Condition001.mat"
# 
# roi_type<-"label_conn"
# 
# id_file<-"CSUB.csv"
# id_type<-"ID_W1_T1QC_rsfMRIexist"
# 
# subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1)
# #subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1, "Sex"=1)
# #subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1, "Sex"=2)
# 
# subject_id_conn<-1:195
# #subject_id_conn<-1:197
# #subject_id_conn<-1
# #subject_id_conn<-1:5
# #ROI_id<-168:177
# #ROI_id<-4:112   #for surface-based data from freesurfer (DK atlas)
# ROI_id<-168:431    #for Power Atlas in CONN
# #ROI_id<-4:135     #for HO and AAL atlas of CONN
# #ROI_id<-4:13
# n_timepoint<-246
# 
# 
# #### Libraries ####
# 
# library(R.matlab)
# 
# 
# #### Functionalities ####
# 
# source(file.path(script_dir,"Functionalities/Functions.R"))
# 
# 
# #### Data Loading ####
# 
# source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
# 
# 
# #### Data Extraction ####
# 
# ExtractTimeseries<-function(input){
#   output<-data.frame(matrix(ncol=length(ROI_id), nrow=n_timepoint))
#   for (i in 1:length(ROI_id)){
#     timeseriesdata<-input[[ROI_id[i]]][[1]]
#     output[,i]<-timeseriesdata
#   }
#   timepoint<-1:dim(output)[1]
#   output<-cbind(timepoint,output)
#   return(output)
# }
# 
# ExtractROINames<-function(input){
#   allrois<-data.frame(matrix(nrow=length(input),ncol=1))
#   for (i in 1:length(input)){
#     roiname<-input[[i]][[1]][[1]]
#     allrois[i,1]<-roiname
#   }
#   return(allrois)
# }
# 
# ExtractData<-function(){
#   dirname<-ExpDir("CONN")
#   output<-data.frame(matrix(ncol=length(ROI_id)+2, nrow=0))
#   for (i in subject_id_conn){
#     datafilename<-sprintf("%03d", i)
#     datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
#     matlabdata<-readMat(file.path(input_dir,datafilename))
#     individualdata<-ExtractTimeseries(matlabdata$data)
#     roinames<-ExtractROINames(matlabdata$names)
#     roinames<-roinames[ROI_id,1]
#     roinames<-as.character(ConvertID(roinames,roi_data, roi_type,"ID_long"))
#     id_pnttc<-as.numeric(ConvertID(i, id_data, id_type, "ID_pnTTC"))
#     individualdata<-cbind(id_pnttc, individualdata)
#     colnames(individualdata)<- c("ID_pnTTC","timeframe",roinames)
#     output<-rbind(output, individualdata)
#   }
#   write.csv(output, file.path(dirname,"CONN_data.csv"),row.names=F)
#   return(output)
# }
# 
# 
# #### ROI Name Inspection ####
# 
# InspectROI<-function(){
#   output<-data.frame((matrix(nrow=1000)))
#   for (i in subject_id_conn){
#     datafilename<-sprintf("%03d", i)
#     datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
#     matlabdata<-readMat(file.path(input_dir,datafilename))
#     roinames<-ExtractROINames(matlabdata$names)
#     roinames<-roinames[,1]
#     length(roinames)<-1000
#     output<-cbind(output,roinames)
#   }
#   output<-output[-1]
#   return(output)
# }# 