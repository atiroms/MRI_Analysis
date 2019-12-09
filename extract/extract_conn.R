#**************************************************
# Description =====================================
#**************************************************
# R script to extract CONN-preprocessed timeseries data


#**************************************************
# Parameters ======================================
#***************************F***********************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/fs"

dir_in<-""
dir_out<-""


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


#**************************************************
# Extract CONN results ============================
#**************************************************

extract_conn<-function(paths_=paths,
                       prefix_in="ROI_Subject",
                       suffix_in="_Condition001.mat",
                       folder_conn="project01"){
  print("Starting extract_conn().")
  nullobj<-func_createdirs(paths_)
  dict_roi <- func_dict_roi(paths_)
  
  df_covar<-read.csv(file.path(paths_$input,"output","df_covar.csv"))
  
  df_ts<-data.frame()
  for (id_img in seq(dim(df_covar)[1])){
    id_subj<-as.numeric(df_covar[id_img,"ID_pnTTC"])
    ses<-as.numeric(df_covar[id_img,"ses"])
    
    # Read CONN result in Matlab file
    file_conn<-paste(prefix_in,sprintf("%03d",id_img),suffix_in,sep="")
    data_conn<-readMat(file.path(paths_$input,"output",filder_conn,file_conn))
    
    # Extract ROI names from Matlab file
    list_roi<-NULL
    for (id_roi_conn in seq(length(data_conn$names))){
      list_roi<-c(list_roi,data_conn$names[[i]][[1]][[1]])
    }
    
    # Extract Timeseries data from Matlab file from individual images
    list_ts_roi<-list()
    ts_roi<-data_conn[[ROI_id[i]]][[1]]
    
    df_ts<-data.frame(matrix(ncol=length(ROI_id), nrow=n_timepoint))
    for (i in 1:length(ROI_id)){
      timeseriesdata<-data_conn[[ROI_id[i]]][[1]]
      output[,i]<-timeseriesdata
    }
    timepoint<-1:dim(output)[1]
    output<-cbind(timepoint,output)
    return(output)
    
    # Combine data
    df_ts<-rbind(df_ts,df_ts_add)
    
      
  
    datafilename<-sprintf("%03d", i)
    datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
    matlabdata<-readMat(file.path(input_dir,datafilename))
    individualdata<-ExtractTimeseries(matlabdata$data)
    roinames<-ExtractROINames(matlabdata$names)
    roinames<-roinames[ROI_id,1]
    roinames<-as.character(ConvertID(roinames,roi_data, roi_type,"ID_long"))
    id_pnttc<-as.numeric(ConvertID(i, id_data, id_type, "ID_pnTTC"))
    individualdata<-cbind(id_pnttc, individualdata)
    colnames(individualdata)<- c("ID_pnTTC","timeframe",roinames)
    output<-rbind(output, individualdata)
  }
  
  
  print("Starting extract_conn().")
}
  



#**************************************************
# OBSOLETE ========================================
#**************************************************
#### Description ####

# R script to extract ROI average BOLD signal from CONN data
# execute ExtractData() for data extraction.
# execute InspectROI() for inspecting ROI labels.


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_HO")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_DK")
output_dir <- file.path(input_dir,"Functional_data")

input_fileprefix<-"ROI_Subject"
input_filesuffix<-"_Condition001.mat"

roi_type<-"label_conn"

id_file<-"CSUB.csv"
id_type<-"ID_W1_T1QC_rsfMRIexist"

subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1)
#subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1, "Sex"=1)
#subject_subset <- data.frame("W1_T1QC_rsfMRIexist"=1, "Sex"=2)

subject_id_conn<-1:195
#subject_id_conn<-1:197
#subject_id_conn<-1
#subject_id_conn<-1:5
#ROI_id<-168:177
#ROI_id<-4:112   #for surface-based data from freesurfer (DK atlas)
ROI_id<-168:431    #for Power Atlas in CONN
#ROI_id<-4:135     #for HO and AAL atlas of CONN
#ROI_id<-4:13
n_timepoint<-246


#### Libraries ####

library(R.matlab)


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))


#### Data Loading ####

source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))


#### Data Extraction ####

ExtractTimeseries<-function(input){
  output<-data.frame(matrix(ncol=length(ROI_id), nrow=n_timepoint))
  for (i in 1:length(ROI_id)){
    timeseriesdata<-input[[ROI_id[i]]][[1]]
    output[,i]<-timeseriesdata
  }
  timepoint<-1:dim(output)[1]
  output<-cbind(timepoint,output)
  return(output)
}

ExtractROINames<-function(input){
  allrois<-data.frame(matrix(nrow=length(input),ncol=1))
  for (i in 1:length(input)){
    roiname<-input[[i]][[1]][[1]]
    allrois[i,1]<-roiname
  }
  return(allrois)
}

ExtractData<-function(){
  dirname<-ExpDir("CONN")
  output<-data.frame(matrix(ncol=length(ROI_id)+2, nrow=0))
  for (i in subject_id_conn){
    datafilename<-sprintf("%03d", i)
    datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
    matlabdata<-readMat(file.path(input_dir,datafilename))
    individualdata<-ExtractTimeseries(matlabdata$data)
    roinames<-ExtractROINames(matlabdata$names)
    roinames<-roinames[ROI_id,1]
    roinames<-as.character(ConvertID(roinames,roi_data, roi_type,"ID_long"))
    id_pnttc<-as.numeric(ConvertID(i, id_data, id_type, "ID_pnTTC"))
    individualdata<-cbind(id_pnttc, individualdata)
    colnames(individualdata)<- c("ID_pnTTC","timeframe",roinames)
    output<-rbind(output, individualdata)
  }
  write.csv(output, file.path(dirname,"CONN_data.csv"),row.names=F)
  return(output)
}


#### ROI Name Inspection ####

InspectROI<-function(){
  output<-data.frame((matrix(nrow=1000)))
  for (i in subject_id_conn){
    datafilename<-sprintf("%03d", i)
    datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
    matlabdata<-readMat(file.path(input_dir,datafilename))
    roinames<-ExtractROINames(matlabdata$names)
    roinames<-roinames[,1]
    length(roinames)<-1000
    output<-cbind(output,roinames)
  }
  output<-output[-1]
  return(output)
}