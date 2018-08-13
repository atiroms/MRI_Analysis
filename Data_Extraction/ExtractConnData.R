#### Description ####

# R script to extract ROI average BOLD signal from CONN data
# execute ExtractData() for data extraction.
# execute InspectROI() for inspecting ROI labels.


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_HO")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_DK")
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
ROI_id<-4:112   #for surface-based data from freesurfer (DK atlas)
#ROI_id<-168:431    #for Power Atlas in CONN
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