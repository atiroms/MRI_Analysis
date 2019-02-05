#### Description ####

# R script to extract ROI average BOLD signal from XCP data
# execute ExtractData() for data extraction.
# execute InspectROI() for inspecting ROI labels.

#### Parameters ####
input_dir <- "02_ts_test_5sub"
output_dir <- "03_ts_extract"

list_parent_dir <- c("D:/atiroms","C:/Users/atiro")

input_fileprefix<-"sub-"
input_filesuffix<-"_power264_ts.1D"

subject_id_xcp<-c(14,19,26,28,29)

roi_type<-"Power"
roi_id<-1:264


#### Directory organization ####

for(dir in list_parent_dir){
  if(file.exists(dir)){
    parent_dir<-dir
  }
}

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/pnTTC/Puberty/Stats/func_XCP",input_dir)
output_dir <- file.path(parent_dir,"DropBox/MRI/pnTTC/Puberty/Stats/func_XCP",output_dir)

list_dircreate<-c(output_dir,file.path(output_dir,"output"))
for(d in list_dircreate){
  if (!file.exists(d)){
    dir.create(d)
  }
}

file.copy(file.path(input_dir,"log"),output_dir,recursive=T)


#### Libraries ####


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))


#### Data Loading ####


#### Data Extraction ####

ExtractData<-function(){
  output<-data.frame(matrix(ncol=length(roi_id)+2, nrow=0))
  for (i in subject_id_xcp){
    datafilename<-sprintf("%05d", i)
    datafilename<-paste(input_fileprefix, datafilename, input_filesuffix, sep="")
    individualdata<-read.csv(file.path(input_dir,"output",datafilename),header=F,sep=" ")
    roi_data<-roi_dict[which(roi_dict$Atlas=="Power"),]
    roinames<-as.character(roi_data$ID_long)
    id_pnttc<-i
    individualdata<-cbind(id_pnttc, seq(dim(individualdata)[1]),individualdata)
    colnames(individualdata)<- c("ID_pnTTC","timeframe",roinames)
    output<-rbind(output, individualdata)
  }
  write.csv(output, file.path(output_dir,"output","timeseries.csv"),row.names=F)
  return(output)
}