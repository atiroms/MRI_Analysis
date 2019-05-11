#**************************************************
# Description =====================================
#**************************************************

# R script to extract ROI average BOLD signal from XCP data
# execute extract_xcp() for data extraction.


#**************************************************
# Parameters ======================================
#**************************************************

path_in <- "/media/veracrypt1/MRI/pnTTC/Preproc"
path_out <- "/media/veracrypt2/MRI/pnTTC/Preproc"

dir_in <-"44_c2_1_xcp_acompcor"
dir_out<-"52_c2_1_ts"


list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-c("aal116")

#list_id_subj<-c(14,19,26,28,29)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms"),
                    path_in_=path_in,
                    path_out_=path_out,
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
  path_in     <- file.path(path_in_,dir_in_)
  path_out    <- file.path(path_out_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,"common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Function library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))


#**************************************************
# Data extraction =================================
#**************************************************
extract_ts_per_atlas<-function(paths__,
                                atlas,
                                dict_roi
                                ){
  
  df_roi<-dict_roi[which(dict_roi$atlas==atlas),]
  list_id_roi<-as.character(df_roi$id)
  
  output<-data.frame()
  list_dir_proc<-list.dirs(file.path(paths__$input,"output"),recursive=F)
  for (dir_proc in list_dir_proc){
    list_dir_subj<-list.dirs(dir_proc,recursive=F,full.names=F)
    list_dir_subj<-list_dir_subj[startsWith(list_dir_subj,'sub-')]
    list_id_subj<-as.integer(substring(list_dir_subj,5,9))
    for (id_subj in list_id_subj){
      file_input<-paste("sub-", sprintf("%05d", id_subj), "_",atlas,"_ts.1D", sep="")
      path_input<-file.path(dir_proc,paste("sub-",sprintf("%05d", id_subj),sep=""),"fcon",atlas,file_input)
      ts_subj<-read.csv(path_input,header=F,sep=" ")
      ts_subj<-cbind(id_subj, seq(dim(ts_subj)[1]),ts_subj)
      colnames(ts_subj)<- c("ID_pnTTC","timeframe",list_id_roi)
      output<-rbind(output, ts_subj)
      print(paste("Finished extracting atlas: ",atlas, ", subject: ",as.character(id_subj),sep=""))
    }
  }
  print(paste("Starting to save timeseries for atlas: ",atlas,sep=""))
  write.csv(output, file.path(paths__$output,"output",paste("ts_",atlas,".csv",sep="")),row.names=F)
  print(paste("Finished saving timeseries for atlas: ",atlas,sep=""))
  #print("Finished extracting all files.")
  #return(output)
}

extract_xcp<-function(paths_=paths,
                      list_atlas_=list_atlas
                      ){
  print("Starting to extract XCP results.")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  for (atlas in list_atlas_){
    print(paste("Starting to extract XCP results for atlas: ",atlas,sep=""))
    extract_ts_per_atlas(paths__=paths_,atlas=atlas,dict_roi=dict_roi)
    print(paste("Finished extracting XCP results for atlas: ",atlas,sep=""))
  }
  print("Finished extracting XCP results.")
}