#**************************************************
# Description =====================================
#**************************************************

# R script to extract ROI average BOLD signal from XCP data
# execute extract_xcp() for data extraction.


#**************************************************
# Parameters ======================================
#**************************************************
#path_in  <- "P:/MRI/pnTTC/Preproc/test_5sub"
#path_out <- "D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
path_in <- "/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub"
path_out <- "/home/atiroms/Documents/Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
#dir_in   <- "30_xcp_36p"
#dir_out  <- "05_ts_temp"
#dir_in   <- "32_xcp_36p_nativein"
#dir_out  <- "06_ts_t1w"
dir_in   <- "33_xcp_36p_templatein"
dir_out  <- "07_ts_temponly"
#dir_in   <- "37_xcp_36p_spkreg_1mm"
#dir_out  <- "08_ts_36p_1mm"
#dir_in   <- "38_xcp_36p_spkreg_2mm"
#dir_out  <- "09_ts_36p_2mm"
#dir_in   <- "39_xcp_36p_spkreg_native"
#dir_out  <- "10_ts_36p_native"
#dir_in   <- "40_xcp_aroma_2mm"
#dir_out  <- "11_ts_aroma_2mm"
#dir_in   <- "41_xcp_acompcor_2mm"
#dir_out  <- "12_ts_acompcor_2mm"

prefix_file_input<-"sub-"
suffix_file_input<-"_power264_ts.1D"

list_id_subj<-c(14,19,26,28,29)

atlas_roi<-"Power"
list_id_roi<-seq(264)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms/Documents"),
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
  path_common <- file.path(path_root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
  path_in     <- file.path(path_in_,dir_in_)
  path_out    <- file.path(path_out_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,"common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Function library ================================
#**************************************************
source(file.path(paths$script,"Functionalities/Functions.R"))


#**************************************************
# Data extraction =================================
#**************************************************
extract_xcp<-function(paths_=paths,
                      dict_roi_=dict_roi,
                      prefix_file_input_=prefix_file_input,
                      suffix_file_input_=suffix_file_input,
                      list_id_subj_=list_id_subj,
                      atlas_roi_=atlas_roi,
                      list_id_roi_=list_id_roi
                      ){
  
  nullobj<-func_createdirs(paths)
  dict_roi<-func_dict_roi(paths)
  
  output<-data.frame(matrix(ncol=length(list_id_roi_)+2, nrow=0))
  for (id_subj in list_id_subj){
    file_input<-paste(prefix_file_input_, sprintf("%05d", id_subj), suffix_file_input_, sep="")
    path_input<-file.path(paths_$input,"output",paste("sub-",sprintf("%05d", id_subj),sep=""),"fcon","power264",file_input)
    if (!file.exists(path_input)){
      print(paste("File absent for subject", as.character(id_subj),sep=" "))
    }else{
      ts_subj<-read.csv(path_input,header=F,sep=" ")
      df_roi<-dict_roi_[which(dict_roi_$Atlas==atlas_roi_),]
      list_id_roi<-as.character(df_roi$ID_long)
      ts_subj<-cbind(id_subj, seq(dim(ts_subj)[1]),ts_subj)
      colnames(ts_subj)<- c("ID_pnTTC","timeframe",list_id_roi)
      output<-rbind(output, ts_subj)
      print(paste("Finished extracting subject",as.character(id_subj),sep=" "))
    }
  }
  write.csv(output, file.path(paths_$output,"output","timeseries.csv"),row.names=F)
  print("Finished extracting all files.")
  return(output)
}