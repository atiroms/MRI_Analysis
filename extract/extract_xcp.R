#**************************************************
# Description =====================================
#**************************************************

# R script to extract ROI average BOLD signal from XCP data
# execute extract_xcp() for data extraction.


#**************************************************
# Parameters ======================================
#**************************************************

#path_in <- "/media/veracrypt2/MRI_img/pnTTC/preproc"
path_in <- "/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc"
#path_in <- "I:/MRI_img/pnTTC/preproc"

#path_out <- "/media/veracrypt2/MRI_img/pnTTC/preproc"
path_out <- "/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc"
#path_out <- "I:/MRI_img/pnTTC/preproc"
#path_out<- '/media/atiroms/MORITA_HDD8/MRI_img/pnTTC/preproc'

#dir_in <-"433_c1_xcp_acompcor"
#dir_out<-"443_c1_ts_acompcor"
#ses<-'ses-01'

#dir_in <-"434_c2_xcp_acompcor"
#dir_out<-"444_c2_ts_acompcor"
#ses<-'ses-02'

#dir_in <-"405_c1_xcp_acompcor_gsr"
#dir_out<-"415_c1_ts_acompcor_gsr"
#ses<-'ses-01'

#dir_in <-"406_c2_xcp_acompcor_gsr"
#dir_out<-"416_c2_ts_acompcor_gsr"
#ses<-'ses-02'

dir_in <-"407_c1_xcp_aroma"
dir_out<-"417_c1_ts_aroma"
ses<-'ses-01'

#dir_in <-"408_c2_xcp_aroma"
#dir_out<-"418_c2_ts_aroma"
#ses<-'ses-02'

#dir_in <-"409_c1_xcp_aroma_gsr"
#dir_out<-"419_c1_ts_aroma_gsr"
#ses<-'ses-01'

#dir_in <-"410_c2_xcp_aroma_gsr"
#dir_out<-"420_c2_ts_aroma_gsr"
#ses<-'ses-02'

#dir_in <-"411_c1_xcp_36p"
#dir_out<-"421_c1_ts_36p"
#ses<-'ses-01'

#dir_in <-"412_c2_xcp_36p"
#dir_out<-"422_c2_ts_36p"
#ses<-'ses-02'

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400","shen268")
#list_atlas<-c("aal116","gordon333","power264","shen268")
list_atlas<-c("aal116","desikanKilliany","glasser360","gordon333","HarvardOxford","power264",
              "schaefer100x7","schaefer100x17","schaefer200x7","schaefer200x17","schaefer400x7","schaefer400x17",
              "shen268")
#list_atlas<-"power264"


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
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
  path_common <- file.path(path_root,"Dropbox/MRI_img/pnTTC/puberty/common")
  path_in     <- file.path(path_in_,dir_in_)
  path_out    <- file.path(path_out_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,"common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Function library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))


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
      if (file.exists(path_input)){
        ts_subj<-read.csv(path_input,header=F,sep=" ")
        ts_subj<-cbind(id_subj, seq(dim(ts_subj)[1]),ts_subj)
        colnames(ts_subj)<- c("ID_pnTTC","timeframe",list_id_roi)
        output<-rbind(output, ts_subj)
        print(paste("Finished extracting atlas: ",atlas, ", subject: ",as.character(id_subj),sep=""))
      }else{
        print(paste("Couled not find atlas: ",atlas, ", subject: ",as.character(id_subj),sep=""))
      }
    }
  }
  print(paste("Starting to save timeseries for atlas: ",atlas,sep=""))
  write.csv(output, file.path(paths__$output,"output",paste("atl-",atlas,"_ts.csv",sep="")),row.names=F)
  print(paste("Finished saving timeseries for atlas: ",atlas,sep=""))
  #print("Finished extracting all files.")
  #return(output)
}

extract_xcp<-function(paths_=paths,
                      list_atlas_=list_atlas
                      ){
  print("Starting extract_xcp().")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  for (atlas in list_atlas_){
    print(paste("Starting to extract XCP results for atlas: ",atlas,sep=""))
    extract_ts_per_atlas(paths__=paths_,atlas=atlas,dict_roi=dict_roi)
    print(paste("Finished extracting XCP results for atlas: ",atlas,sep=""))
  }
  print("Finished extract_xcp().")
}


#**************************************************
# Combine timeseries data from 2 sessions =========
#**************************************************

combine_ts<-function(#path_exp="C:/Users/NICT_WS/Dropbox/temp",
                     path_exp="I:/MRI_img/pnTTC/preproc",
                     
                     #list_src=list(list("dir"="443_c1_ts_acompcor","ses"=1),
                     #              list("dir"="444_c2_ts_acompcor","ses"=2)),
                     #dir_dst="450_ts_acompcor",
                     
                     #list_src=list(list("dir"="415_c1_ts_acompcor_gsr","ses"=1),
                     #              list("dir"="416_c2_ts_acompcor_gsr","ses"=2)),
                     #dir_dst="410_ts_acompcor_gsr",
                     
                     list_src=list(list("dir"="417_c1_ts_aroma","ses"=1),
                                   list("dir"="418_c2_ts_aroma","ses"=2)),
                     dir_dst="420_ts_aroma",
                     
                     #list_src=list(list("dir"="419_c1_ts_aroma_gsr","ses"=1),
                     #              list("dir"="420_c2_ts_aroma_gsr","ses"=2)),
                     #dir_dst="430_ts_aroma_gsr",
                     
                     #list_src=list(list("dir"="421_c1_ts_36p","ses"=1),
                     #              list("dir"="422_c2_ts_36p","ses"=2)),
                     #dir_dst="440_ts_36p",
                     
                     list_atlas_=list_atlas){

  print("Starting combine_ts().")
  for (atlas in list_atlas_){
    print(paste("Calculating atlas: ",atlas,sep=""))
    df_out<-data.frame()
    for (src in list_src){
      dir_src<-src$dir
      path_file_in<-file.path(path_exp,dir_src,"output",paste("atl-",atlas,"_ts.csv",sep=""))
      #print(path_file_in)
      df_out_add<-read.csv(path_file_in)
      df_out_add<-cbind(ses=src$ses,df_out_add)
      df_out<-rbind(df_out,df_out_add)
    }
    file_out<-paste("atl-",atlas,"_ts.csv",sep="")
    path_file_out<-file.path(path_exp,dir_dst,"output",file_out)
    #print(path_file_out)
    write.csv(df_out,path_file_out,row.names=F)
  }
  print("Finished combine_ts().")
}