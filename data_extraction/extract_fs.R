#**************************************************
# Description =====================================
#**************************************************
# R script to extract FreeSurfer-preprocessed data


#**************************************************F
# Parameters ======================================
#***************************F***********************

path_in<-c("Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/00_source/w1/13_meas",
           "Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/00_source/w2/18_meas")
path_out<-"Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS/01_extract"


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms"),
                    path_in_=path_in,
                    path_out_=path_out){
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
  path_in     <- file.path(path_root,path_in_)
  path_out    <- file.path(path_root,path_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
  #               "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
                 "common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))


#**************************************************
# GLM of structural measures ======================
#**************************************************

extract_fs<-function(paths_=paths){
  dict_roi <- func_dict_roi(paths_)
  df_out<-data.frame()
  df_out_global<-data.frame()
  for (p_in in paths_$input){
    for (file_in in list.files(file.path(p_in,"output"))){
      path_file_in<-file.path(p_in,"output",file_in)
      df_in<-read.csv(path_file_in,sep="\t")
      colnames_in<-colnames(df_in)[-1]
      if (grepl("w1",file_in,fixed=T)){
        wave<-1
      }else if (grepl("w2",file_in,fixed=T)){
        wave<-2
      }
      for (meas in c("volume","thickness","area")){
        if(grepl(meas,file_in,fixed=T)){
          measure<-meas
        }
        colnames_in<-gsub(paste("_",meas,sep=""),"",colnames_in)
      }
      df_out_sub<-data.frame("ID_pnTTC"=df_in[,1])
      df_out_global_sub<-data.frame("ID_pnTTC"=df_in[,1])
      for (i in seq(length(colnames_in))){
        if (is.element(colnames_in[i],dict_roi$label_fs)){
          df_out_sub<-cbind(df_out_sub,df_in[,i+1])
          colnames(df_out_sub)[ncol(df_out_sub)]<-as.character(dict_roi[which(dict_roi$label_fs==colnames_in[i]),"id"])
        }else{
          df_out_global_sub<-cbind(df_out_global_sub,df_in[,i+1])
          colnames(df_out_global_sub)[ncol(df_out_global_sub)]<-colnames_in[i]
        }
      }
      df_out_sub<-gather(df_out_sub,key=roi,value=value, -ID_pnTTC)
      df_out_sub<-cbind("ID_pnTTC"=df_out_sub$ID_pnTTC,"wave"=wave,"measure"=measure,df_out_sub[,-which(colnames(df_out_sub)=="ID_pnTTC")])
      df_out_global_sub<-gather(df_out_global_sub,key=roi,value=value,-ID_pnTTC)
      df_out_global_sub<-cbind("ID_pnTTC"=df_out_global_sub$ID_pnTTC,"wave"=wave,"measure"="global",df_out_global_sub[,-which(colnames(df_out_global_sub)=="ID_pnTTC")])
      df_out<-rbind(df_out,df_out_sub)
      df_out_global<-rbind(df_out_global,df_out_global_sub)
      print(paste("Finished extracting ",file_in,sep=""))
    }
  }
  df_out_global<-unique(df_out_global)
  df_out<-rbind(df_out,df_out_global)
  write.csv(df_out,file.path(paths_$output,"output","fs_measure.csv"),row.names=F)
  print("Finished saving results.")
  return(df_out)
}