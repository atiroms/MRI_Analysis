#**************************************************
# Description =====================================
#**************************************************
# R script to extract FreeSurfer-preprocessed data


#**************************************************
# Parameters ======================================
#***************************F***********************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/fs"

dir_in<-"30_meas_fs"
dir_out<-"31_meas"

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
# Extract FS results ==============================
#**************************************************

extract_fs<-function(paths_=paths){
  print("Starting extract_fs().")
  nullobj<-func_createdirs(paths_)
  
  dict_roi <- func_dict_roi(paths_)
  df_out<-data.frame()
  df_out_global<-data.frame()

  for (file_in in list.files(file.path(paths_$input,"output"))){
    path_file_in<-file.path(paths_$input,"output",file_in)
    df_in<-read.csv(path_file_in,sep="\t")
    colnames_in<-colnames(df_in)[-1]
    if (grepl("w1_",file_in,fixed=T)){
      wave<-1
    }else if (grepl("w2_",file_in,fixed=T)){
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
    df_out_sub<-cbind("ses"=wave,"ID_pnTTC"=df_out_sub$ID_pnTTC,"measure"=measure,df_out_sub[,-which(colnames(df_out_sub)=="ID_pnTTC")])
    df_out_global_sub<-gather(df_out_global_sub,key=roi,value=value,-ID_pnTTC)
    df_out_global_sub<-cbind("ses"=wave,"ID_pnTTC"=df_out_global_sub$ID_pnTTC,"measure"="global",df_out_global_sub[,-which(colnames(df_out_global_sub)=="ID_pnTTC")])
    df_out<-rbind(df_out,df_out_sub)
    df_out_global<-rbind(df_out_global,df_out_global_sub)
    print(paste("Finished extracting ",file_in,sep=""))
  }
  
  df_out_global<-unique(df_out_global)
  df_out<-rbind(df_out,df_out_global)
  write.csv(df_out,file.path(paths_$output,"output","fs_measure.csv"),row.names=F)
  print("Finished extract_fs().")
  return(df_out)
}