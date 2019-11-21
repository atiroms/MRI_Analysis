#**************************************************
# Description =====================================
#**************************************************
# R script to extract XCP-preprocessed functional correlation data for use in GraphVar


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in<-"102_fc_acompcor"
dir_in<-"201_fc_acompcor"
dir_out<-"310_corrmat"

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100","schaefer200","schaefer400","shen268")
list_atlas<-"power264"


#**************************************************
# Libraries =======================================
#**************************************************
library(R.matlab)
library(tidyr)
library(parallel)


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
# Data extraction =================================
#**************************************************

prep_graphvar_core<-function(data_src){
  path_file_out<-data_src$path_file_out
  df_fc<-data_src$df_fc
  corrmat<-matrix(nrow=n_node,ncol=n_node)
  pmat<-matrix(nrow=n_node,ncol=n_node)
  for (i_edge in 1:dim(df_fc)[1]){
    i_x<-which(list_node==df_fc[i_edge,"from"])
    i_y<-which(list_node==df_fc[i_edge,"to"])
    corrmat[i_y,i_x]<-corrmat[i_x,i_y]<-df_fc[i_edge,"r"]
    pmat[i_y,i_x]<-pmat[i_x,i_y]<-df_fc[i_edge,"p"]
    
  }
  for (i_node in 1:n_node){
    corrmat[i_node,i_node]<-1
    pmat[i_node,i_node]<-0
  }
  writeMat(path_file_out,CorrMatrix=corrmat,PValMatrix=pmat)
  return(path_file_out)
}

prep_graphvar<-function(paths_=paths,
                        list_atlas_=list_atlas
                        ){
  print("Starting prep_graphvar()")
  nullobj<-func_createdirs(paths_)
  
  for (atlas in list_atlas_){
    df_fc<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    path_out_atlas<-file.path(paths_$output,"output",atlas)
    if (!file.exists(path_out_atlas)){
      dir.create(path_out_atlas)
    }
    list_node<-sort(unique(c(unique(as.character(df_fc$from)),unique(as.character(df_fc$to)))))
    n_node<-length(list_node)
    list_ses<-sort(unique(df_fc$ses))
    
    list_src_corrmat<-list()
    df_ses_subj<-data.frame()
    for (ses in list_ses){
      list_subj<-sort(unique(df_fc[df_fc$ses==ses,"ID_pnTTC"]))
      for (subj in list_subj){
        df_ses_subj[nrow(df_ses_subj)+1,"Subj_ID"]<-sprintf("%02d_%05d",ses,subj)
        df_fc_subj<-df_fc[df_fc$ses==ses & df_fc$ID_pnTTC==subj,c("from","to","r","p")]
        file_out<-sprintf("%02d_%05d.mat",ses,subj)
        list_src_corrmat<-c(list_src_corrmat,
                            list(list("df_fc"=df_fc_subj,
                                      "path_file_out"=file.path(paths_$output,"output",atlas,file_out))))
      }
    }
    write.csv(df_ses_subj,file.path(paths_$output,"output",atlas,paste("Variables_",atlas,".csv",sep="")),row.names=F)
    
    # Parallel computing of correlation matrices
    print(paste("Atlas: ",atlas,", calculating correlation matrices in parallel.",sep=""))
    clust<-makeCluster(floor(detectCores()*3/4))
    clusterExport(clust,
                  varlist=c("list_node","n_node","writeMat"),
                  envir=environment())
    list_path_tmp<-parSapply(clust,list_src_corrmat,prep_graphvar_core)
    stopCluster(clust)
  }
  print("Finished prep_graphvar()")
}
