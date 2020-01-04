#**************************************************
# Description =====================================
#**************************************************
# R script to extract XCP-preprocessed functional correlation data for use in GraphVar


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

list_id_dir<-list("acompcor"=201,
                  "aroma"=211,
                  "acompcor_gsr"=231,
                  "aroma_gsr"=241,
                  "acompcor"=301,
                  "aroma"=311,
                  "acompcor_gsr"=331,
                  "aroma_gsr"=341)

dir_in<-"201_fc_acompcor"
dir_out<-"400_fc_acompcor"

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
  list_node<-data_src$list_node
  n_node<-data_src$n_node
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
                        list_atlas_=list_atlas,
                        key_grp="group_3"
                        ){
  print("Starting prep_graphvar()")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  
  for (atlas in list_atlas_){
    print(paste("Atlas: ",atlas,", loading FC data.",sep=""))
    path_out_atlas<-file.path(paths_$output,"output",atlas)
    if (!file.exists(path_out_atlas)){
      dir.create(path_out_atlas)
    }
    df_fc<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    list_ses<-sort(unique(df_fc$ses))
    
    # Create list of nodes within each group
    list_node<-sort(unique(c(unique(as.character(df_fc$from)),unique(as.character(df_fc$to)))))
    list_grp<-sort(unique(as.character(dict_roi[dict_roi$id %in% list_node,key_grp])))
    list_node_grp<-list("whole"=list_node)
    dir.create(file.path(path_out_atlas,"whole"))
    for (grp in list_grp){
      list_node_grp_add<-list(sort(unique(as.character(dict_roi[dict_roi$atlas==atlas & dict_roi[key_grp]==grp,"id"]))))
      names(list_node_grp_add)<-grp
      list_node_grp<-c(list_node_grp,list_node_grp_add)
      dir.create(file.path(path_out_atlas,grp))
    }
    
    # Prepare dataset
    list_src_corrmat<-list()
    df_ses_subj<-data.frame()
    for (ses in list_ses){
      list_subj<-sort(unique(df_fc[df_fc$ses==ses,"ID_pnTTC"]))
      for (subj in list_subj){
        #print(paste("Atlas: ",atlas,", Session: ",ses,", Subject: ",subj,", preparing data.",sep=""))
        df_ses_subj[nrow(df_ses_subj)+1,"Subj_ID"]<-sprintf("%02d_%05d",ses,subj)
        df_fc_subj<-df_fc[df_fc$ses==ses & df_fc$ID_pnTTC==subj,c("from","to","r","p")]
        file_out<-sprintf("%02d_%05d.mat",ses,subj)
        for (grp in names(list_node_grp)){
          if (grp=="whole"){
            df_fc_subj_grp<-df_fc_subj
          }else{
            df_fc_subj_grp<-df_fc_subj[(df_fc_subj$from %in% list_node_grp[[grp]]) & (df_fc_subj$to %in% list_node_grp[[grp]]),]
          }
          n_node<-length(list_node_grp[[grp]])
          path_file_out<-file.path(path_out_atlas,grp,file_out)
          list_src_corrmat<-c(list_src_corrmat,
                              list(list("df_fc"=df_fc_subj_grp,"list_node"=list_node_grp[[grp]],"n_node"=n_node,
                                        "path_file_out"=path_file_out)))
        }
      }
    }
    df_ses_subj$dummy<-0
    write.csv(df_ses_subj,file.path(paths_$output,"output",atlas,paste("Variables_",atlas,".csv",sep="")),row.names=F)
    
    # Load inter-group FC data
    print(paste("Atlas: ",atlas,", loading group-wise FC data.",sep=""))
    dir.create(file.path(path_out_atlas,"group"))
    df_fc<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc_grp.csv",sep="")))
    list_ses<-sort(unique(df_fc$ses))
    list_node<-sort(unique(c(unique(as.character(df_fc$from)),unique(as.character(df_fc$to)))))
    n_node<-length(list_node)
    
    # Prepare dataset
    df_ses_subj<-data.frame()
    for (ses in list_ses){
      list_subj<-sort(unique(df_fc[df_fc$ses==ses,"ID_pnTTC"]))
      for (subj in list_subj){
        #print(paste("Atlas: ",atlas,", Session: ",ses,", Subject: ",subj,", preparing data.",sep=""))
        df_ses_subj[nrow(df_ses_subj)+1,"Subj_ID"]<-sprintf("%02d_%05d",ses,subj)
        df_fc_subj<-df_fc[df_fc$ses==ses & df_fc$ID_pnTTC==subj,c("from","to","r","p")]
        file_out<-sprintf("%02d_%05d.mat",ses,subj)
        path_file_out<-file.path(path_out_atlas,"group",file_out)
        list_src_corrmat<-c(list_src_corrmat,
                            list(list("df_fc"=df_fc_subj,"list_node"=list_node,"n_node"=n_node,
                                      "path_file_out"=path_file_out)))
      }
    }
    
    # Parallel computing of correlation matrices
    print(paste("Atlas: ",atlas,", creating correlation matrices in parallel.",sep=""))
    #clust<-makeCluster(floor(detectCores()*3/4))
    clust<-makeCluster(floor(detectCores()*1/4))
    clusterExport(clust,
                  varlist=c("writeMat"),
                  envir=environment())
    list_path_tmp<-parSapply(clust,list_src_corrmat,prep_graphvar_core)
    stopCluster(clust)
  }
  print("Finished prep_graphvar()")
}


#**************************************************
# Data extraction from multiple preprocessing =====
#**************************************************
prep_graphvar_multi<-function(path_exp_=path_exp,
                              list_id_dir_=list_id_dir,list_atlas_=list_atlas){
  
  print("Starting prep_graphvar_multi()")
  for (suffix_dir in names(list_id_dir_)){
    id_dir_fc<-list_id_dir_[[suffix_dir]]
    dir_fc<-paste(as.character(id_dir_fc),"fc",suffix_dir,sep="_")
    dir_fc_mat<-paste(as.character(id_dir_fc),"fc_mat",suffix_dir,sep="_")
    paths_=func_path(path_exp_=path_exp_,
                     dir_in_=dir_fc,
                     dir_out_=dir_fc_mat)
    nullobj<-prep_graphvar(paths_=paths_,list_atlas_=list_atlas_)
  }
  print("Finished prep_graphvar_multi()")
}