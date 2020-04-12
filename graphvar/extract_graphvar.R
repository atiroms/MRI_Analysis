#**************************************************
# Description =====================================
#**************************************************
# R script to extract XGraphVar-calculated graph metrices


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/graph_GV"

list_id_dir<-list("acompcor"=201,
                  "aroma"=211,
                  "acompcor_gsr"=231,
                  "aroma_gsr"=241,
                  "acompcor"=301,
                  "aroma"=311,
                  "acompcor_gsr"=331,
                  "aroma_gsr"=341)

#dir_in<-"402.1_graph_acompcor"
#dir_out<-"403.1_graph_acompcor"
dir_in<-"422.2_graph_aroma"
dir_out<-"423.2_graph_aroma"

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100","schaefer200","schaefer400","shen268")
list_atlas<-"power264"
#list_metric_local=c("degrees_und","efficiency_local_bin","rich_club_bu2")
#list_metric_global=c("efficiency_bin","charpath_B_radius",
#                      "charpath_B_diameter","smallworldness_bu")
#list_metric_local=c("degrees_und","efficiency_local_bin")
#list_metric_global=c("efficiency_bin","charpath_B_radius",
#                      "charpath_B_diameter")
#list_metric_local="efficiency_local_bin"
#list_metric_global="smallworldness_bu"
list_metric_local="efficiency_local_bin"
list_metric_global=c("modularity_QOut_und","modularity_louvain_QOut_und","efficiency_bin","small_world_propensity_bin")


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
mean_cost<-function(data_src){
  df_metric_roi<-data_src$df_metric
  roi<-data_src$roi
  df_metric_out<-data.frame()
  for(ses in names(list_subj)){
    for(id_subj in list_subj[[ses]]){
      df_metric_subj<-df_metric_roi[df_metric_roi$ses==ses & df_metric_roi$ID_pnTTC==id_subj,]
      df_metric_add<-df_metric_subj[1,]
      df_metric_add[1,"value"]<-mean(df_metric_subj[,"value"])
      df_metric_add[1,"cost"]<-(-1)
      df_metric_out<-rbind(df_metric_out,df_metric_add)
    }
  }
  
  file_tmp<-paste("TMP_roi-",roi,"_mean.csv",sep="")
  path_file_tmp<-file.path(paths_$output,"output",file_tmp)
  #list_path_tmp<-c(list_path_tmp,path_file_tmp)
  write.csv(df_metric_out,path_file_tmp,row.names=F)
  return(path_file_tmp)
}


extract_gv<-function(paths_=paths,list_atlas_=list_atlas,
                     list_metric_local_=list_metric_local,
                     list_metric_global_=list_metric_global){
  print("Starting extract_graphvar()")
  nullobj<-func_createdirs(paths_)

  for (atlas in list_atlas_){
    path_graph_in<-file.path(paths_$input,"output",atlas)
    list_grp<-list.dirs(path_graph_in,full.names=F,recursive=F)
    
    # Collect local graph variables
    df_metric_local_raw<-data.frame()
    df_metric_local<-data.frame()
    df_metric_global_raw<-data.frame()
    df_metric_global<-data.frame()
    for (grp in list_grp){
      path_grp<-file.path(path_graph_in,grp)
      for (m_local in list_metric_local_){
        print(paste("Atlas: ",atlas,", Group: ",grp,", Metric: ",m_local,", extracting data.",sep=""))
        # Collect within group and metric
        df_metric_grp<-data.frame()
        for (cost in seq(0.1,0.5,0.01)){
          file_in<-paste(m_local,"-",as.character(cost),".txt",sep="")
          path_file_in<-file.path(path_grp,file_in)
          df_local<-read.csv(path_file_in,sep="\t")
          df_metric_add<-gather(df_local,key="roi",value="value",-X)
          df_metric_add$ses<-as.numeric(substr(df_metric_add$X,1,2))
          df_metric_add$ID_pnTTC<-as.numeric(substr(df_metric_add$X,4,8))
          df_metric_add$group<-grp
          df_metric_add$cost<-cost
          df_metric_add$metric<-m_local
          df_metric_add<-df_metric_add[,c("ses","ID_pnTTC","group","roi","metric","cost","value")]
          df_metric_grp<-rbind(df_metric_grp,df_metric_add)
        }
        df_metric_local_raw<-rbind(df_metric_local_raw,df_metric_grp)
        
        # Create list of subjects and list of ROIs using data loaded last
        list_ses<-sort(unique(df_metric_add$ses))
        list_subj<-list()
        for (ses in list_ses){
          list_subj_ses<-list(sort(unique(df_metric_add[df_metric_add$ses==ses,"ID_pnTTC"])))
          names(list_subj_ses)<-as.character(ses)
          list_subj<-c(list_subj,list_subj_ses)
        }
        list_roi<-sort(unique(df_metric_add$roi))
        
        # Calculate average over costs
        list_src<-list()
        print(paste("Atlas: ",atlas,", Group: ",grp,", Metric: ",m_local,", averaging over costs.",sep=""))
        for(roi in list_roi){
          df_metric_roi<-df_metric_grp[df_metric_grp$roi==roi,]
          list_src<-c(list_src,
                      list(list("df_metric"=df_metric_roi,"roi"=roi)))
        }
        clust<-makeCluster(floor(detectCores()*3/4))
        clusterExport(clust,
                      varlist=c("paths_","list_subj","mean"),
                      envir=environment())
        list_path_tmp<-parSapply(clust,list_src,mean_cost)
        stopCluster(clust)
        
        # Bind results of average over costs
        print(paste("Atlas: ",atlas,", Group: ",grp,", Metric: ",m_local,", binding average results.",sep=""))
        df_metric_grp_mean<-data.frame()
        for (path_tmp in list_path_tmp){
          df_tmp<-read.csv(path_tmp)
          df_metric_grp_mean<-rbind(df_metric_grp_mean,df_tmp)
          file.remove(path_tmp)
        }
        df_metric_local<-rbind(df_metric_local,df_metric_grp_mean)
        
        # Calculate average over rois (into global)
        print(paste("Atlas: ",atlas,", Group: ",grp,", Metric: ",m_local,", averaging over ROIs.",sep=""))
        df_metric_grp_mean_global<-data.frame()
        for(ses in names(list_subj)){
          for(id_subj in list_subj[[ses]]){
            #print(paste("Atlas: ",atlas,", Group: ",grp,", Metric: ",m_local,", Ses: ",ses,", Subj: ",
            #            as.character(id_subj),", averaging over ROIs.",sep=""))
            df_metric_subj<-df_metric_grp_mean[df_metric_grp_mean$ses==ses & df_metric_grp_mean$ID_pnTTC==id_subj,]
            df_metric_add<-df_metric_subj[1,]
            df_metric_add[1,"value"]<-mean(df_metric_subj[,"value"])
            df_metric_grp_mean_global<-rbind(df_metric_grp_mean_global,df_metric_add)
          }
        }
        df_metric_grp_mean_global$roi<-"global"
        df_metric_global<-rbind(df_metric_global,df_metric_grp_mean_global)
      }# finish loop over metric
    }# finish loop over group
    
    # Save raw and local metrics
    write.csv(df_metric_local_raw,file.path(paths_$output,"output",paste("atl-",atlas,"_graph_local_raw.csv",sep="")),row.names=F)
    write.csv(df_metric_local,file.path(paths_$output,"output",paste("atl-",atlas,"_graph_local.csv",sep="")),row.names=F)
    
    # Collect global graph variables
    for (grp in list_grp){
      path_grp<-file.path(path_graph_in,grp)
      print(paste("Atlas: ",atlas,", Group: ",grp,", extracting data.",sep=""))
      # Collect within group and metric
      df_metric_grp<-data.frame()
      for (cost in seq(0.1,0.5,0.01)){
        file_in<-paste("globalVariables-",as.character(cost),".txt",sep="")
        path_file_in<-file.path(path_grp,file_in)
        df_global<-read.csv(path_file_in,sep="\t")
        df_metric_add<-df_global[,c("X",list_metric_global_)]
        df_metric_add<-gather(df_metric_add,key="metric",value="value",-X)
        df_metric_add$ses<-as.numeric(substr(df_metric_add$X,1,2))
        df_metric_add$ID_pnTTC<-as.numeric(substr(df_metric_add$X,4,8))
        df_metric_add$group<-grp
        df_metric_add$cost<-cost
        df_metric_add$roi<-"global"
        df_metric_add<-df_metric_add[,c("ses","ID_pnTTC","group","roi","metric","cost","value")]
        df_metric_global_raw<-rbind(df_metric_global_raw,df_metric_add)
      }
    }
    write.csv(df_metric_global_raw,file.path(paths_$output,"output",paste("atl-",atlas,"_graph_global_raw.csv",sep="")),row.names=F)
    
    # Average global metrics over costs
    df_metric_global_mean<-data.frame()
    for (grp in list_grp){
      for (m_global in list_metric_global_){
        df_metric_grp<-df_metric_global_raw[df_metric_global_raw$group==grp
                                            & df_metric_global_raw$metric==m_global,]
        for(ses in names(list_subj)){
          for(id_subj in list_subj[[ses]]){
            df_metric_subj<-df_metric_grp[df_metric_grp$ses==ses
                                          & df_metric_grp$ID_pnTTC==id_subj,]
            df_metric_add<-df_metric_subj[1,]
            df_metric_add$value<-mean(df_metric_subj$value)
            df_metric_add$cost<-(-1)
            df_metric_global_mean<-rbind(df_metric_global_mean,df_metric_add)
          }
        }
      }
    }
    df_metric_global<-rbind(df_metric_global,df_metric_global_mean)
    write.csv(df_metric_global,file.path(paths_$output,"output",paste("atl-",atlas,"_graph_global.csv",sep="")),row.names=F)

  }# finish loop over atlas
  print("Finished extract_graphvar()")
}