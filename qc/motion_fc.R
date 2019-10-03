#**************************************************
# Description =====================================
#**************************************************
# R script to examine relationships between rsfMRI motion parameter and
# FC divergence from FC average over subjects.


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"201_fc_acompcor"
dir_out<-"205_fc_motion_acompcor"
dir_motion<-c("69_c1_motion","70_c2_motion")

list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"


#**************************************************
# Libraries =======================================
#**************************************************
library(dplyr)
library(ggplot2)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
                    path_exp_=path_exp,
                    dir_in_=dir_in,
                    dir_out_=dir_out,
                    dir_motion_=dir_motion){
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
  path_motion <- file.path(path_root,path_exp_,dir_motion_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,"motion"=path_motion,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_,"dir_motion"=dir_motion_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# FC dist - motion correlation ====================
#**************************************************
fc_motion<-function(paths_=paths,
                    list_atlas_=list_atlas){
  
  print("Starting fc_motion().")
  nullobj<-func_createdirs(paths_)
 
  for (atlas in list_atlas_){
    # Load connection data
    df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    df_edge$from<-as.character(df_edge$from)
    df_edge$to<-as.character(df_edge$to)
    
    # Examine existing subject IDs and sessions in connection data
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_conn_ses$ID_pnTTC))
    }
    # Create combined dataframe of Z-transformed correlation coefficients
    # Also create dataframe of sessions and subjects
    n_edge<-dim(df_edge)[1]
    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
    df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
    colnames(df_ses_subj)<-c("ses","ID_pnTTC")
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        df_conn_subj<-df_conn[df_conn$ID_pnTTC==id_subj & df_conn$ses==ses,"z_r"]
        df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj)
        df_ses_subj<-rbind(df_ses_subj,data.frame(ses=ses,ID_pnTTC=id_subj))
      }
    }
    colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
    rownames(df_conn_cbind)<-NULL
    
    # Calculate average FC
    mean_fc<-rowMeans(df_conn_cbind)
    
    # Calculate distance between average FC and each subject FC
    df_dist_fc<-df_ses_subj
    for (idx_column in seq(dim(df_conn_cbind)[2])){
      df_dist_fc[idx_column,"dist"]<-sqrt(sum((mean_fc-df_conn_cbind[,idx_column])^2))
    }
    
    df_motion<-data.frame()
    for (ses in c(1,2)){
      df_motion_ses<-read_tsv(file.path(paths_$motion[ses],"output","motion.tsv"))
      df_motion_ses<-data.frame(ses=ses,df_motion_ses)
      df_motion<-rbind(df_motion,df_motion_ses)
    }
    
    df_plot<-left_join(df_dist_fc,df_motion,by=c("ses","ID_pnTTC"))
    df_plot$max_trans<-pmax(df_plot$trans_x_max,df_plot$trans_y_max,df_plot$trans_z_max)
    df_plot$max_rot<-pmax(df_plot$rot_x_max,df_plot$rot_y_max,df_plot$rot_z_max)
    df_plot<-df_plot[,c("ses","ID_pnTTC","dist","max_trans","max_rot","max_max")]
    
    for (type_max in c("max_trans","max_rot","max_max")){
      df_plot_subset<-df_plot[,c("dist",type_max)]
      colnames(df_plot_subset)<-c("x","y")
      plot<-(ggplot(df_plot_subset)
             + aes(x=x,y=y)
             + geom_point(size=2)
             + ggtitle(paste("FC-motion correlation,",type_max,sep=" "))
             + xlab("Distance from average")
             + ylab(type_max)
             + theme_light()
             + theme(plot.title = element_text(hjust = 0.5))
            )
      ggsave(paste("atl-",atlas,"_type-",type_max,"_fc_motion.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    }
  }
  print("Finished fc_motion().")
}