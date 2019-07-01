#**************************************************
# Description =====================================
#**************************************************
# R script to analyze fingerprint data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for gta_bin() and gta_weight()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in<-"55_fingerprint"
#dir_out<-"55_gta_bin"
dir_out<-"56_fp_identification"

#list_wave <- c(1,2)

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
list_atlas<-"aal116"
#list_atlas<-"schaefer400"
#list_atlas<-c("glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#thr_pvalue <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************
library(dplyr)
#library(ggplot2)
#library(GGally)
#library(igraph)
#library(qgraph)


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
source(file.path(paths$script,"functionality/function.R"))
#source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
#source(file.path(paths$script,"functionality/gta_function.R"))


#**************************************************
# Fingerprint identification ======================
#**************************************************
fp_identify<-function(paths_=paths,
                      list_atlas_=list_atlas
                      ){
  print("Starting fingerprint identification.")
  nullobj<-func_createdirs(paths_)
  for (atlas in list_atlas_){
    print(paste("Calculate atlas: ",atlas,sep=""))
    file_fp<-paste("atl-",atlas,"_fingerprint.csv",sep="")
    df_fp<-read.csv(file.path(paths_$input,"output",file_fp))
    list_id_subj_exist<-c(df_fp[["from_ID_pnTTC"]],df_fp[["to_ID_pnTTC"]])
    list_id_subj_exist<-unique(list_id_subj_exist)
    list_id_subj_exist<-list_id_subj_exist[order(list_id_subj_exist)]
    list_id_subj_exist_twice<-NULL
    for (id_subj in list_id_subj_exist){
      df_fp_subset_from<-df_fp[df_fp['from_ID_pnTTC']==id_subj,]
      list_subset_ses_from<-unique(df_fp_subset_from[['from_ses']])
      df_fp_subset_to<-df_fp[df_fp['to_ID_pnTTC']==id_subj,]
      list_subset_ses_to<-unique(df_fp_subset_to[['to_ses']])
      list_subset_ses<-unique(c(list_subset_ses_from,list_subset_ses_to))
      if (length(list_subset_ses)==2){
        list_id_subj_exist_twice<-c(list_id_subj_exist_twice,id_subj)
      }
    }
    print(paste("Nuber of subjects with two sessions: ",as.character(length(list_id_subj_exist_twice)),sep=""))
    df_fp_exist_twice<-df_fp[is.element(df_fp$from_ID_pnTTC,list_id_subj_exist_twice),]
    df_fp_exist_twice<-df_fp_exist_twice[is.element(df_fp_exist_twice$to_ID_pnTTC,list_id_subj_exist_twice),]
    df_fp_exist_twice<-df_fp_exist_twice[(df_fp_exist_twice$from_ses==1 & df_fp_exist_twice$to_ses==2),]
    
    for(ses in c(1,2)){
      if(ses==1){
        df_fp_pool<-df_fp_exist_twice[c("from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
        colnames(df_fp_pool)<-c("target_ses","target_ID_pnTTC","pool_ses","pool_ID_pnTTC","r")
      }else{
        df_fp_pool<-df_fp_exist_twice[c("to_ses","to_ID_pnTTC","from_ses","from_ID_pnTTC","r")]
        colnames(df_fp_pool)<-c("target_ses","target_ID_pnTTC","pool_ses","pool_ID_pnTTC","r")
      }
      n_identified<-0
      for (id_subj in list_id_subj_exist_twice){
        df_fp_subset<-df_fp_pool[df_fp_pool$target_ID_pnTTC==id_subj,]
        #df_fp_subset_from<-df_fp_exist_twice[(df_fp_exist_twice$from_ses==ses & df_fp_exist_twice$from_ID_pnTTC==id_subj),]
        #df_fp_subset_from<-df_fp_subset_from[df_fp_subset_from$to_ses==(3-ses),]
        #df_fp_subset_from<-df_fp_subset_from[,c("to_ses","to_ID_pnTTC","r")]
        #df_fp_subset_to<-df_fp_exist_twice[(df_fp_exist_twice$to_ses==ses & df_fp_exist_twice$to_ID_pnTTC==id_subj),]
        #df_fp_subset_to<-df_fp_subset_to[df_fp_subset_to$from_ses==(3-ses),]
        #df_fp_subset_to<-df_fp_subset_to[,c("from_ses","from_ID_pnTTC","r")]
        #colnames(df_fp_subset_to)<-c("to_ses","to_ID_pnTTC","r")
        #df_fp_subset<-rbind(df_fp_subset_from,df_fp_subset_to)
        if(all(is.na(df_fp_subset$r))){
          id_subj_max<-0
        }else{
          id_subj_max<-df_fp_subset[which.max(df_fp_subset$r),"pool_ID_pnTTC"]
        }
        if(id_subj_max==id_subj){
          #print(paste("Subject: ",as.character(id_subj)," identified from Session: ",as.character(3-ses)," using Session: ",as.character(ses)))
          n_identified<-n_identified+1
        }else{
          #print(paste("Subject: ",as.character(id_subj)," not identified from Session: ",as.character(3-ses)," using Session: ",as.character(ses)))
        }
      }
      print(paste("Number of subjects identified: ",as.character(n_identified),sep=""))
      
      for (i in seq(1,1000)){
        df_rand<-data.frame(pool_ID_pnTTC=list_id_subj_exist_twice,rand_ID_pnTTC=sample(list_id_subj_exist_twice,length(list_id_subj_exist_twice)))
        df_fp_rand<-left_join(df_fp_subset,df_rand,by="pool_ID_pnTTC")
        n_identified<-0
        for (id_subj in list_id_subj_exist_twice){
          df_fp_rand<-df_fp_rand[df_fp_rand$target_ID_pnTTC==id_subj,]
          if(all(is.na(df_fp_rand$r))){
            id_subj_max<-0
          }else{
            id_subj_max<-df_fp_rand[which.max(df_fp_rand$r),"rand_ID_pnTTC"]
          }
          if(id_subj_max==id_subj){
            n_identified<-n_identified+1
          }
        }
        print(paste("Iteration: ",as.character(i),", subjects identified: ",as.character(n_identified),sep=""))
      }
      
    }
    
    
  }
}
