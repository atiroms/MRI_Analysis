#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************

# parameters for fc()
path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in <-"53_ts_acompcor"
dir_out <-"54_fc_acompcor"


list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"

#subset_roi  <- c("Uncertain","Default mode","Sensory/somatomotor Hand",
#                 "Sensory/somatomotor Mouth","Fronto-parietal Task Control",
#                 "Cingulo-opercular Task Control","Subcortical","Salience",
#                 "Auditory","Visual","Dorsal attention","Ventral attention",
#                 "Memory retrieval?","Cerebellar")


#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)
library(ggplot2)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms"),
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
  path_common <- file.path(path_root,"Dropbox/MRI/pnTTC/Puberty/Stats/CommonData")
  path_in     <- file.path(path_root,path_exp_,dir_in_)
  path_out    <- file.path(path_root,path_exp_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))
source(file.path(paths$script,"functionality/graph.R"))


#**************************************************
# Functional data loading =========================
#**************************************************
#func_data_timeseries<-function(paths,data_clinical,subset_roi,atlas_=atlas){
func_data_timeseries<-function(paths__,atlas_=atlas){
  print("Starting to load timeseries data.")
  df_timeseries <- read.csv(file.path(paths__$input,"output",paste("atl-",atlas_,"_ts.csv",sep="")))
  #df_timeseries <- df_timeseries[is.element(df_timeseries$ID_pnTTC,
  #                                          data_clinical_$list_id_subj),]
  list_id_roi <- colnames(df_timeseries)[c(-1,-2,-3)]
  dict_roi <- func_dict_roi(paths__)
  dict_roi <- dict_roi[is.element(dict_roi$id,list_id_roi),]
  #dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- as.character(dict_roi$id)
  list_id_roi<-sort(list_id_roi)
  n_roi <- length(list_id_roi)
  df_timeseries <- cbind(df_timeseries[,c(1,2,3)],
                         df_timeseries[,is.element(colnames(df_timeseries),
                                                   list_id_roi)])
  list_ses_exist <- sort(unique(df_timeseries$ses))
  list_id_subj_exist<-list()
  n_subj_exist<-NULL
  for (ses in list_ses_exist){
    df_timeseries_ses<-df_timeseries[df_timeseries$ses==ses,]
    list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_timeseries_ses$ID_pnTTC))
    n_subj_exist<-c(n_subj_exist,length(list_id_subj_exist[[as.character(ses)]]))
  }
  #n_subj_exist <- length(list_id_subj_exist)
  list_ses_exist <- sort(unique(df_timeseries$ses))
  output <- list("df_timeseries"=df_timeseries,"list_id_roi"=list_id_roi,
                 "dict_roi"=dict_roi,"n_roi"=n_roi,
                 "list_id_subj_exist"=list_id_subj_exist,
                 "n_subj_exist"=n_subj_exist,
                 "list_ses_exist"=list_ses_exist)
  print("Finished loading timeseries data.")
  return(output)
}


#**************************************************
# Functional correlation of each subs =============
#**************************************************
fc<-function(paths_=paths,
             list_atlas_=list_atlas){
  print("Starting to calculate FC.")
  nullobj<-func_createdirs(paths_)
  
  for (atlas in list_atlas_){
    print(paste("Starting to calculate for atlas: ",atlas,sep=""))
    data_timeseries<-func_data_timeseries(paths__=paths_,atlas_=atlas)
    
    #df_fc_stack<-data.frame()
    list_path_tmp<-NULL
    for (ses in data_timeseries$list_ses){
      for (id_subj in data_timeseries$list_id_subj_exist[[as.character(ses)]]){
        df_timeseries_ses_subj<-data_timeseries$df_timeseries[which(data_timeseries$df_timeseries$ID_pnTTC==id_subj),]
        df_timeseries_ses_subj<-df_timeseries_ses_subj[which(df_timeseries_ses_subj$ses==ses),c(-1,-2,-3)]
        data_fc<-func_cor(input=df_timeseries_ses_subj)
        
        df_fc_flat<-data_fc$cor_flat
        df_fc_flat<-cbind(ses=ses,ID_pnTTC=id_subj,df_fc_flat)
        colnames(df_fc_flat)<-c("ses","ID_pnTTC","from","to","r","p")
        file_tmp<-paste("TMP_atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d", id_subj),"_fc.csv",sep="")
        path_file_tmp<-file.path(paths_$output,"output",file_tmp)
        list_path_tmp<-c(list_path_tmp,path_file_tmp)
        write.csv(df_fc_flat,path_file_tmp,row.names=F)
        #df_fc_stack<-rbind(df_fc_stack,cbind(ses=ses,ID_pnTTC=id_subj,df_fc_flat))
        
        # Convert 'id' to 'label' for heatmap plotting
        df_fc_roilabel<-data.frame(data_fc$cor$r)
        dict_roi<-data_timeseries$dict_roi
        for(i in seq(ncol(df_fc_roilabel))){
          colnames(df_fc_roilabel)[i]<-as.character(dict_roi[which(dict_roi$id==colnames(df_fc_roilabel)[i]),"label"])
        }
        df_fc_roilabel<-rownames_to_column(df_fc_roilabel, "row")
        for(i in seq(nrow(df_fc_roilabel))){
          df_fc_roilabel$row[i]<-as.character(dict_roi[which(dict_roi$id==df_fc_roilabel$row[i]),"label"])
        }
        
        # Heatmap plot of FC correlation matrix
        fig_fc_heatmap<-cor_heatmap(input=df_fc_roilabel)
        fig_fc_heatmap<-fig_fc_heatmap + ggtitle(paste(sprintf("%05d", id_subj),"Wave",as.character(ses),"Functional Connectivity",sep=" "))+ theme(plot.title = element_text(hjust = 0.5))
        
        # Save heatmap plot
        ggsave(paste("atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d", id_subj),"_fc.eps",sep=""),plot=fig_fc_heatmap,device=cairo_ps,
               path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
        
        print(paste("Finished Wave: ",as.character(ses),", Subject: ",as.character(id_subj),sep=""))
      }
    }
    
    print("Starting to bind all subject results.")
    df_fc_stack<-data.frame()
    for (path_tmp in list_path_tmp){
      df_tmp<-read.csv(path_tmp)
      df_fc_stack<-rbind(df_fc_stack,df_tmp)
      file.remove(path_tmp)
      print(paste("Finished binding: ",path_tmp,sep=""))
    }
    colnames(df_fc_stack)<-c("ses","ID_pnTTC","from","to","r","p")
    write.csv(df_fc_stack, file.path(paths_$output,"output",paste("atl-",atlas,"_fc.csv",sep="")),row.names = F)
    df_fc_stack<-NULL
    print("Finished saving all subject results.")
    print(paste("Finished calculating for atlas: ",atlas, sep=""))
  }
  print("Finished calculating all FCs.")
  return(df_fc_stack)
}


#**************************************************
# Functional correlation of all subs ==============
#**************************************************
fc_all<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_timeseries<-func_data_timeseries(paths_,data_clinical,subset_roi_)
  data_clinical$list_id_subj<-data_timeseries$list_id_subj_exist
  data_clinical$n_subj<-data_timeseries$n_subj_exist
  nullobj<-func_createdirs(paths_)
  fc<-func_corr(input=data_timeseries$df_timeseries[,c(-1,-2)],
                dict_roi=data_timeseries$dict_roi,
                paths_,prefix_outputfile="FC",
                plot=T,save=T)
  graph<-Corr2Graph(fc)
  fig_circ<-CircularPlot(graph,
                         pvalue_type="p_Benjamini_Hochberg",
                         input_title = "Functional Correlation for All Subjects")
  output<-list("fc"=fc$corr,"fc_flat"=corr$corr_flat,"fig_corrmat"=graph,"fig_circ"=fig_circular)
  return(output)
}