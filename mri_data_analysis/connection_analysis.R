#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for glm_fc()
list_covar<-c("W1_Tanner_Stage","W1_Age_at_MRI")

# parameters for fc_corr()
path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
dir_in <- c("13_fc_temp","14_fc_t1w","15_fc_temponly","17_fc_36p_2mm",
            "18_fc_36p_native","19_fc_aroma_2mm","20_fc_acompcor_2mm")
#dir_in <- c("17_fc_36p_2mm","19_fc_aroma_2mm","20_fc_acompcor_2mm")
dir_out <- "22_fc_corr"
#dir_out <- "23_fc_corr_heatmap"
subset_subj <- list(list("column"="W1_5sub","value"=1))


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(GGally)


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
  path_common <- file.path(path_root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
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
source(file.path(paths$script,"functionality/graph.R"))


#**************************************************
# General linear model of FCs =====================
#**************************************************
glm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar){
  print("Starting to calculate GLM of FCs.")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=T)
  df_fc<-read.csv(file.path(paths_$input,"output","fc.csv"))
  df_fc<-df_fc[,-which(colnames(df_fc)=="p")]
  colnames(df_fc)[colnames(d_fc)=="r"]<-"value"
  
  df_glm<-func_glm(df_fc,data_clinical,list_covar_)
  path_file_glm<-file.path(paths_$output,"output","glm.csv")
  write.csv(df_glm,path_file_glm)
  print(paste("  GLM of FCs saved in:",path_file_glm,sep=" "))
  print("Finished calculating GLM of FCs.")
  
  # DF of GLM models and explainig variables
  df_model_expvar<-df_glm[intersect(which(df_glm[,"from"]==df_glm[1,"from"]),
                                    which(df_glm[,"to"]==df_glm[1,"to"])),
                          c("model","exp_var")]
  
  fig<-NULL
  glm_ordered<-NULL
  # iterate over model / explaining variable pairs
  for (i in 1:nrow(df_model_expvar)){
    # Subset of df_glm of all connections with the i'th model / expvar pair
    df_glm_subset<-df_glm[intersect(which(glm[,"model"]==df_model_expvar[i,"model"]),
                                    which(glm[,"exp_var"]==df_model_expvar[i,"exp_var"])),]
    # Add columns for all connection-level multiple comparison-corrected p values
    df_glm_subset<-cbind(df_glm_subset,mltcomp_corr(df_glm_subset))
    
    # for each ROI, calculate seed-level multiple comparison-corrected p values
    for (j in rois){
      id_obs<-union(which(glm_subset$from==j),which(glm_subset$to==j))
      id_obs<-id_obs[order(id_obs)]
      glm_subsubset<-glm_subset[id_obs,]  # subset of glm_subset, starts or ends ad ROI j
      pvalues<-MultCompCorr(glm_subsubset)  # multiple comparison-corrected p values
      for (k in 1:length(id_obs)){  # iterate over connections which starts / ends at ROI j
        for (l in colnames(pvalues)){
          if (is.null(glm_subset[id_obs[k],paste("seed",l,sep="_")])){
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
          }else if (is.na(glm_subset[id_obs[k],paste("seed",l,sep="_")])){
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
          }else{
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-min(glm_subset[id_obs[k],paste("seed",l,sep="_")],
                                                               pvalues[k,l])
          }
        }
      }
    }
    nodes_edges<-GLM_FC2Graph(glm_subset,rois)
    fig_title<-paste("GLM Beta of Model:",df_model_expvar[i,"model"],
                     ", Explanatory Variable:",df_model_expvar[i,"exp_var"],sep=" ")
    fig<-c(fig,list(CircularPlot(nodes_edges,
                                 pvalue_type="seed_p_Benjamini_Hochberg",
                                 input_title=fig_title)))
    #    fig<-c(fig,list(CircularPlot(nodes_edges,pvalue_type="p","GLM Beta Values")))
    glm_ordered<-rbind(glm_ordered,glm_subset)
  }
  write.csv(glm_ordered, file.path(dirname,"GLM_FC_ordered.csv"),row.names=F)
  output<-list(glm_ordered,fig)
  names(output)<-c("GLM","Figures")
  return(output)
}



#**************************************************
# FC-FC correlation ===============================
#**************************************************
# for comparison of preprocessing methods

fc_corr<-function(paths_=paths,subset_subj_=subset_subj){
  print("Starting to calculate FC-FC correlation.")
  df_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=F)
  figs<-list()
  for (id_subj in df_clinical$list_id_subj){
    print(paste("Starting to calculate",as.character(id_subj),sep=" "))
    list_path_file_input<-NULL
    id_study_first<-0
    for (id_study in seq(length(paths_$dir_in))){
      file_input<-paste("fc_",sprintf("%05d", id_subj),"_rp.csv",sep="")
      path_file_input<-file.path(paths_$input[id_study],"output",file_input)
      if (file.exists(path_file_input)){
        list_path_file_input<-c(list_path_file_input,path_file_input)
        if(id_study_first==0){
          id_study_first<-id_study
        }
      }else{
        list_path_file_input<-c(list_path_file_input,NA)
      }
    }
    if(id_study_first==0){
      print("No input available.")
    }
    
    for (id_study in seq(length(paths_$dir_in))){
      path_file_input<-list_path_file_input[id_study]
      if(!is.na(path_file_input)){
        df_fc<-read.csv(path_file_input)
        if(id_study==id_study_first){
          df_fc_allstudy<-data.frame(matrix(ncol=length(paths_$dir_in)+2,nrow=nrow(df_fc)))
          colnames(df_fc_allstudy)<-c("from","to",paths_$dir_in)
        }
        df_fc_allstudy[,c("from","to",paths_$dir_in[id_study])]<-df_fc[,c("from","to","r")]
      }
    }
    
    fig<-ggpairs(df_fc_allstudy[,c(-1,-2)],
                 upper=list(continuous=custom_corr_heatmap),
                 #lower=list(continuous=wrap("points",alpha=0.01,size=0.001,stroke = 0, shape = ".")),
                 lower=list(continuous=custom_smooth),
                 diag=list(continuous=custom_densityDiag),
                 title=paste(sprintf("%05d",id_subj),"fc_corr",sep="_"))
    ggsave(paste(sprintf("%05d",id_subj),"fc_corr.eps",sep="_"),plot=fig,device=cairo_ps,
           path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    figs<-c(figs,list(fig))
    print(paste("Finished calculating subject",as.character(id_subj),sep=" "))
  }
  names(figs)<-as.character(df_clinical$list_id_subj)
  print("Finished calculating FC-FC correlation.")
  return(figs)
}

