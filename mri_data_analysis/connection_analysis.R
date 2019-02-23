#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for glm_fc()
list_covar<-c("W1_Tanner_Max","W1_Age_at_MRI")

# parameters for fc_corr()
#path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
#dir_in <- c("13_fc_temp","14_fc_t1w","15_fc_temponly","17_fc_36p_2mm",
#            "18_fc_36p_native","19_fc_aroma_2mm","20_fc_acompcor_2mm")
#dir_in <- c("17_fc_36p_2mm","19_fc_aroma_2mm","20_fc_acompcor_2mm")
#dir_out <- "22_fc_corr"
#dir_out <- "23_fc_corr_heatmap"
#dir_in <- "25_fc_acompcor_2mm"
#dir_out <- "26_glm_fc_acompcor"
#dir_out <- "27_glm_fc_acompcor2"
dir_in<-"28_fc_acompcor"
dir_out<-"29_glm_fc_acompcor"

#dir_in <- "07_fc_acompcor"
#dir_out <- "11_glm_fc_acompcor"
subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_5sub","value"=1),list("column"="Sex","value"=1))

list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-c("aal116")

thr_pvalue <- 0.05


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
source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))


#**************************************************
# General linear model of FCs =====================
#**************************************************
glm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                 thr_pvalue_=thr_pvalue,list_atlas_=list_atlas){
  print("Starting glm_fc()")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  for (atlas in list_atlas_){
    print(paste("  Starting to calculate for atlas: ",atlas,sep=""))
    df_fc<-read.csv(file.path(paths_$input,"output",paste(atlas,"_fc.csv",sep="")))
    # Convert NaN's to zero, delete p column and change column name r to value
    df_fc$r[which(is.nan(df_fc$r))]<-0
    df_fc<-df_fc[,-which(colnames(df_fc)=="p")]
    colnames(df_fc)[colnames(df_fc)=="r"]<-"value"
    # Calculate GLM
    print("    Starting to calculate GLM of FCs.")
    data_glm<-func_glm(df_mri=df_fc,data_clinical,list_covar=list_covar_)
    print("    Finished calculating GLM of FCs.")
    print("    Starting to save GLM of FCs.")
    write.csv(data_glm$glm,file.path(paths_$output,"output",paste(atlas,"_glm.csv",sep="")),row.names = F)
    write.csv(data_glm$ic,file.path(paths_$output,"output",paste(atlas,"_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$min_ic,file.path(paths_$output,"output",paste(atlas,"_min_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$vif,file.path(paths_$output,"output",paste(atlas,"_vif.csv",sep="")),row.names = F)
    print("    Finished saving GLM of FCs.")
    
    print("    Starting to calculate seed-level multiple comparison correction and graphs.")
    
    list_roi <- sort(unique(c(as.character(unique(df_fc$from)),as.character(unique(df_fc$to)))))
    dict_roi <- func_dict_roi(paths_)
    dict_roi <- dict_roi[is.element(dict_roi$id,list_roi),]
    
    for (i in seq(length(data_glm$list_model))){
      model<-data_glm$list_model[[i]]
      for (var_exp in model){
        name_model<-paste(list_covar[model],collapse="_")
        name_var_exp<-list_covar[var_exp]
        print(paste("      Starting to calculate for model: ",name_model,", expvar: ",name_var_exp,sep=""))
        df_glm_subset<-data_glm$glm[intersect(which(data_glm$glm[,"model"]==name_model),
                                              which(data_glm$glm[,"var_exp"]==name_var_exp)),]
        # Add columns for global-level multiple comparison-corrected p values
        df_glm_subset<-cbind(df_glm_subset,mltcomp_corr(df_glm_subset))
        
        for (j in list_roi){
          id_obs<-sort(union(which(df_glm_subset$from==j),which(df_glm_subset$to==j)))
          df_glm_subsubset<-df_glm_subset[id_obs,]  # subset of df_glm_subset, starts or ends ad ROI j
          pvalues<-mltcomp_corr(df_glm_subsubset)  # multiple comparison-corrected p values
          for (k in seq(length(id_obs))){  # iterate over connections which starts / ends at ROI j
            for (l in colnames(pvalues)){
              if (is.null(df_glm_subset[id_obs[k],paste("seed",l,sep="_")])){
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
              }else if (is.na(df_glm_subset[id_obs[k],paste("seed",l,sep="_")])){
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
              }else{
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-min(df_glm_subset[id_obs[k],paste("seed",l,sep="_")],
                                                                      pvalues[k,l])
              }
            }
          }
        }
        
        
        # For each model / expvar pair, convert df into nodes/edges data
        graph<-glm_fc2graph(df_glm_subset,list_roi)
        for (j in 1:nrow(graph$node)){
          graph$node[j,"label"]<-as.character(dict_roi[which(dict_roi$id==graph$node[j,"label"]),"label"])
        }
        
        fig_circular<-graph_circular(input=graph,type_pvalue="seed_p_benjamini_hochberg",thr_pvalue=thr_pvalue_)
        fig_circular<-fig_circular +
          ggtitle(paste("GLM Beta\nModel: ",name_model,"\nExplanatory Variable: ",name_var_exp,sep=" ")) +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(paste(atlas,"_model-",name_model,"_expvar-",name_var_exp,"_glm_graph.eps",sep=""),
               plot=fig_circular,device=cairo_ps,path=file.path(paths_$output,"output"),
               dpi=300,height=10,width=10,limitsize=F)
        
        print(paste("      Finished calculating for model: ",name_model,", expvar: ",name_var_exp,sep=""))
      }
    }
    print("    Finished calculating seed-level multiple comparison correction and graphs.")
    print(paste("  Finished calculating for atlas: ",atlas,sep=""))
  }
  
  print("Finished all calculations of glm_fc()")
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

