#**************************************************
# Description =====================================
#**************************************************

# R script for common MRI analysis functions


#**************************************************
# Libraries =======================================
#**************************************************
library(tidyverse)

#**************************************************
# Factor to numeric function ======================
#**************************************************
as.numeric.factor <- function(x) {
  if (class(x)=="factor"){
    return(as.numeric(levels(x))[x])
  }else{
    return(x)
  }
}


#**************************************************
# Experiment folder preparation ===================
#**************************************************
func_createdirs<-function(paths,copy_log=T){
  list_createdirs<-c(paths$output,file.path(paths$output,"output"))
  for(d in list_createdirs){
    if (!file.exists(d)){
      dir.create(d)
    }
  }
  file.copy(file.path(paths$input,"log"),paths$output,recursive=T)
}


#**************************************************
# Returns ROI dictionary ==========================
#**************************************************
func_dict_roi<-function(paths,
                        file_roi="ROI.csv"){
  output<-read.csv(file.path(paths$common,file_roi))
  return(output)
}


#**************************************************
# Clinical data loading ===========================
#**************************************************
func_clinical_data<-function(paths,
                             subset_subj,
                             file_clinical= "CSUB.csv"
                             ){
  df_clinical <- read.csv(file.path(paths$common,file_clinical))
  for (list_subset in subset_subj){
    df_clinical <- df_clinical[which(df_clinical[,list_subset[["column"]]]==list_subset[["value"]]),]
  }
  list_id_subj<-df_clinical$ID_pnTTC
  df_id_subj<-df_clinical[,1:(which(colnames(df_clinical)=="Clinical")-1)]
  df_clinical<-df_clinical[,(-2):(-which(colnames(df_clinical)=="Clinical"))]
  n_subj<-length(list_id_subj)
  n_data_clinical<-ncol(df_clinical)-1
  
  output<-list("df_clinical"=df_clinical,"list_id_subj"=list_id_subj,"df_id_subj"=df_id_subj,
               "n_subj"=n_subj,"n_data_clinical"=n_data_clinical)
  return(output)
}


#**************************************************
# General correlation calculation =================
#**************************************************
func_corr<-function(input, dict_roi, paths, prefix_outputfile, plot=T,save=T,save_plot=T){
  corr <-rcorr(as.matrix(input), type="pearson")
  n_node<-ncol(input)
  corr_flat<-data.frame(matrix(nrow=n_node*(n_node-1)/2,ncol=4))
  colnames(corr_flat)<-c("from","to","r","p")
  k<-0
  for (i in 1:(n_node-1)){
    for (j in (i+1):n_node){
      k<-k+1
      corr_flat[k,1:4]<-c(rownames(corr$r)[i],
                          colnames(corr$r)[j],
                          corr$r[i,j],
                          corr$P[i,j])
    }
  }
  if (plot){
    fig<-plot_corrmat(input=corr$r,dict_roi,title=paste(prefix_outputfile,"correlation matrix"))
    if(save_plot){
      ggsave(paste(prefix_outputfile,"mat.eps",sep="_"),plot=fig,device=cairo_ps,
             path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    }
  }else{
    fig<-NULL
  }
  if (save){
    write.csv(corr_flat, file.path(paths$output,"output",paste(prefix_outputfile,"rp.csv",sep="_")),row.names=F)
  }
  output<-list("corr"=corr, "corr_flat"=corr_flat,"fig"=fig)
  return(output)
}


#**************************************************
# Multiple comparison correction of p values ======
#**************************************************

mltcomp_corr<-function(input){
  output<-data.frame(p_Bonferroni=p.adjust(input$p,method = "bonferroni"))
  output$p_Holm_Bonferroni<-p.adjust(input$p,method = "holm")
  output$p_Hochberg<-p.adjust(input$p,method = "hochberg")
  output$p_Hommel<-p.adjust(input$p,method = "hommel")
  output$p_Benjamini_Hochberg<-p.adjust(input$p,method="BH")
  output$p_Benjamini_Yekutieli<-p.adjust(input$p,method="BY")
  return(output)
}
