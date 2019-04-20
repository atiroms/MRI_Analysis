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
# Longitudinal clinical data loading ==============
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
func_cor<-function(input){
  cor <-rcorr(as.matrix(input), type="pearson")
  n_node<-ncol(input)
  cor_flat<-data.frame(matrix(nrow=n_node*(n_node-1)/2,ncol=4))
  colnames(cor_flat)<-c("from","to","r","p")
  k<-0
  for (i in 1:(n_node-1)){
    for (j in (i+1):n_node){
      k<-k+1
      cor_flat[k,1:4]<-c(rownames(cor$r)[i],
                          colnames(cor$r)[j],
                          cor$r[i,j],
                          cor$P[i,j])
    }
  }
  output<-list("cor"=cor, "cor_flat"=cor_flat)
  return(output)
}


#**************************************************
# Multiple comparison correction of p values ======
#**************************************************

mltcomp_corr<-function(input){
  output<-data.frame("p_bonferroni"=p.adjust(input$p,method = "bonferroni"),
                     "p_holm_bonferroni"=p.adjust(input$p,method = "holm"),
                     "p_hockberg"=p.adjust(input$p,method = "holm"),
                     "p_hommel"=p.adjust(input$p,method = "hommel"),
                     "p_benjamini_hochberg"=p.adjust(input$p,method="BH"),
                     "p_benjamini_yukutieli"=p.adjust(input$p,method="BY"))
  return(output)
}
