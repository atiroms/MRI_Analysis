#**************************************************
# Description =====================================
#**************************************************

# R script to analyze ROI average BOLD signal.
# Execute fc_all() to execute


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
dir_in   <- "03_ts_extract"
dir_out  <- "04_fc"
subset_subj <- list(list("column"="W1_5sub","value"=1))
subset_roi  <- c("Uncertain","Default mode","Sensory/somatomotor Hand","Sensory/somatomotor Mouth",
                 "Fronto-parietal Task Control","Cingulo-opercular Task Control","Subcortical",
                 "Salience","Auditory","Visual","Dorsal attention","Ventral attention",
                 "Memory retrieval?","Cerebellar")


#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro"),
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
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,"common"=path_common)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"Functionalities/Functions.R"))
source(file.path(paths$script,"Functionalities/Graphs.R"))


#**************************************************
# Functional data loading =========================
#**************************************************
func_data_functional<-function(paths,data_clinical,subset_roi){
  df_functional <- read.csv(file.path(paths$input,"output","timeseries.csv"))
  df_functional <- df_functional[is.element(df_functional$ID_pnTTC,data_clinical$list_id_subj),]
  list_id_roi <- colnames(df_functional)[c(-1,-2)]
  dict_roi <- func_dict_roi(paths)
  dict_roi <- dict_roi[is.element(dict_roi$ID_long,list_id_roi),]
  dict_roi <- dict_roi[is.element(dict_roi$group,subset_roi),]
  list_id_roi <- as.character(dict_roi$ID_long)
  n_roi <- length(list_id_roi)
  df_functional <- cbind(df_functional[,c(1,2)],df_functional[,is.element(colnames(df_functional),list_id_roi)])
  output <- list("df_functional"=df_functional,"list_id_roi"=list_id_roi,"dict_roi"=dict_roi,"n_roi"=n_roi)
  return(output)
}


#**************************************************
# Functional correlation of each subs =============
#**************************************************
fc<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_functional<-func_data_functional(paths_,data_clinical,subset_roi_)
  nullobj<-func_createdirs(paths_)
  fc_stack<-data.frame(matrix(ncol=5,nrow=0))
  for (id in data_clinical$list_id_subj){
    fc<-func_corr(input=data_functional$df_functional[which(data_functional$df_functional$ID_pnTTC==id),c(-1,-2)],
                  dict_roi=data_functional$dict_roi,
                  paths_,
                  prefix_outputfile=paste("fc",sprintf("%05d", id),sep="_"),
                  plot=T,save=T)$corr_flat
    fc_stack<-rbind(fc_stack,cbind(ID_pnTTC=rep(id,nrow(fc)),fc))
  }
  colnames(fc_stack)<-c("ID_pnTTC","from","to","r","p")
  write.csv(fc_stack, file.path(paths_$output,"output","fc.csv"),row.names = F)
  return(fc_stack)
}


#**************************************************
# Functional correlation of all subs ==============
#**************************************************
fc_all<-function(paths_=paths,subset_subj_=subset_subj,subset_roi_=subset_roi){
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  data_functional<-func_data_functional(paths_,data_clinical,subset_roi_)
  nullobj<-func_createdirs(paths_)
  fc<-func_corr(input=data_functional$df_functional[,c(-1,-2)],
                dict_roi=data_functional$dict_roi,
                paths_,prefix_outputfile="FC",
                plot=T,save=T)
  graph<-Corr2Graph(fc)
  fig_circ<-CircularPlot(graph,
                         pvalue_type="p_Benjamini_Hochberg",
                         input_title = "Functional Correlation for All Subjects")
  output<-list("fc"=fc$corr,"fc_flat"=corr$corr_flat,"fig_corrmat"=graph,"fig_circ"=fig_circular)
  return(output)
}


#**************************************************
# OBSOLETE ========================================
#**************************************************


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_HO")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_Power")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_DK")
output_dir <- file.path(input_dir,"Functional_data")

#functional_file <- "W1_CONN_BOLD_HO.csv"
#functional_file <- "W1_CONN_BOLD_Power.csv"
functional_file <- "W1_CONN_BOLD_DK.csv"

#for CONN Harvard-Oxford+AAL atlas and FreeSurfer Desikan-Killiany atlas
#roi_subset <- NULL
#roi_subset <- "cortex"
#roi_subset <- "subcortex"
#roi_subset <- "cerebellum"
#roi_subset <- "white matter"
#roi_subset <- "global"
#roi_subset <- "misc"
#roi_subset <- c("cortex","subcortex")
roi_subset <- c("cortex","subcortex","celebellum")
#for Power Atlas
#roi_subset<-c("Uncertain","Default mode","Sensory/somatomotor Hand","Sensory/somatomotor Mouth",
#              "Fronto-parietal Task Control","Cingulo-opercular Task Control","Subcortical",
#              "Salience","Auditory","Visual","Dorsal attention","Ventral attention",
#              "Memory retrieval?","Cerebellar")
#roi_subset<-c("Default mode","Sensory/somatomotor Hand","Sensory/somatomotor Mouth",
#              "Fronto-parietal Task Control","Cingulo-opercular Task Control","Subcortical",
#              "Salience","Auditory","Visual","Dorsal attention","Ventral attention",
#              "Memory retrieval?","Cerebellar")


#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1,Sex=1,W1_Tanner_Stage=1)
subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1)


input_roi_type <- "label_conn"

p_uncorrected<-0.001
p_corrected<-0.05

n_components<-10
#n_components<-4


#### Libraries ####

#library(ggplot2)
library(Hmisc)
#library(corrplot)
#library(gplots)
library(FactoMineR)
library(ica)
library(tidyverse)
#library(ggpubr)


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))
source(file.path(script_dir,"Functionalities/Figures.R"))


#### Data Loading ####
source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
            colnames(clinical_data)[-1],scale_data = T)
functional_data <-read.csv(file.path(input_dir,functional_file))
functional_data$flag<-F
for (i in subject_id){
  functional_data[which(functional_data$ID_pnTTC==i),"flag"]<-T
}
functional_data<-functional_data[which(functional_data$flag),-which(colnames(functional_data)=="flag")]

if (!is.null(roi_subset)){
  roi_group<-ConvertID(colnames(functional_data)[c(-1,-2)],roi_data,"ID_long","group")
  roi_selected<-NULL
  for (i in roi_subset){
    roi_selected<-c(roi_selected,which(roi_group==i))
  }
  functional_data<-functional_data[,c(1,2,roi_selected+2)]
}
n_ROI<-ncol(functional_data)-2


#### Group Functional Correlation ####

DoFC_All<-function(){
  dirname<-ExpDir("FC_All")
  corr<-CalcCorr(functional_data[,c(-1,-2)], dirname,"FC")
  fig1<-corr[[3]]
  graph<-Corr2Graph(corr)
  fig2<-CircularPlot(graph,
                     pvalue_type="p_Benjamini_Hochberg",
                     input_title = "Functional Correlation for All Subjects")
  return(list(corr,fig1,fig2))
}


#### Individual Functional Correlation ####

DoFC<-function(){
  dirname<-ExpDir("FC")
  fc_stack<-data.frame(matrix(ncol=7,nrow=0))
  for (i in 1:n_subject){
    fc<-CalcCorr(functional_data[which(functional_data$ID_pnTTC==subject_id[i]),c(-1,-2)],
                 dirname,sprintf("%05d", subject_id[i]),plot=F,save=F)[[2]]
    fc_stack<-rbind(fc_stack,cbind(ID_pnTTC=rep(subject_id[i],nrow(fc)),fc))
  }
  colnames(fc_stack)<-c("ID_pnTTC","from","from_label","to","to_label","r","p")
  write.csv(fc_stack, file.path(dirname,"FC.csv"),row.names = F)
  return(fc_stack)
}


#### Group Principal Component Analysis ####

DoPCA_All<-function(){
  dirname<-ExpDir("PCA_All")
  data<-functional_data[c(-1,-2)]
  indexcolumn<-functional_data[c(1,2)]
  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
  varfactor<-data.frame(pca$var$coord)
  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_label",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(pca$ind$coord)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-pca$eig
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))

  return(list(varfactor,indfactor,varianceaccounted))
}


#### Group Independent Component Analysis ####

DoICA_All<-function(){
  dirname<-ExpDir("ICA_All")
  data<-functional_data[c(-1,-2)]
  indexcolumn<-functional_data[c(1,2)]
  data<-data.matrix(data)
  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
  varfactor<-data.frame(ica$M)
  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_label",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(ica$S)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-ica$vafs
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  
  return(list(varfactor,indfactor,varianceaccounted))
}


#### Jackknife Estimate One subject at a time ####

# not yet checked after update
DoJK<-function(jktype="JKPE"){
  dirname<-ExpDir(paste(jktype,sep=""))
  for (i in subject_id){
    output<-data.frame(matrix(ncol=7,nrow=0))
    individual_timeseries<-functional_data[which(functional_data$ID_pnTTC==i),]
#    overall_corr<-CalcCorr(individual_timeseries[c(-1,-2)],plot=F,save=F)[[2]]
    overall_corr<-data.frame(ID_pnTTC=rep(i,length(overall_corr)),timeframe=rep(NA,length(overall_corr)),
                           CalcCorr(individual_timeseries[c(-1,-2)],plot=F,save=F)[[2]][,-6])
#    corr<-c(i,NA,overall_corr$r)
#    output<-rbind(output, corr)
    for (j in 1:nrow(individual_timeseries)){
#      print(paste("calculating subject", as.character(i),"time",as.character(j)))
      jk<-CalcCorr(individual_timeseries[-j,c(-1,-2)],plot=F,save=F)[[2]]$r
      if (jktype=="JKPV"){
        jk<-n_subject*overall_corr$r-(n_subject-1)*jk
      }
      jk<-data.frame(ID_pnTTC=rep(i,length(jk)),timeframe=rep(j,length(jk)),
                             from=overall_corr$row,to=overall_corr$column,r=jk)

      output<-rbind(output, overall_corr,jk)
    }
    colnames(output)<-c("ID_pnTTC", "timeframe","from","from_label","to","to_label","r")
#    print("saving")
    write.csv(output, file.path(dirname,sprintf("JK_%05d.csv",i)),row.names=F)
  }
}