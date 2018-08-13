#### Description ####

# R script to analyze ROI average BOLD signal.
# Execute DoFC_All(), DoFC(), DoPCA_All(), DoICA_All(), DoJK()


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_HO")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_Power")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Functional_CONN_DK")
output_dir <- file.path(input_dir,"Functional_data")

functional_file <- "W1_CONN_BOLD_HO.csv"
#functional_file <- "W1_CONN_BOLD_Power.csv"
#functional_file <- "W1_CONN_BOLD_DK.csv"

roi_subset<- ""
#roi_subset<- "cortex"
#roi_subset<- "subcortex"
#roi_subset<- "cerebellum"
#roi_subset<- "global"
#roi_subset<- "misc"

subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)


input_roi_type <- "label_conn"

p_uncorrected<-0.001
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
source(file.path(script_dir,"Functionalities/Graphs.R"))


#### Data Loading ####
source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
            colnames(clinical_data)[-1],scale_data = T)
functional_data <-read.csv(file.path(input_dir,functional_file))
functional_data$flag<-F
for (i in subject_id){
  functional_data[which(functional_data$ID_pnTTC==i),"flag"]<-T
}
functional_data<-functional_data[which(functional_data$flag),-ncol(functional_data)]
if (roi_subset!=""){
  functional_data<-cbind(functional_data[,c(1,2)],
                         functional_data[,which(ConvertID(colnames(functional_data)[c(-1,-2)],roi_data,"ID_long","group")==roi_subset)+2])
}
n_ROI<-ncol(functional_data)-2


#### Group Functional Correlation ####

DoFC_All<-function(){
  dirname<-ExpDir("FC_All")
  corr<-CalcCorr(functional_data[,c(-1,-2)], dirname,"FC")
  fig1<-corr[[3]]
  graph<-Corr2Graph(corr)
  fig2<-CircularPlot(graph)
  return(list(corr,fig1,fig2))
}


#### Individual Functional Correlation ####

DoFC<-function(){
  dirname<-ExpDir("FC")
  output<-data.frame(matrix(ncol=7,nrow=0))
  for (i in 1:n_subject){
    fc<-CalcCorr(functional_data[which(functional_data$ID_pnTTC==subject_id[i]),c(-1,-2)],
                 dirname,sprintf("%05d", subject_id[i]),plot=F)[[2]]
    output<-rbind(output,cbind(ID_pnTTC=rep(subject_id[i],nrow(fc)),fc))
  }
  colnames(output)<-c("ID_pnTTC","from","from_label","to","to_label","r","p")
  write.csv(output, file.path(dirname,"FC.csv"),row.names = F)
  return(output)
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