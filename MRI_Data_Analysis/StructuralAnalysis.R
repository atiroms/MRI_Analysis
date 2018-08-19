#### Description ####

# R script to analyze structural properties from FreeSurfer data.
# Execute DoSCA(), DoPCA(), DoICA(), DoTSNE(), DoJK() to perform each analysis.


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Structural_FS")
output_dir <- file.path(input_dir,"Structural_data")

#structural_file <- "W1_FS_Volume_Cortex.csv"
structural_file <- "W1_FS_Volume_Subcortex.csv"
#structural_file <- "W1_FS_Volume_WM.csv"
#structural_file <- "W1_FS_Volume_Cerebellum.csv"
#structural_file <- "W1_FS_Volume_Hippocampus.csv"
#structural_file <- "W1_FS_Thickness.csv"
#structural_file <- "W1_FS_Area.csv"
#structural_file <- "W1_FS_Curv.csv"
#structural_file <- "W1_FS_Global.csv"


global_covariate_file<-"W1_FS_Global.csv"
global_covariate_label<-"BrainSegVolNotVent"

#subject_subset <- data.frame(W1_T1QC=1)
subject_subset <- data.frame(W1_T1QC=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)

input_roi_type <- "label_fs"

p_uncorrected<-0.001
n_components<-10
#n_components<-30
#n_components<-5
tsne_dims<-2
tsne_perplexity<-30
tsne_max_itr<-1000


#### Libraries ####

#library(corrplot)
#library(gplots)
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)
library(Rtsne)
#library(ggpubr)
library(multcomp)
library(car)
#library(fastICA)


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))
source(file.path(script_dir,"Functionalities/Figures.R"))


#### Data Loading ####

source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
            colnames(clinical_data)[-1],scale_data = T)
structural_data<-read.csv(file.path(input_dir,structural_file))
structural_data$flag<-F
for (i in subject_id){
  structural_data[which(structural_data$ID_pnTTC==i),"flag"]<-T
}
structural_data<-structural_data[which(structural_data$flag),-(ncol(structural_data))]
colnames(structural_data)[-1]<-ConvertID(colnames(structural_data)[-1],roi_data,input_roi_type,"ID_long")
n_rois<-ncol(structural_data)-1
HeatmapPlot(structural_data,
            "Scaled ROI Measures",
            "ROI",
            ConvertID(colnames(structural_data)[-1],roi_data,"ID_long","label_proper"),
            scale_data=T)


#### General Linear Model Analysis ####

GLMroutine<-function(input_structural_data,input_covar,id_covar,n_expvar){
  output<-data.frame(matrix(ncol=2+5*n_expvar,nrow=n_rois))
  collabel<-colnames(input_covar)[id_covar+1]
  input_covar<-data.frame(input_covar[,id_covar+1])
  colnames(input_covar)<-collabel
  for (i in 1:n_rois){
    if (length(id_covar)==1){
      glmfit<-lm(input_structural_data[,i+1]~input_covar[,1])
    }else if (length(id_covar)==2){
      glmfit<-lm(input_structural_data[,i+1]~input_covar[,1]+input_covar[,2])
    }else if (length(id_covar)==3){
      glmfit<-lm(input_structural_data[,i+1]~input_covar[,1]+input_covar[,2]+input_covar[,3])
    }else if (length(id_covar)==4){
      glmfit<-lm(input_structural_data[,i+1]~input_covar[,1]+input_covar[,2]+input_covar[,3]+input_covar[,4])
    }
    if (length(id_covar)>=2){
      vifactor<-vif(glmfit)
    }else{
      vifactor<-NA
    }
    stats<-c(AIC(glmfit),BIC(glmfit))
    for (j in 1:n_expvar){
      contrast<-matrix(0L,nrow=1, ncol=length(id_covar)+1)
      contrast[1,j+1]<-1
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(stats, ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],vifactor[j])
    }
    output[i,]<-stats
  }
  collabel<-NULL
  for (j in 1:n_expvar){
    collabel<-c(collabel,paste(colnames(input_covar)[j],c("beta","sigma","t","p","VIF"),sep="_"))
  }
  collabel<-c("AIC","BIC",collabel)
  model_name<-paste(colnames(input_covar),collapse="_")
  collabel<-paste(collabel,"of",model_name,"model",sep="_")
  colnames(output)<-collabel
  return(output)
}

DoGLM<-function(covariates_label=c("W1_Tanner_Stage","W1_Age_at_MRI"),global_covariate=F){
  dirname<-ExpDir("GLM")
  n_covariates<-length(covariates_label)
  global_covariate_data<-read.csv(file.path(input_dir,global_covariate_file))
  output<-data.frame(matrix(ncol=2, nrow=n_rois))
  output[,1]<-colnames(structural_data)[-1]
  output[,2]<-ConvertID(colnames(structural_data)[-1],roi_data,"ID_long","label_proper")
  colnames(output)<-c("ROI_ID","ROI_name")
  clinical_data_subset<-clinical_data
  for (i in 1:n_covariates){
    clinical_data_subset<-clinical_data_subset[which(!is.na(clinical_data_subset[,covariates_label[i]])),]
  }
  subject_id_subset<-clinical_data_subset$ID_pnTTC
  
  covariates_data_subset<-data.frame(ID_pnTTC=clinical_data_subset$ID_pnTTC)
  for (i in 1:n_covariates){
    covariates_data_subset<-cbind(covariates_data_subset,clinical_data_subset[,covariates_label[i]])
  }
  for (i in 2:ncol(covariates_data_subset)){
    ave<-mean(covariates_data_subset[,i])
    covariates_data_subset[,i]<-covariates_data_subset[,i]-ave
  }
  colnames(covariates_data_subset)[-1]<-covariates_label
  
  structural_data_subset<-data.frame(matrix(ncol=ncol(structural_data),nrow=0))
  for (i in subject_id_subset){
    structural_data_subset<-rbind(structural_data_subset,structural_data[which(structural_data$ID_pnTTC==i),])
  }
  colnames(structural_data_subset)<-colnames(structural_data)
  
  if (global_covariate){
    global_covariate_data_subset<-data.frame(matrix(ncol=ncol(global_covariate_data),nrow=0))
    for (i in subject_id_subset){
      global_covariate_data_subset<-rbind(global_covariate_data_subset,global_covariate_data[which(global_covariate_data$ID_pnTTC==i),])
    }
    global_covariate_data_subset<-cbind(global_covariate_data_subset$ID_pnTTC,global_covariate_data_subset[,global_covariate_label])
    ave<-mean(global_covariate_data_subset[,2])
    global_covariate_data_subset[,2]<-global_covariate_data_subset[,2]-ave
    all_covariates_data<-cbind(covariates_data_subset,global_covariate_data_subset[,-1])
    colnames(all_covariates_data)[ncol(all_covariates_data)]<-global_covariate_label
  }else{
    all_covariates_data<-covariates_data_subset
  }
  
  for (i in n_covariates:1){
    n_expvar<-i
    for (j in 1:dim(combn(n_covariates,i))[2]){
      id_covar<-combn(n_covariates,i)[,j]
      if (global_covariate){
        id_covar<-c(id_covar,n_covariates+1)
      }
      output<-cbind(output,GLMroutine(structural_data_subset, all_covariates_data,id_covar,n_expvar))
    }
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(structural_data_subset)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in 1:(ncol(structural_data_subset)-1)){
      best_model[j,1]<-which.min(xic[j,])
    }
    colnames(best_model)<-paste(i,"best_model",sep="_")
    output<-cbind(output,best_model)
  }
  write.csv(output, file.path(dirname,"GLM.csv"),row.names=F)
  return(output)
}


#### Structural Correlation Analysis ####

DoSCA<-function(){
  dirname<-ExpDir("SCA")
  corr<-CalcCorr(structural_data[-1], dirname,"SCA")
  graph<-Corr2Graph(corr)
  fig1<-corr[[3]]
  fig2<-CircularPlot(graph)
  return(list(corr,fig1,fig2))
}


#### Principal Component Analysis ####

DoPCA<-function(){
  dirname<-ExpDir("PCA")
  data<-structural_data[-1]
  indexcolumn<-structural_data[1]
  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
  varfactor<-data.frame(pca$var$coord)
  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(pca$ind$coord)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-pca$eig
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  clincorr<-MeasClinicalCorr(indfactor,dirname)
  
  return(list(varfactor,indfactor,varianceaccounted,clincorr))
}


#### Independent Component Analysis ####

DoICA<-function(){
  dirname<-ExpDir("ICA")
  data<-structural_data[-1]
  indexcolumn<-structural_data[1]
  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
  varfactor<-data.frame(ica$M)
  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(ica$S)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-ica$vafs
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  clincorr<-MeasClinicalCorr(indfactor,dirname)
  
  return(list(varfactor,indfactor,varianceaccounted,clincorr))
}


#### t-SNE Analysis ####

DoTSNE<-function(){
  dirname<-ExpDir("tSNE")
  data<-structural_data[-1]
  indexcolumn<-structural_data[1]
  #  tsne <- Rtsne(data[-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
  tsne<-Rtsne(data, dims = tsne_dims, perplexity=tsne_perplexity, verbose=TRUE, max_iter = tsne_max_itr)
  indfactor<-data.frame(tsne$Y)
  colnames(indfactor)<-sprintf("Dim_%02d",1:tsne_dims)
  #  pairs(indfactor)
  indfactor<-cbind(indexcolumn,indfactor)
  graph<-ComponentPlot(indfactor[-1],"t-SNE")
  write.csv(indfactor, file.path(dirname,"Coordinates.csv"),row.names=F)
  return(list(tsne,graph))
}


#### Jackknife Estimate of Structural Covariance ####

DoJK<-function(){
  dirname<-ExpDir("JK")
  overall_corr<-CalcCorr(structural_data[-1],plot=F,save=F)[[2]]
  jkstats<-data.frame(matrix(nrow=0,ncol=6))
  jkstats<-rbind(jkstats,cbind(metric=rep("overall_corr",nrow(overall_corr)),overall_corr[,-6]))
  jkpe<-jkpv<-jkz<-data.frame(matrix(ncol=6,nrow=0))
  for (i in 1:n_subject){
    jk<-CalcCorr(structural_data[-i,-1],plot=F,save=F)[[2]]
    jkpe<-rbind(jkpe,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,-6]))
    jkpv<-rbind(jkpv,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,1:4],
                           n_subject*overall_corr$r-(n_subject-1)*jk$r))
  }
  for (j in 1:nrow(overall_corr)){
    mean_jkpe<-mean(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
                                              which(jkpe$to==overall_corr[j,"to"])),"r"]))
    sd_jkpe<-sd(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
                                          which(jkpe$to==overall_corr[j,"to"])),"r"]))
    jkstats<-rbind(jkstats,cbind(metric="mean_jkpe",overall_corr[j,1:4],r=mean_jkpe))
    jkstats<-rbind(jkstats,cbind(metric="sd_jkpe",overall_corr[j,1:4],r=sd_jkpe))
    for (i in 1:n_subject){
      jkz<-rbind(jkz,
                 cbind(ID_pnTTC=subject_id[i], overall_corr[j,1:4],
                       r=(-1)*(jkpe[intersect(intersect(which(jkpe$ID_pnTTC==subject_id[i]),
                                                        which(jkpe$from==overall_corr[j,"from"])),
                                              which(jkpe$to==overall_corr[j,"to"])),"r"]-mean_jkpe)/sd_jkpe))
    }
  }
  jkz<-jkz[order(jkz$ID_pnTTC),]
  jkstats<-jkstats[c(which(jkstats$metric=="overall_corr"),
                     which(jkstats$metric=="mean_jkpe"),
                     which(jkstats$metric=="sd_jkpe")),]
  colnames(jkstats)<-c("metric","from","from_label","to","to_label","value")
  colnames(jkpe)<-colnames(jkpv)<-colnames(jkz)<-c("ID_pnTTC","from","from_label","to","to_label","r")
  write.csv(jkstats, file.path(dirname,"JKStats.csv"),row.names=F)
  write.csv(jkpe, file.path(dirname,"JKPE.csv"),row.names=F)
  write.csv(jkpv, file.path(dirname,"JKPV.csv"),row.names=F)
  write.csv(jkz, file.path(dirname,"JKZ.csv"),row.names=F)
  return(list(jkstats,jkpe,jkpv,jkz))
}


#### Laterality Index Calculation ####

DoLI<-function(){
  dirname<-ExpDir("LI")
  roi_id<-colnames(structural_data)[-1]
  roi_label<-ConvertID(roi_id,roi_data,"ID_long","label_proper")
  roi_id_left<-roi_id[grep("^L",roi_label)]
  roi_label_left<-roi_label[grep("^L",roi_label)]
  roi_label_bilateral<-substr(roi_label_left,3,1000)
  n_roi_bilateral<-length(roi_label_bilateral)
  roi_label_right<-paste("R",roi_label_bilateral,sep=" ")
  roi_id_right<-NULL
  for (i in 1:n_roi_bilateral){
    roi_id_right<-c(roi_id_right,roi_id[which(roi_label==roi_label_right[i])])
  }
  output<-data.frame(matrix(ncol=n_roi_bilateral+1,nrow=n_subject+2))
  colnames(output)<-c("ID_pnTTC",roi_label_bilateral)
  output[1,]<-c("left_ROI_ID",roi_id_left)
  output[2,]<-c("right_ROI_ID",roi_id_right)
  for (i in 1:n_subject){
    left_measure<-structural_data[which(structural_data$ID_pnTTC==subject_id[i]),roi_id_left]
    right_measure<-structural_data[which(structural_data$ID_pnTTC==subject_id[i]),roi_id_right]
    li<-(left_measure-right_measure)/(left_measure+right_measure)
    output[i+2,]<-c(subject_id[i],li)
  }
  write.csv(output, file.path(dirname,"LateralityIndex.csv"),row.names=F)
  clincorr<-MeasClinicalCorr(output[c(-1,-2),],dirname)
  return(list(output,clincorr))
}
