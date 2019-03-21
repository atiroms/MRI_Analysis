#**************************************************
# Description =====================================
#**************************************************
# R script to analyze structural properties from FreeSurfer data.
# Execute DoSCA(), DoPCA(), DoICA(), DoTSNE(), DoJK() to perform each analysis.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS"
dir_in <-"01_extract"
dir_out <-"02_glm"

wave <- 1
measure <-c("volume","thickness","area")
list_covar<-c("W1_Tanner_Max","W1_Age_at_MRI")

#key_global_covar<-"BrainSegVolNotVent"
key_global_covar<-"eTIV"

file_input<-"fs_measure.csv"

subset_subj <- list(list("column"="W1_T1QC","value"=1),list("column"="Sex","value"=1))
#subset_subj <- list(list("column"="W1_T1QC","value"=1),list("column"="Sex","value"=2))
#subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_5sub","value"=1),list("column"="Sex","value"=1))


#**************************************************
# Libraries =======================================
#**************************************************
#library(corrplot)
#library(gplots)
library(Hmisc)
library(FactoMineR)
library(ica)
library(tidyverse)
library(Rtsne)
library(tidyr)
#library(ggpubr)
#library(multcomp)
#library(car)
#library(fastICA)


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
#source(file.path(script_dir,"Functionalities/LI_Functions.R"))


#**************************************************
# GLM of structural measures ======================
#**************************************************
glm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                  wave_=wave,key_global_covar_=key_global_covar
                  ){
  print("Starting glm_str()")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=T)
  df_str<-read.csv(file.path(paths_$input,"output",file_input_))
  df_str$value[which(is.nan(df_str$value))]<-0
  df_global_covar<-df_str[which(df_str$measure=="global" & df_str$wave==wave & df_str$roi==key_global_covar_),]
  for (meas in measure){
    print(paste("    Starting to calculate GLM of ",meas,sep=""))
    df_str_meas<-df_str[which(df_str$measure==meas & df_str$wave==wave_),]
    if (meas=="volume"){
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_,df_global_covar=df_global_covar,key_global_covar=key_global_covar_)
    }else{
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_)
    }
    write.csv(data_glm$glm,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_glm.csv",sep="")),row.names = F)
    write.csv(data_glm$ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$min_ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_min_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$vif,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_vif.csv",sep="")),row.names = F)
    print(paste("    Finished calculating GLM of ",meas,sep=""))
  }
  print("Finished glm_str().")
  
}


#**************************************************
# OBSOLETE ========================================
#**************************************************

##### Parameters ####
#p_uncorrected<-0.001
#p_corrected<-0.05
#
#n_components<-10
##n_components<-30
##n_components<-5
#tsne_dims<-2
#tsne_perplexity<-30
#tsne_max_itr<-1000
#
##### Data Loading ####
#
#source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
#HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
#            colnames(clinical_data)[-1],scale_data = T)
#structural_data<-read.csv(file.path(input_dir,structural_file))
#structural_data$flag<-F
#for (i in subject_id){
#  structural_data[which(structural_data$ID_pnTTC==i),"flag"]<-T
#}
#structural_data<-structural_data[which(structural_data$flag),-(ncol(structural_data))]
#colnames(structural_data)[-1]<-ConvertID(colnames(structural_data)[-1],roi_data,input_roi_type,"ID_long")
#n_rois<-ncol(structural_data)-1
#HeatmapPlot(structural_data,
#            "Scaled ROI Measures",
#            "ROI",
#            ConvertID(colnames(structural_data)[-1],roi_data,"ID_long","label_proper"),
#            scale_data=T)
#
#
##### General Linear Model Analysis ####
#
#DoGLM<-function(){
#  dirname<-ExpDir("GLM")
#  structural_data_tidy<-gather(structural_data,key=ROI,value=value,-ID_pnTTC)
#  structural_data_tidy$ROI_label<-ConvertID(structural_data_tidy$ROI,roi_data,"ID_long","label_proper")
#  glm<-CommonGLM(structural_data_tidy,covariate_label,global_covariate=F,dirname,"GLM_Structure.csv")
#  models_expvars<-glm[which(glm[,"ROI"]==glm[1,"ROI"]),
#                      c("model","exp_var")]
#  glm_ordered<-NULL
#  for (i in 1:nrow(models_expvars)){
#    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
#    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
#    glm_subset<-glm[id_obs,]
#    glm_subset<-MultCompCorr(glm_subset)
#    glm_ordered<-rbind(glm_ordered,glm_subset)
#  }
#  write.csv(glm_ordered, file.path(dirname,"GLM_ordered.csv"),row.names=F)
#  output<-list(glm,glm_ordered)
#  names(output)<-c("GLM","GLM_ordered")
#  return(output)
#}
#
#
##### Structural Correlation Analysis ####
#
#DoSCA<-function(){
#  dirname<-ExpDir("SCA")
#  corr<-CalcCorr(structural_data[-1], dirname,"SCA")
#  graph<-Corr2Graph(corr)
#  fig1<-corr[[3]]
#  fig2<-CircularPlot(graph,
#                     pvalue_type="p_Benjamini_Hochberg",
#                     input_title="Structural Covariance for All Subjects")
#  return(list(corr,fig1,fig2))
#}
#
#
##### Principal Component Analysis ####
#
#DoPCA<-function(){
#  dirname<-ExpDir("PCA")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
#  varfactor<-data.frame(pca$var$coord)
#  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
#  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
#  rownames(varfactor)<-NULL
#  indfactor<-data.frame(pca$ind$coord)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  indfactor<-cbind(indexcolumn,indfactor)
#  varianceaccounted<-pca$eig
#  
#  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
#  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
#  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
#  clincorr<-MeasClinicalCorr(indfactor,dirname)
#  
#  return(list(varfactor,indfactor,varianceaccounted,clincorr))
#}
#
#
##### Independent Component Analysis ####
#
#DoICA<-function(){
#  dirname<-ExpDir("ICA")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
#  varfactor<-data.frame(ica$M)
#  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
#  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
#  rownames(varfactor)<-NULL
#  indfactor<-data.frame(ica$S)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  indfactor<-cbind(indexcolumn,indfactor)
#  varianceaccounted<-ica$vafs
#  
#  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
#  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
#  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
#  clincorr<-MeasClinicalCorr(indfactor,dirname)
#  
#  return(list(varfactor,indfactor,varianceaccounted,clincorr))
#}
#
#
##### t-SNE Analysis ####
#
#DoTSNE<-function(){
#  dirname<-ExpDir("tSNE")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  #  tsne <- Rtsne(data[-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
#  tsne<-Rtsne(data, dims = tsne_dims, perplexity=tsne_perplexity, verbose=TRUE, max_iter = tsne_max_itr)
#  indfactor<-data.frame(tsne$Y)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:tsne_dims)
#  #  pairs(indfactor)
#  indfactor<-cbind(indexcolumn,indfactor)
#  graph<-ComponentPlot(indfactor[-1],"t-SNE")
#  write.csv(indfactor, file.path(dirname,"Coordinates.csv"),row.names=F)
#  return(list(tsne,graph))
#}
#
#
##### Jackknife Estimate of Structural Covariance ####
#
#DoJK<-function(){
#  dirname<-ExpDir("JK")
#  overall_corr<-CalcCorr(structural_data[-1],plot=F,save=F)[[2]]
#  jkstats<-data.frame(matrix(nrow=0,ncol=6))
#  jkstats<-rbind(jkstats,cbind(metric=rep("overall_corr",nrow(overall_corr)),overall_corr[,-6]))
#  jkpe<-jkpv<-jkz<-data.frame(matrix(ncol=6,nrow=0))
#  for (i in 1:n_subject){
#    jk<-CalcCorr(structural_data[-i,-1],plot=F,save=F)[[2]]
#    jkpe<-rbind(jkpe,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,-6]))
#    jkpv<-rbind(jkpv,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,1:4],
#                           n_subject*overall_corr$r-(n_subject-1)*jk$r))
#  }
#  for (j in 1:nrow(overall_corr)){
#    mean_jkpe<-mean(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
#                                              which(jkpe$to==overall_corr[j,"to"])),"r"]))
#    sd_jkpe<-sd(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
#                                          which(jkpe$to==overall_corr[j,"to"])),"r"]))
#    jkstats<-rbind(jkstats,cbind(metric="mean_jkpe",overall_corr[j,1:4],r=mean_jkpe))
#    jkstats<-rbind(jkstats,cbind(metric="sd_jkpe",overall_corr[j,1:4],r=sd_jkpe))
#    for (i in 1:n_subject){
#      jkz<-rbind(jkz,
#                 cbind(ID_pnTTC=subject_id[i], overall_corr[j,1:4],
#                       r=(-1)*(jkpe[intersect(intersect(which(jkpe$ID_pnTTC==subject_id[i]),
#                                                        which(jkpe$from==overall_corr[j,"from"])),
#                                              which(jkpe$to==overall_corr[j,"to"])),"r"]-mean_jkpe)/sd_jkpe))
#    }
#  }
#  jkz<-jkz[order(jkz$ID_pnTTC),]
#  jkstats<-jkstats[c(which(jkstats$metric=="overall_corr"),
#                     which(jkstats$metric=="mean_jkpe"),
#                     which(jkstats$metric=="sd_jkpe")),]
#  colnames(jkstats)<-c("metric","from","from_label","to","to_label","value")
#  colnames(jkpe)<-colnames(jkpv)<-colnames(jkz)<-c("ID_pnTTC","from","from_label","to","to_label","r")
#  write.csv(jkstats, file.path(dirname,"JKStats.csv"),row.names=F)
#  write.csv(jkpe, file.path(dirname,"JKPE.csv"),row.names=F)
#  write.csv(jkpv, file.path(dirname,"JKPV.csv"),row.names=F)
#  write.csv(jkz, file.path(dirname,"JKZ.csv"),row.names=F)
#  return(list(jkstats,jkpe,jkpv,jkz))
#}
#
#
##### Laterality Index Calculation ####
#
#DoLI<-function(){
#  dirname<-ExpDir("LI")
#  structural_data_tidy<-gather(structural_data,key=ROI,value=value,-ID_pnTTC)
#  li<-CommonLI(structural_data_tidy,"ROI",dirname,"LI_Structure.csv")
##  clincorr<-MeasClinicalCorr(output[c(-1,-2),],dirname)
#  li_tidy<-li[,c("ID_pnTTC","ROI","L_ROI_ID","R_ROI_ID","Laterality_Index")]
#  colnames(li_tidy)[5]<-"value"
#  glm<-CommonGLM(li_tidy,covariate_label,F,dirname,"GLM_LI_Structure.csv")
#  models_expvars<-glm[which(glm[,"ROI"]==glm[1,"ROI"]),
#                      c("model","exp_var")]
#  glm_ordered<-NULL
#  for (i in 1:nrow(models_expvars)){
#    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
#    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
#    glm_subset<-glm[id_obs,]
#    glm_subset<-cbind(glm_subset,MultCompCorr(glm_subset))
#    glm_ordered<-rbind(glm_ordered,glm_subset)
#  }
#  write.csv(glm_ordered, file.path(dirname,"GLM_LI_Structure_ordered.csv"),row.names=F)
#  output<-list(li,glm,glm_ordered)
#  names(output)<-c("Laterality_Index","GLM_of_LI","GLM_of_LI_ordered")
#  return(output)
#}
#