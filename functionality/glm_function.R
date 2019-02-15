#**************************************************
# Description =====================================
#**************************************************
# R script for performing general linear model fitting calculations


#**************************************************
# Libraries =======================================
#**************************************************

library(multcomp)
library(car)


#**************************************************
# GLM analysis ====================================
#**************************************************

# GLM modelling and statistics calculation per one GLM
glm_routine<-function(df_mri,df_meas_mri,df_covar_sub){
  name_model<-paste(colnames(df_covar_sub)[,-which(colnames(df_covar_sub)=="ID_pnTTC")],sep="_")
  output<-data.frame()
  for (i in seq(nrow(df_meas_mri))){
    df_mri_sub<-df_mri
    for (col_meas in colnames(df_meas_mri)){
      df_mri_sub<-df_mri_sub[which(df_mri_sub[,col_meas]==df_meas_mri[i,col_meas]),]
    }
    
    # Objective variable (some MRI measure)
    var_obj<-as.numeric.factor(df_mri_sub["value"])
    
    # GLM calculation
    formula_lm<-var_obj~1
    for (j in seq(ncol(df_covar_sub))){
      formula_lm<-update(formula_lm,~.+df_covar_sub[,j])
    }
    glmfit<-lm(formula_lm)
    
    # Variance inflation factor
    vifactor<-NaN
    if(ncol(df_covar_sub)>1){
      suppressWarnings(vifactor<-vif(glmfit))
    }
    
    # Calculate stats for each explanatory variables
    for (j in seq(ncol(df_covar_sub))){
      contrast<-matrix(0L,nrow=1, ncol=ncol(df_covar_sub)+1)
      contrast[1,j+1]<-1
      ttest<-tryCatch(summary(glht(glmfit, linfct = contrast))$test,
                      error=function(e){return(NaN)},
                      warning=function(e){return(NaN)},
                      silent=T)
      if (length(ttest)==1){
        output_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j],
                          NaN,NaN,NaN,NaN,
                          vifactor[j],xic[1],xic[2])
      }else{
        output_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j],
                          ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],
                          vifactor[j],xic[1],xic[2])
      }
      colnames(output_meas)<-c(colnames(df_meas_mri),"model","var_exp","beta","sigma","t","p","VIF","AIC","BIC")
      output<-rbind(output,output_add)
    }
  }
  return(output)
}


func_glm<-function(df_mri,data_clinical,list_covar){
  # subset clinical data according to clinical data availability
  df_clinical<-data_clinical$df_clinical
  for (covar in list_covar){
    df_clinical<-df_clinical[which(!is.na(df_clinical[,covar])),]
  }
  list_id_subj<-df_clinical$ID_pnTTC
  
  # Set up covariate dataframe according to subsetting
  df_covar<-df_clinical[,c("IC_pnTTC",list_covar)]
  for (covar in list_covar){
    df_covar[,covar]<-df_covar[,covar]-mean(df_covar[,covar])
  }
  
  # Subset MRI datarame
  df_mri<-df_mri[df_mri$ID_pnTTC %in% list_id_subj]
  
  # df of measurements per subject
  df_meas_mri<-df_mri[df_mri$ID_pnTTC==df_mri[1,ID_pnTTC],]
  df_meas_mri[,c("ID_pnTTC","value")]<-list(NULL)
  
  n_covar<-length(list_covar)
  
  output<-NULL
  for (i in n_covar:1){
    combn_covar<-combn(n_covar,i)
    n_combn_covar<-dim(combn_covar)[2]
    for (j in 1:combn_covar){
      id_covar<-combn_covar[,j]
      df_covar_sub<-df_covar[,c("ID_pnTTC",list_covar[id_covar])]
      output<-rbind(output,glm_routine(df_mri,df_meas_mri,df_covar_sub))
    }
  }
  output$AIC_best<-F
  output$BIC_best<-F
  for (i in 1:n_measures){
    obs_id<-1:nrow(output)
    for (j in colnames(measures)){
      obs_id<-intersect(obs_id,which(output[,j]==measures[i,j]))
    }
    id_minAIC<-which.min(output[obs_id,"AIC"])
    output[obs_id[id_minAIC],"AIC_best"]<-T
    id_minBIC<-which.min(output[obs_id,"BIC"])
    output[obs_id[id_minBIC],"BIC_best"]<-T
  }
  
  #write.csv(output, file.path(dirname,filename),row.names=F)
  return(output)
}
