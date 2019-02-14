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
#glm_routine<-function(input_MRI_data,input_measures,input_covar,id_covar,n_expvar){
glm_routine<-function(df_mri,df_meas_mri,df_covar_sub){
  name_model<-paste(colnames(df_covar_sub)[,-which(colnames(df_covar_sub)=="ID_pnTTC")],sep="_")
  for (i in seq(nrow(df_meas_mri))){
    df_mri_sub<-df_mri
    for (col_meas in colnames(df_meas_mri)){
      df_mri_sub<-df_mri_sub[which(df_mri_sub[,col_meas]==df_meas_mri[i,col_meas]),]
    }
    
    # objective variable (some MRI measure)
    var_obj<-as.numeric.factor(df_mri_sub["value"])
    
    formula_lm<-var_obj~1
    for (j in seq(ncol(df_covar_sub))){
      formula_lm<-update(formula_lm,~.+df_covar_sub[,j])
    }
    glmfit<-lm(formula_lm)
    if(ncol(df_covar_sub)>1){
      suppressWarnings(vifactor<-vif(glmfit))
    }else{
      vifactor<-NaN
    }
    
  }
  
  
  
  
  
  n_measures<-nrow(input_measures)
  output<-data.frame(matrix(nrow=0,ncol=(9+ncol(input_measures))))
#  output<-data.frame(matrix(ncol=2+5*n_expvar,nrow=n_measures))
  collabel<-colnames(input_covar)[id_covar+1]
  input_covar<-data.frame(input_covar[,id_covar+1])
  colnames(input_covar)<-collabel
  model_name<-paste(colnames(input_covar),collapse="_")
  for (i in 1:n_measures){
    id_MRI_measure<-1:nrow(input_MRI_data)
    for (j in colnames(input_measures)){
      if (is.na(input_measures[i,j])){
        id_MRI_measure<-intersect(id_MRI_measure,which(is.na(input_MRI_data[,j])))
      }else{
        id_MRI_measure<-intersect(id_MRI_measure,which(input_MRI_data[,j]==input_measures[i,j]))
      }
    }
    MRI_measure<-as.numeric.factor(input_MRI_data[id_MRI_measure,"value"])
    if (length(id_covar)==1){
      glmfit<-lm(MRI_measure~input_covar[,1])
    }else if (length(id_covar)==2){
      glmfit<-lm(MRI_measure~input_covar[,1]+input_covar[,2])
    }else if (length(id_covar)==3){
      glmfit<-lm(MRI_measure~input_covar[,1]+input_covar[,2]+input_covar[,3])
    }else if (length(id_covar)==4){
      glmfit<-lm(MRI_measure~input_covar[,1]+input_covar[,2]+input_covar[,3]+input_covar[,4])
    }
    if (length(id_covar)>=2){
      suppressWarnings(vifactor<-vif(glmfit))
    }else{
      vifactor<-NaN
    }
#    stats<-c(AIC(glmfit),BIC(glmfit))
    xic<-c(AIC(glmfit),BIC(glmfit))
    for (j in 1:n_expvar){
      contrast<-matrix(0L,nrow=1, ncol=length(id_covar)+1)
      contrast[1,j+1]<-1
      ttest<-tryCatch(summary(glht(glmfit, linfct = contrast))$test,
                      error=function(e){return(NaN)},
                      warning=function(e){return(NaN)},
                      silent=T)
      if (length(ttest)==1){
#        stats <-c(stats, rep(NaN,4),vifactor[j])
        output_add<-cbind(input_measures[i,],model_name,colnames(input_covar)[j],
                          NaN,NaN,NaN,NaN,
                          vifactor[j],xic[1],xic[2])
      }else{
#        stats <-c(stats, ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],vifactor[j])
        output_add<-cbind(input_measures[i,],model_name,colnames(input_covar)[j],
                          ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],
                          vifactor[j],xic[1],xic[2])
      }
      colnames(output_add)<-c(colnames(input_measures),"model","exp_var","beta","sigma","t","p","VIF","AIC","BIC")
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
  
  # set up covariate dataframe according to subsetting
  df_covar<-df_clinical[,c("IC_pnTTC",list_covar)]
  for (covar in list_covar){
    df_covar[,covar]<-df_covar[,covar]-mean(df_covar[,covar])
  }
  
  # subset MRI datarame
  df_mri<-df_mri[df_mri$ID_pnTTC %in% list_id_subj]
  
  # df of measurements per subject
  df_meas_mri<-df_mri[df_mri$ID_pnTTC==df_mri[1,ID_pnTTC],]
  df_meas_mri[,c("ID_pnTTC","value")]<-list(NULL)
  
  n_covar<-length(list_covar)
  
  for (i in n_covar:1){
    combn_covar<-combn(n_covar,i)
    n_combn_covar<-dim(combn_covar)[2]
    for (j in 1:combn_covar){
      id_covar<-combn_covar[,j]
      df_covar_sub<-df_covar[,c("ID_pnTTC",list_covar[id_covar])]
      glm_routine(df_mri,df_meas_mri,df_covar_sub)
    }
  }

  
  
#  MRI_data,input_covariate_label=covariate_label,global_covariate=F,dirname,filename){
#  n_covariates<-length(input_covariate_label)


#  output<-MRI_data[which(MRI_data$ID_pnTTC==MRI_data[1,"ID_pnTTC"]),
#                   -c(which(colnames(MRI_data)=="ID_pnTTC"),which(colnames(MRI_data)=="value"))]
#  output<-data.frame(output)
#  colnames(output)<-colnames(MRI_data)[-c(which(colnames(MRI_data)=="ID_pnTTC"),
#                                                  which(colnames(MRI_data)=="value"))]
  clinical_data_subset<-clinical_data
  for (i in 1:n_covariates){
    clinical_data_subset<-clinical_data_subset[which(!is.na(clinical_data_subset[,input_covariate_label[i]])),]
  }
  subject_id_subset<-clinical_data_subset$ID_pnTTC
  
  covariates_data_subset<-data.frame(ID_pnTTC=clinical_data_subset$ID_pnTTC)
  for (i in 1:n_covariates){
    covariates_data_subset<-cbind(covariates_data_subset,clinical_data_subset[,input_covariate_label[i]])
  }
  for (i in 2:ncol(covariates_data_subset)){
    ave<-mean(covariates_data_subset[,i])
    covariates_data_subset[,i]<-covariates_data_subset[,i]-ave
  }
  colnames(covariates_data_subset)[-1]<-input_covariate_label
  
  MRI_data_subset<-data.frame(matrix(ncol=ncol(MRI_data),nrow=0))
  for (i in subject_id_subset){
    MRI_data_subset<-rbind(MRI_data_subset,MRI_data[which(MRI_data$ID_pnTTC==i),])
  }
  colnames(MRI_data_subset)<-colnames(MRI_data)
  
  measures<-MRI_data_subset[which(MRI_data_subset$ID_pnTTC==MRI_data_subset[1,"ID_pnTTC"]),
                            -c(which(colnames(MRI_data_subset)=="ID_pnTTC"),which(colnames(MRI_data_subset)=="value"))]
  measures<-data.frame(measures)
  colnames(measures)<-colnames(MRI_data_subset)[-c(which(colnames(MRI_data_subset)=="ID_pnTTC"),
                                                   which(colnames(MRI_data_subset)=="value"))]
  n_measures<-nrow(measures)
  if (global_covariate){
    global_covariate_data<-read.csv(file.path(input_dir,global_covariate_file))
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
  output<-NULL
  for (i in n_covariates:1){
    n_expvar<-i
    for (j in 1:dim(combn(n_covariates,i))[2]){
      id_covar<-combn(n_covariates,i)[,j]
      if (global_covariate){
        id_covar<-c(id_covar,n_covariates+1)
      }
      output<-rbind(output,GLMroutine(MRI_data_subset,measures,all_covariates_data,id_covar,n_expvar))
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
  
  write.csv(output, file.path(dirname,filename),row.names=F)
  return(output)
}
