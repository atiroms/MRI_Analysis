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
  name_model<-paste(colnames(df_covar_sub)[-which(colnames(df_covar_sub)=="ID_pnTTC")],collapse="_")
  df_output<-data.frame()
  for (i in seq(nrow(df_meas_mri))){
    df_mri_sub<-df_mri
    for (col_meas in colnames(df_meas_mri)){
      df_mri_sub<-df_mri_sub[which(df_mri_sub[,col_meas]==df_meas_mri[i,col_meas]),]
    }
    
    # Objective variable (some MRI measure)
    var_obj<-as.numeric.factor(df_mri_sub$value)
    
    # Define GLM formula as string first
    for (j in seq(ncol(df_covar_sub)-1)){
      assign(paste("covar",as.character(j),sep="_"),df_covar_sub[,j+1])
      if (j==1){
        str_formula_lm<-"var_obj~covar_1"
      }else{
        str_formula_lm<-paste(str_formula_lm,"+covar_",as.character(j),sep="")
      }
    }
    # Convert model string to formula ovject
    formula_lm<-as.formula(str_formula_lm)
    glmfit<-lm(formula_lm)
    
    # Variance inflation factor
    vifactor<-NaN
    if(ncol(df_covar_sub)>1){
      suppressWarnings(vifactor<-vif(glmfit))
    }
    
    # Calculate stats for each explanatory variables
    for (j in seq(ncol(df_covar_sub)-1)){
      contrast<-matrix(0L,nrow=1, ncol=ncol(df_covar_sub))
      contrast[1,j+1]<-1
      ttest<-tryCatch(summary(glht(glmfit, linfct = contrast))$test,
                      error=function(e){return(NaN)},
                      warning=function(e){return(NaN)},
                      silent=T)
      if (length(ttest)==1){
        df_output_per_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j],
                                  NaN,NaN,NaN,NaN,
                                  vifactor[j],xic[1],xic[2])
      }else{
        df_output_per_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j],
                                  ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],
                                  vifactor[j],xic[1],xic[2])
      }
      colnames(df_output_per_meas)<-c(colnames(df_meas_mri),"model","var_exp","beta","sigma","t","p","VIF","AIC","BIC")
      df_output<-rbind(df_output,df_output_per_meas)
    }
  }
  return(df_output)
}


func_glm<-function(df_mri,data_clinical,list_covar){
  # subset clinical data according to clinical data availability
  df_clinical<-data_clinical$df_clinical
  for (covar in list_covar){
    df_clinical<-df_clinical[which(!is.na(df_clinical[,covar])),]
  }
  list_id_subj<-df_clinical$ID_pnTTC
  
  # Set up covariate dataframe according to subsetting
  df_covar<-df_clinical[,c("ID_pnTTC",list_covar)]
  for (covar in list_covar){
    df_covar[,covar]<-df_covar[,covar]-mean(df_covar[,covar])
  }
  
  # Subset MRI datarame
  df_mri<-df_mri[df_mri$ID_pnTTC %in% list_id_subj,]
  
  # df of measurements per subject
  df_meas_mri<-df_mri[df_mri$ID_pnTTC==df_mri[1,"ID_pnTTC"],]
  df_meas_mri[,c("ID_pnTTC","value")]<-list(NULL)
  
  df_glm<-NULL
  for (i in length(list_covar):1){   # Iterate over number of covariates
    combn_covar<-combn(length(list_covar),i)
    n_combn_covar<-dim(combn_covar)[2]
    for (j in 1:n_combn_covar){   # Iterate over combination of covariates
      id_covar<-combn_covar[,j]
      df_covar_sub<-df_covar[,c("ID_pnTTC",list_covar[id_covar])]
      # Calculate GLM model for each model (set of variables)
      df_glm<-rbind(df_glm,glm_routine(df_mri,df_meas_mri,df_covar_sub)) 
    }
  }
  
  # For each MRI measure, calculate which model exhibits smallest AIC or BIC
  df_glm$AIC_min<-F
  df_glm$BIC_min<-F
  for (i in 1:nrow(df_meas_mri)){
    id_obs<-1:nrow(df_glm)
    for (j in colnames(df_meas_mri)){
      id_obs<-intersect(id_obs,which(df_glm[,j]==df_meas_mri[i,j]))
    }
    id_minAIC<-which.min(df_glm[id_obs,"AIC"])
    df_glm[id_obs[id_minAIC],"AIC_min"]<-T
    id_minBIC<-which.min(df_glm[id_obs,"BIC"])
    df_glm[id_obs[id_minBIC],"BIC_min"]<-T
  }
  
  return(df_glm)
}


#**************************************************
# GLM of FC results into nodes and edges ==========
#**************************************************

glm_fc2graph<-function(df_input,list_node){
  node<-data.frame(label=as.character(list_node),stringsAsFactors = F)
  node<-rowid_to_column(node, "id")
  df_input$from<-as.character(df_input$from)
  df_input$to<-as.character(df_input$to)
  edge<-left_join(df_input, node, by = c("from" = "label")) 
  edge<-edge[,-which(colnames(edge)=="from")]
  edge<-rename(edge, from = id)
  edge<-left_join(edge, node, by = c("to" = "label"))
  edge<-edge[,-which(colnames(edge)=="to")]
  edge<-rename(edge, to = id)
  edge<-rename(edge, weight=beta)
  collabel<-colnames(edge)
  collabel<-collabel[-c(which(collabel=="from"),which(collabel=="to"))]
  collabel<-c("from","to",collabel)
  edge<-edge[,collabel]
  output<-list("node"=node,"edge"=edge)
  return(output)
}

