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
#glm_routine<-function(df_mri,df_meas_mri,df_covar_sub){
#  name_model<-paste(colnames(df_covar_sub)[-which(colnames(df_covar_sub)=="ID_pnTTC")],collapse="_")
#  print(paste("  Calculating model",name_model,sep=" "))
#  df_output<-data.frame()
#  for (i in seq(nrow(df_meas_mri))){
#    df_mri_sub<-df_mri
#    for (col_meas in colnames(df_meas_mri)){
#      df_mri_sub<-df_mri_sub[which(df_mri_sub[,col_meas]==df_meas_mri[i,col_meas]),]
#    }
#    
#    # Objective variable (some MRI measure)
#    var_obj<-as.numeric.factor(df_mri_sub$value)
#    
#    # Define GLM formula as string first
#    for (j in seq(ncol(df_covar_sub)-1)){
#      assign(paste("covar",as.character(j),sep="_"),df_covar_sub[,j+1])
#      if (j==1){
#        str_formula_lm<-"var_obj~covar_1"
#      }else{
#        str_formula_lm<-paste(str_formula_lm,"+covar_",as.character(j),sep="")
#      }
#    }
#    # Convert model string to formula ovject
#    formula_lm<-as.formula(str_formula_lm)
#    glmfit<-lm(formula_lm)
#    
#    # Variance inflation factor
#    vifactor<-NaN
#    if(ncol(df_covar_sub)>2){
#      suppressWarnings(vifactor<-vif(glmfit))
#    }
#    
#    # Calculate information criteria
#    xic<-c(AIC(glmfit),BIC(glmfit))
#    
#    # Calculate stats for each explanatory variables
#    for (j in seq(ncol(df_covar_sub)-1)){
#      contrast<-matrix(0L,nrow=1, ncol=ncol(df_covar_sub))
#      contrast[1,j+1]<-1
#      ttest<-tryCatch(summary(glht(glmfit, linfct = contrast))$test,
#                      error=function(e){return(NaN)},
#                      warning=function(e){return(NaN)},
#                      silent=T)
#      if (length(ttest)==1){
#        df_output_per_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j+1],
#                                  NaN,NaN,NaN,NaN,
#                                  vifactor[j],xic[1],xic[2])
#      }else{
#        df_output_per_meas<-cbind(df_meas_mri[i,],name_model,colnames(df_covar_sub)[j+1],
#                                  ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],
#                                  vifactor[j],xic[1],xic[2])
#      }
#      colnames(df_output_per_meas)<-c(colnames(df_meas_mri),"model","var_exp","beta","sigma",
#                                      "t","p","VIF","AIC","BIC")
#      df_output<-rbind(df_output,df_output_per_meas)
#    }
#    #print(paste("  Finished calculating MRI measure ",
#    #            as.character(i)," / ",as.character(nrow(df_meas_mri)),sep=""))
#  }
#  return(df_output)
#}


func_glm<-function(df_mri,data_clinical,list_covar){
  # subset clinical data according to clinical data availability
  df_clinical<-data_clinical$df_clinical
  for (covar in list_covar){
    n_subj_pre<-dim(df_clinical)[1]
    df_clinical<-df_clinical[which(!is.na(df_clinical[,covar])),]
    n_subj_post<-dim(df_clinical)[1]
    print(paste("      Covariate: ",covar," exists in ",as.character(n_subj_post)," / ",
                as.characer(n_subj_pre)," subjects.",sep=""))
  }
  list_id_subj<-df_clinical$ID_pnTTC
  
  # Set up covariate dataframe according to subsetting
  df_covar<-df_clinical[,c("ID_pnTTC",list_covar)]
  df_covar_demean<-df_covar
  for (covar in list_covar){
    df_covar_demean[,covar]<-df_covar_demean[,covar]-mean(df_covar_demean[,covar])
  }
  
  # Subset MRI dataframe
  df_mri<-df_mri[df_mri$ID_pnTTC %in% list_id_subj,]
  
  # df of measurements per subject
  df_meas_mri<-df_mri[df_mri$ID_pnTTC==df_mri[1,"ID_pnTTC"],]
  df_meas_mri[,c("ID_pnTTC","value")]<-list(NULL)
  print(paste("      ",as.characer(dim(df_meas_mri)[1])," measures per subject.",sep=""))
  
  list_model<-list()
  for (i in length(list_covar):1){
    combn_covar<-combn(length(list_covar),i)
    n_combn_covar<-dim(combn_covar)[2]
    for (j in seq(n_combn_covar)){
      list_model<-c(list_model,list(combn_covar[,j]))
    }
  }
  print(paste("  ",as.character(length(list_model)),"      GLM models will be calculated.",sep=""))
  
  df_output_per_model<-data.frame()               # VIF is stored
  df_output_per_meas<-data.frame()                # AIC_min and BIC_min are stored
  df_output_per_model_meas<-data.frame()          # AIC and BIC are stored
  df_output_per_model_expvar_meas<-data.frame()   # beta, sigma, t and p are stored
  
  for (i in seq(length(list_model))){  # Iterate over models (set of covariates)
    model<-list_model[[i]]
    list_covar_sub<-list_covar[model]
    name_model<-paste(list_covar_sub,collapse="_")
    print(paste("      Model with covarite: ",paste(list_covar_sub,collapse=" ")," will be calculated.",sep=""))
    for (j in seq(length(list_covar_sub))){
      assign(paste("covar",as.character(j),sep="_"),df_covar[,list_covar_sub[j]])
      if (j==1){
        str_formula_lm<-"var_obj~covar_1"
      }else{
        str_formula_lm<-paste(str_formula_lm,"+covar_",as.character(j),sep="")
      }
    }
    
    # Convert model string to formula object
    # This formula is iteratively used for all mri measures
    formula_lm<-as.formula(str_formula_lm)
    
    for (j in seq(length(model))){  # Iterate over explanatory variables
      # contrast matrix
      print(paste("      Explanaory variable: ",list_covar_sub[j],sep=""))
      contrast<-matrix(0L,nrow=1, ncol=length(model)+1)
      contrast[1,j+1]<-1
      for (k in seq(nrow(df_meas_mri))){  # Iterate over MRI measures
        df_mri_sub<-df_mri
        for (col_meas in colnames(df_meas_mri)){
          df_mri_sub<-df_mri_sub[which(df_mri_sub[,col_meas]==df_meas_mri[k,col_meas]),]
        }
        # df_mri_sub is df of the i'th measure for all subjects 
        var_obj<-as.numeric.factor(df_mri_sub$value)
        
        glmfit<-lm(formula_lm)
        
        if (j==1){  # Asessed once for a model and MRI measure 
          # Calculate variance information factor
          # VIF is defined one for each model (set of covariates)
          if (k==1){
            vifactor<-NaN
            if(length(model)>1){
              suppressWarnings(vifactor<-vif(glmfit))
            }
            df_output_per_model_row<-data.frame("model"=name_model, "vif"=vifactor)
            df_output_per_model<-rbind(df_output_per_model,df_output_per_model_row)
          }
          
          # Calculate information criteria
          # AIC and BIC are defined for each model and each MRI measure
          df_output_per_model_meas_row<-data.frame(df_meas_mri[k,],"model"=name_model,"aic"=AIC(glmfit),"bic"=BIC(glmfit))
          df_output_per_model_meas<-rbind(df_output_per_model_meas,df_output_per_model_meas_row)
        }
        
        ttest<-tryCatch(summary(glht(glmfit, linfct = contrast))$test,
                        error=function(e){return(NaN)},
                        warning=function(e){return(NaN)},
                        silent=T)
        if (length(ttest)==1){
          df_output_per_model_expvar_meas_row<-cbind(df_meas_mri[k,],name_model,list_covar_sub[j],
                                                     NaN,NaN,NaN,NaN)
        }else{
          df_output_per_model_expvar_meas_row<-cbind(df_meas_mri[k,],name_model,list_covar_sub[j],
                                                     ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1])
        }
        colnames(df_output_per_model_expvar_meas_row)<-c(colnames(df_meas_mri),"model","var_exp","beta","sigma","t","p")
        df_output_per_model_expvar_meas<-rbind(df_output_per_model_expvar_meas,df_output_per_model_expvar_meas_row)
      }
      
    }
  }
  
  #df_glm<-NULL
  #for (i in length(list_covar):1){   # Iterate over number of covariates
  #  combn_covar<-combn(length(list_covar),i)
  #  n_combn_covar<-dim(combn_covar)[2]
  #  for (j in 1:n_combn_covar){   # Iterate over combination of covariates
  #    id_covar<-combn_covar[,j]
  #    df_covar_sub<-df_covar[,c("ID_pnTTC",list_covar[id_covar])]
  #    # Calculate GLM model for each model (set of variables)
  #    df_glm<-rbind(df_glm,glm_routine(df_mri,df_meas_mri,df_covar_sub))
  #  }
  #}
  
  # For each MRI measure, calculate which model exhibits smallest AIC or BIC
  print("      Starting to compare information criteria.")
  for (i in seq(nrow(df_meas_mri))){
    df_output_per_model_meas_sub<-df_output_per_model_meas
    for (col_meas in colnames(df_meas_mri)){
      df_output_per_model_meas_sub<-df_output_per_model_meas_sub[which(df_output_per_model_meas_sub[,col_meas]==df_meas_mri[i,col_meas]),]
    }
    df_output_per_meas_row<-data.frame(df_meas_mri[i,],
                                       "min_aic"=df_output_per_model_meas_sub[which.min(df_output_per_model_meas_sub$aic),"model"],
                                       "min_bic"=df_output_per_model_meas_sub[which.min(df_output_per_model_meas_sub$bic),"model"])
    df_output_per_meas<-rbind(df_output_per_meas,df_output_per_meas_row)
  }
  print("      Finished comparing information criteria")
  
  output<-list("glm"=df_output_per_model_expvar_meas,"ic"=df_output_per_model_meas,
               "min_ic"=df_output_per_meas,"vif"=df_output_per_model,"list_model"=list_modelf)
  
  return(output)
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

