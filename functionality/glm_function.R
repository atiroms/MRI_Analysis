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

func_glm<-function(df_mri,data_clinical,list_covar,df_global_covar=NA,key_global_covar=NA){
  # subset clinical data according to clinical data availability
  df_clinical<-data_clinical$df_clinical
  for (covar in list_covar){
    n_subj_pre<-dim(df_clinical)[1]
    df_clinical<-df_clinical[which(!is.na(df_clinical[,covar])),]
    n_subj_post<-dim(df_clinical)[1]
    print(paste("      Covariate: ",covar," exists in ",as.character(n_subj_post)," / ",
                as.character(n_subj_pre)," subjects.",sep=""))
  }
  list_id_subj<-df_clinical$ID_pnTTC
  
  # Set up covariate dataframe according to subsetting
  df_covar<-df_clinical[,c("ID_pnTTC",list_covar)]
  df_covar_demean<-df_covar
  for (covar in list_covar){
    df_covar_demean[,covar]<-df_covar_demean[,covar]-mean(df_covar_demean[,covar])
  }
  
  # SUbset global covariate dataframe
  if(!is.na(key_global_covar)){
    df_global_covar<-df_global_covar[df_global_covar$ID_pnTTC %in% list_id_subj,]
    colnames(df_global_covar)[2]<-"value"
    df_global_covar_demean<-df_global_covar
    df_global_covar_demean[,"value"]<-df_global_covar_demean[,"value"]-mean(df_global_covar_demean[,"value"])
  }
  
  # Subset MRI dataframe
  df_mri<-df_mri[df_mri$ID_pnTTC %in% list_id_subj,]
  
  # df of measurements per subject
  df_meas_mri<-df_mri[df_mri$ID_pnTTC==df_mri[1,"ID_pnTTC"],]
  df_meas_mri[,c("ID_pnTTC","value")]<-list(NULL)
  print(paste("      There are ",as.character(dim(df_meas_mri)[1])," measures per subject.",sep=""))
  
  list_model<-list()
  for (i in length(list_covar):1){
    combn_covar<-combn(length(list_covar),i)
    n_combn_covar<-dim(combn_covar)[2]
    for (j in seq(n_combn_covar)){
      list_model<-c(list_model,list(combn_covar[,j]))
    }
  }
  print(paste("      ",as.character(length(list_model))," GLM models will be calculated.",sep=""))
  
  df_output_per_model<-data.frame()               # VIF is stored
  df_output_per_meas<-data.frame()                # AIC_min and BIC_min are stored
  df_output_per_model_meas<-data.frame()          # AIC and BIC are stored
  df_output_per_model_expvar_meas<-data.frame()   # beta, sigma, t and p are stored
  
  for (i in seq(length(list_model))){  # Iterate over models (set of covariates)
    model<-list_model[[i]]
    list_covar_sub<-list_covar[model]
    
    name_model<-paste(list_covar_sub,collapse="_")
    name_model_pint<-paste(list_covar_sub,collapse=" and ")

    for (j in seq(length(list_covar_sub))){
      assign(paste("covar",as.character(j),sep="_"),df_covar_demean[,list_covar_sub[j]])
      if (j==1){
        str_formula_lm<-"var_obj~covar_1"
      }else{
        str_formula_lm<-paste(str_formula_lm,"+covar_",as.character(j),sep="")
      }
    }
    if (!is.na(key_global_covar)){
      j<-j+1
      assign(paste("covar",as.character(j),sep="_"),df_global_covar_demean[,"value"])
      str_formula_lm<-paste(str_formula_lm,"+covar_",as.character(j),sep="")
      name_model<-paste(name_model,key_global_covar,sep="_")
      name_model_print<-paste(name_model_print,key_global_covar,sep=" and ")
    }
    print(paste("      Model with covarite: ",name_model_print," will be calculated.",sep=""))
    
    # Convert model string to formula object
    # This formula is iteratively used for all mri measures
    formula_lm<-as.formula(str_formula_lm)
    
    for (j in seq(length(model))){  # Iterate over explanatory variables
      # contrast matrix
      print(paste("      Evaluating explanaory variable: ",list_covar_sub[j],sep=""))
      if(!is.na(key_global_covar)){
        contrast<-matrix(0L,nrow=1, ncol=length(model)+1)
      }else{
        contrast<-matrix(0L,nrow=1, ncol=length(model)+2)
      }
      
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
            if(!is.na(key_global_covar)){
              if(length(model)>2){
                suppressWarnings(vifactor<-vif(glmfit))
              }
            }else{
              if(length(model)>1){
                suppressWarnings(vifactor<-vif(glmfit))
              }
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
               "min_ic"=df_output_per_meas,"vif"=df_output_per_model,"list_model"=list_model)
  
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

