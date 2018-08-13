#### Description ####

# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs are those processed in DoGTA() of ConnectionAnalysis.R
# 

#### Parameters ####

working_dir <- "G:/MRI/Statistics/GraphMetric"

commondata_dir<-"G:/MRI/Statistics/CommonData"

data_subdir<-"GraphMetric_data"

#metric_file <- "GTA_HO_FC_TSexist_male.csv"
#metric_file <- "GTA_HO_FC_TSexist_female.csv"
#metric_file <- "GTA_Power_FC_TSexist_male.csv"
metric_file <- "GTA_Power_FC_TSexist_female.csv"

roiid_file<-"ROI_All.csv"
roiid_type<-"label_conn"

clinical_file <- "CSUB_Clinical_Data.csv"


#### Libraries ####

library(multcomp)


#### Data Loading ####

metric_data<-read.csv(file.path(working_dir, metric_file))
metrics<-colnames(metric_data)[-1]
n_metrics<-length(metrics)
subject_id<-unique(metric_data$ID_pnTTC)
n_subject<-length(subject_id)

clinical_data <-read.csv(file.path(commondata_dir,clinical_file))
clinical_row<-NULL
for (i in 1:n_subject){
  clinical_row<-c(clinical_row,which(clinical_data$ID_pnTTC==subject_id[i]))
}
clinical_data<-clinical_data[clinical_row,]
rownames(clinical_data)<-NULL


##### Directory Organization ####

if (!file.exists(file.path(working_dir,data_subdir))){
  dir.create(file.path(working_dir, data_subdir))
}

ExpDir<-function(exptype){
  timestamp <- strftime(Sys.time(),"%Y%m%d_%H%M%S")
  data_dir<-file.path(working_dir, data_subdir, paste(timestamp,exptype,sep="_"))
  dir.create(data_dir)
  return(data_dir)
}


#### ID Converter ####

roiid_data<-read.csv(file.path(commondata_dir,roiid_file))
roiid_data$ID_long<-as.character(roiid_data$ID_long)
roiid_data$label<-as.character(roiid_data$label)
roiid_data$label_proper<-as.character(roiid_data$label_proper)
roiid_data$label_long<-as.character(roiid_data$label_long)
roiid_data$label_conn<-as.character(roiid_data$label_conn)
roiid_data$label_fs<-as.character(roiid_data$label_fs)

ConvertID<-function(input,iddata,inputcolumn,outputcolumn){
  idcolname<-data.frame(label=colnames(iddata))
  input_col<-iddata[,which(idcolname$label==inputcolumn)]
  output_col<-iddata[,which(idcolname$label==outputcolumn)]
  output<-data.frame(matrix(nrow=length(input),ncol=1))
  for (i in 1:length(input)){
    input_value<-input[[i]]
    output_row<-which(input_col==input_value)
    output[i,1]<-output_col[output_row]
  }
  output<-output[,1]
  return(output)
}


#### GLM Analysis  ####

DoGLM_GraphMetric<-function(){
  dirname<-ExpDir("GLMClinical")
  output<-data.frame(matrix(ncol=1, nrow=n_metrics))
  output[,1]<-colnames(metric_data)[-1]
  colnames(output)<-c("GTA_metric")
  
  TS<-clinical_data$Tanner_Stage
  age<-clinical_data$Age_at_MRI
  centered_TS<-TS-mean(TS)
  centered_age<-age-mean(age)
  
  output_add<-data.frame(matrix(ncol=6,nrow=n_metrics))
  
  for (i in 1:n_metrics){
    ObjV<-metric_data[,i+1]
    glmfit <- lm(ObjV ~ centered_TS + centered_age)
    contrast <- matrix(c(0, 1, 0), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("TS_Age_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  for (i in 1:n_metrics){
    ObjV<-metric_data[,i+1]
    glmfit <- lm(ObjV ~ centered_TS)
    contrast <- matrix(c(0, 1), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("TS_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  for (i in 1:n_metrics){
    ObjV<-metric_data[,i+1]
    glmfit <- lm(ObjV ~ centered_age)
    contrast <- matrix(c(0, 1), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("Age_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  
  best_model<-data.frame(matrix(ncol=1, nrow=n_metrics))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (k in 1:n_metrics){
      best_id<-which.min(xic[k,])
      if (best_id==1){
        best_model[k,1]<-"TS_Age"
      }else if (best_id==2){
        best_model[k,1]<-"TS"
      }else if (best_id==3){
        best_model[k,1]<-"Age"
      }
    }
    colnames(best_model)<-paste(i,"best_model",sep="_")
    output<-cbind(output,best_model)
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}
