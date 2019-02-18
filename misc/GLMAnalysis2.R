#### Description ####

# R script to general linear model analysis on FreeSurfer data.
# Execute **** to perform each analysis
# 


#### Parameters ####

working_dir <- "D:/MRI/Statistics/GLM2"
#working_dir <- "G:/MRI/Statistics/GLM2"
commondata_dir<-"D:/MRI/Statistics/CommonData"
#commondata_dir<-"G:/MRI/Statistics/CommonData"

data_subdir<-"GLM2_data"

volume_file <- "Volume_208.csv"
thickness_file <- "Thickness_208.csv"
area_file <- "Area_208.csv"

clinical_file <- "Clinical_208.csv"

roiid_file<-"ROI_All.csv"
roiid_type<-"label_fs"

p_uncorrected<-0.05


#### Libraries ####

library(multcomp)


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


#### Data Loading ####

volume_data<-read.csv(file.path(working_dir,volume_file))
thickness_data<-read.csv(file.path(working_dir,thickness_file))
area_data<-read.csv(file.path(working_dir,area_file))
clinical_data <-read.csv(file.path(working_dir,clinical_file))


#colnames(structural_data)[-1]<-ConvertID(colnames(structural_data)[-1],roiid_data,roiid_type,"ID_long")
subject_id<-unique(volume_data$ID_pnTTC)
n_subject<-length(subject_id)


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


#### Data centering ####
sex<-clinical_data$Sex


#### Volume Analysis ####

VolumeAnalysis<-function(){
  dirname<-ExpDir("VolumeAnalysis")
  output<-data.frame(matrix(ncol=2, nrow=(ncol(volume_data)-1)))
  output[,1]<-ConvertID(colnames(volume_data)[-1],roiid_data,roiid_type,"ID_long")
  output[,2]<-ConvertID(colnames(volume_data)[-1],roiid_data,roiid_type,"label_proper")
#  colnames(output)<-c("ROI_ID","ROI_proper","male_beta","male_sigma","male_t","male_p","female_beta","female_sigma","female_t","female_p")
  colnames(output)<-c("ROI_ID","ROI_proper")
  for (j in 1:2){
    if (j==1){
      sex_label<-"male"
    }else{
      sex_label<-"female"
    }
    clinical_data_subset<-clinical_data[which(sex==j),]
    TS<-clinical_data_subset$Tanner_Stage
    age<-clinical_data_subset$Age_at_MRI
    TBV<-clinical_data_subset$SPM_TBV
    centered_TS<-TS-mean(TS)
    centered_age<-age-mean(age)
    centered_TBV<-TBV-mean(TBV)
    
    sexwise_data<-volume_data[which(sex==j),]
    output_add<-data.frame(matrix(ncol=6,nrow=(ncol(volume_data)-1)))
    colnames(output_add)<-paste(sex_label,c("beta","sigma","t","p","AIC","BIC"),sep="_")
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      glmfit <- lm(ObjV ~ centered_TS + centered_age + centered_TBV)
#      summary(glmfit)
      contrast <- matrix(c(0, 1, 0, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      glmfit <- lm(ObjV ~ centered_TS + centered_TBV)
      #      summary(glmfit)
      contrast <- matrix(c(0, 1, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]

      glmfit <- lm(ObjV ~ centered_age + centered_TBV)
      #      summary(glmfit)
      contrast <- matrix(c(0, 1, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(volume_data)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in c('male','female')){
      xmale_xic<-xic[, grep(paste('^',j,sep=""), names(xic))]
      for (k in 1:(ncol(volume_data)-1)){
        best_id<-which.min(xmale_xic[k,])
        if (best_id==1){
          best_model[k,1]<-"TS_Age"
        }else if (best_id==2){
          best_model[k,1]<-"TS"
        }else if (best_id==3){
          best_model[k,1]<-"Age"
        }
      }
      colnames(best_model)<-paste(j,i,"best_model",sep="_")
      output<-cbind(output,best_model)
    }
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}


#### Volume Analysis without TBV ####

VolumeAnalysis2<-function(){
  dirname<-ExpDir("VolumeAnalysis2")
  output<-data.frame(matrix(ncol=2, nrow=(ncol(volume_data)-1)))
  output[,1]<-ConvertID(colnames(volume_data)[-1],roiid_data,roiid_type,"ID_long")
  output[,2]<-ConvertID(colnames(volume_data)[-1],roiid_data,roiid_type,"label_proper")
  #  colnames(output)<-c("ROI_ID","ROI_proper","male_beta","male_sigma","male_t","male_p","female_beta","female_sigma","female_t","female_p")
  colnames(output)<-c("ROI_ID","ROI_proper")
  for (j in 1:2){
    if (j==1){
      sex_label<-"male"
    }else{
      sex_label<-"female"
    }
    clinical_data_subset<-clinical_data[which(sex==j),]
    TS<-clinical_data_subset$Tanner_Stage
    age<-clinical_data_subset$Age_at_MRI
#    TBV<-clinical_data_subset$SPM_TBV
    centered_TS<-TS-mean(TS)
    centered_age<-age-mean(age)
#    centered_TBV<-TBV-mean(TBV)
    
    sexwise_data<-volume_data[which(sex==j),]
    output_add<-data.frame(matrix(ncol=6,nrow=(ncol(volume_data)-1)))
    colnames(output_add)<-paste(sex_label,c("beta","sigma","t","p","AIC","BIC"),sep="_")
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
#      glmfit <- lm(ObjV ~ centered_TS + centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS + centered_age)
      #      summary(glmfit)
#      contrast <- matrix(c(0, 1, 0, 0), 1)
      contrast <- matrix(c(0, 1, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
#      glmfit <- lm(ObjV ~ centered_TS + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS)
      #      summary(glmfit)
#      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
#      glmfit <- lm(ObjV ~ centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_age)
      #      summary(glmfit)
#      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(volume_data)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in c('male','female')){
      xmale_xic<-xic[, grep(paste('^',j,sep=""), names(xic))]
      for (k in 1:(ncol(volume_data)-1)){
        best_id<-which.min(xmale_xic[k,])
        if (best_id==1){
          best_model[k,1]<-"TS_Age"
        }else if (best_id==2){
          best_model[k,1]<-"TS"
        }else if (best_id==3){
          best_model[k,1]<-"Age"
        }
      }
      colnames(best_model)<-paste(j,i,"best_model",sep="_")
      output<-cbind(output,best_model)
    }
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}


#### Thickness Analysis  ####

ThicknessAnalysis<-function(){
  dirname<-ExpDir("ThicknessAnalysis")
  output<-data.frame(matrix(ncol=2, nrow=(ncol(thickness_data)-1)))
  output[,1]<-ConvertID(colnames(thickness_data)[-1],roiid_data,roiid_type,"ID_long")
  output[,2]<-ConvertID(colnames(thickness_data)[-1],roiid_data,roiid_type,"label_proper")
  #  colnames(output)<-c("ROI_ID","ROI_proper","male_beta","male_sigma","male_t","male_p","female_beta","female_sigma","female_t","female_p")
  colnames(output)<-c("ROI_ID","ROI_proper")
  for (j in 1:2){
    if (j==1){
      sex_label<-"male"
    }else{
      sex_label<-"female"
    }
    clinical_data_subset<-clinical_data[which(sex==j),]
    TS<-clinical_data_subset$Tanner_Stage
    age<-clinical_data_subset$Age_at_MRI
    #    TBV<-clinical_data_subset$SPM_TBV
    centered_TS<-TS-mean(TS)
    centered_age<-age-mean(age)
    #    centered_TBV<-TBV-mean(TBV)
    
    sexwise_data<-thickness_data[which(sex==j),]
    output_add<-data.frame(matrix(ncol=6,nrow=(ncol(thickness_data)-1)))
    colnames(output_add)<-paste(sex_label,c("beta","sigma","t","p","AIC","BIC"),sep="_")
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_TS + centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS + centered_age)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0, 0), 1)
      contrast <- matrix(c(0, 1, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_TS + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_age)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(thickness_data)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in c('male','female')){
      xmale_xic<-xic[, grep(paste('^',j,sep=""), names(xic))]
      for (k in 1:(ncol(thickness_data)-1)){
        best_id<-which.min(xmale_xic[k,])
        if (best_id==1){
          best_model[k,1]<-"TS_Age"
        }else if (best_id==2){
          best_model[k,1]<-"TS"
        }else if (best_id==3){
          best_model[k,1]<-"Age"
        }
      }
      colnames(best_model)<-paste(j,i,"best_model",sep="_")
      output<-cbind(output,best_model)
    }
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}


#### Area Analysis  ####

AreaAnalysis<-function(){
  dirname<-ExpDir("AreaAnalysis")
  output<-data.frame(matrix(ncol=2, nrow=(ncol(area_data)-1)))
  output[,1]<-ConvertID(colnames(area_data)[-1],roiid_data,roiid_type,"ID_long")
  output[,2]<-ConvertID(colnames(area_data)[-1],roiid_data,roiid_type,"label_proper")
  #  colnames(output)<-c("ROI_ID","ROI_proper","male_beta","male_sigma","male_t","male_p","female_beta","female_sigma","female_t","female_p")
  colnames(output)<-c("ROI_ID","ROI_proper")
  for (j in 1:2){
    if (j==1){
      sex_label<-"male"
    }else{
      sex_label<-"female"
    }
    clinical_data_subset<-clinical_data[which(sex==j),]
    TS<-clinical_data_subset$Tanner_Stage
    age<-clinical_data_subset$Age_at_MRI
    #    TBV<-clinical_data_subset$SPM_TBV
    centered_TS<-TS-mean(TS)
    centered_age<-age-mean(age)
    #    centered_TBV<-TBV-mean(TBV)
    
    sexwise_data<-area_data[which(sex==j),]
    output_add<-data.frame(matrix(ncol=6,nrow=(ncol(area_data)-1)))
    colnames(output_add)<-paste(sex_label,c("beta","sigma","t","p","AIC","BIC"),sep="_")
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_TS + centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS + centered_age)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0, 0), 1)
      contrast <- matrix(c(0, 1, 0), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_TS + centered_TBV)
      glmfit <- lm(ObjV ~ centered_TS)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"TS_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
    for (i in 1:(ncol(sexwise_data)-1)){
      ObjV<-sexwise_data[,i+1]
      #      glmfit <- lm(ObjV ~ centered_age + centered_TBV)
      glmfit <- lm(ObjV ~ centered_age)
      #      summary(glmfit)
      #      contrast <- matrix(c(0, 1, 0), 1)
      contrast <- matrix(c(0, 1), 1)
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
      output_add[i,]<-stats
      collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
      collabel<-c(collabel,"AIC","BIC")
      collabel<-paste(sex_label,"Age_model",collabel,sep="_")
      colnames(output_add)<-collabel
    }
    output<-cbind(output,output_add)
    
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(area_data)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in c('male','female')){
      xmale_xic<-xic[, grep(paste('^',j,sep=""), names(xic))]
      for (k in 1:(ncol(area_data)-1)){
        best_id<-which.min(xmale_xic[k,])
        if (best_id==1){
          best_model[k,1]<-"TS_Age"
        }else if (best_id==2){
          best_model[k,1]<-"TS"
        }else if (best_id==3){
          best_model[k,1]<-"Age"
        }
      }
      colnames(best_model)<-paste(j,i,"best_model",sep="_")
      output<-cbind(output,best_model)
    }
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}

