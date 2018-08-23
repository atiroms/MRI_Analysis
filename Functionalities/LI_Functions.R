#### Description ####

# R script for performing general linear model fitting calculations


#### Libraries ####


### Laterality Index Calculations ####

CommonLI<-function(MRI_data,roiid_colname,dirname,filename){
  # DO NOT pass MRI_data including ROI information in columns other than roiid_dolname.
  measures<-MRI_data[which(MRI_data$ID_pnTTC==MRI_data[1,"ID_pnTTC"]),
                     -c(which(colnames(MRI_data)=="ID_pnTTC"),which(colnames(MRI_data)=="value"))]
  measures<-data.frame(measures)
  colnames(measures)<-colnames(MRI_data)[-c(which(colnames(MRI_data)=="ID_pnTTC"),
                                            which(colnames(MRI_data)=="value"))]
  roi_id<-as.character(measures[,which(colnames(measures)==roiid_colname)])
  roi_label<-ConvertID(roi_id,roi_data,"ID_long","label_proper")
  roi_id_left<-grep("^L ",roi_label)
  measures_bilateral<-measures[roi_id_left,]
  measures_bilateral<-data.frame(measures_bilateral)
  colnames(measures_bilateral)<-colnames(measures)
  colnames(measures_bilateral)[which(colnames(measures_bilateral)==roiid_colname)]<-"L_ROI_ID"
  measures_bilateral$roi_label_left<-roi_label[roi_id_left]
  measures_bilateral$ROI_label<-substr(measures_bilateral$roi_label_left,3,1000)
  measures_bilateral$roi_label_right<-paste("R",measures_bilateral$ROI_label,sep=" ")
  measures_bilateral$R_ROI_ID<-NA
  for (i in 1:nrow(measures_bilateral)){    
    measures_bilateral[i,"R_ROI_ID"]<-roi_id[which(roi_label==measures_bilateral[i,"roi_label_right"])[[1]]]
  }
  if (ncol(measures)>1){
    colindex_notroi<-2:ncol(measures)
  }else{
    colindex_notroi<-NULL
  }
  measures_bilateral<-cbind(measures_bilateral[,"ROI_label"],
                            measures_bilateral[,colindex_notroi],
                            measures_bilateral[,c("L_ROI_ID","R_ROI_ID")])
  colnames(measures_bilateral)[1:ncol(measures)]<-colnames(measures)
  output<-data.frame(matrix(nrow=0,ncol=ncol(measures_bilateral)+2))
  for (i in subject_id){
    for (j in 1:nrow(measures_bilateral)){
      left_id<-intersect(which(MRI_data[,"ID_pnTTC"]==i),
                         which(MRI_data[,roiid_colname]==measures_bilateral[j,"L_ROI_ID"]))
      right_id<-intersect(which(MRI_data[,"ID_pnTTC"]==i),
                          which(MRI_data[,roiid_colname]==measures_bilateral[j,"R_ROI_ID"]))
      for (k in colindex_notroi){
        left_id<-intersect(left_id,which(MRI_data[,k+1]==measures_bilateral[j,k]))
        right_id<-intersect(right_id,which(MRI_data[,k+1]==measures_bilateral[j,k]))
      }
      left_measure<-as.numeric.factor(MRI_data[left_id,"value"])
      right_measure<-as.numeric.factor(MRI_data[right_id,"value"])
      li<-(left_measure-right_measure)/(left_measure+right_measure)
      output_add<-cbind(ID_pnTTC=i,measures_bilateral[j,],L_value=left_measure,R_value=right_measure,Laterality_Index=li)
      output<-rbind(output,output_add)
    }
  }
  
  write.csv(output, file.path(dirname,filename),row.names=F)
  return(output)
}