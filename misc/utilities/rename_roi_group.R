

path_root<-"D:/atiroms/Dropbox/MRI_img/pnTTC/puberty/common"
file_src<-"ROI.csv"
file_convert<-"convertROI.csv"
file_dst<-"ROI_new.csv"


extract_roigroup<-function(){
  df_src<-read.csv(file.path(path_root,file_src))
  
  list_atlas<-sort(unique(df_src$atlas))
  
  df_convert<-data.frame(matrix(nrow=0,ncol=4))
  colnames(df_convert)<-c("atlas","group_1","group_2","group_3")
  
  for (atlas in list_atlas){
    id_roi_atlas<-which(df_src$atlas==atlas)
    list_group_1<-as.character(sort(unique(df_src[id_roi_atlas,"group_1"])))
    if(any(is.na(df_src[id_roi_atlas,"group_1"]))){
      list_group_1<-c(list_group_1,"NA")
    }
    df_convert<-rbind(df_convert,
                      data.frame(atlas=atlas,group_1=list_group_1,group_2=NA))
  }
  
  write.csv(df_convert,file.path(path_root,file_convert),row.names=F)
}


convert_roigroup<-function(){
  df_src<-read.csv(file.path(path_root,file_src))
  
  list_atlas<-as.character(sort(unique(df_src$atlas)))
  
  df_convert<-read.csv(file.path(path_root,file_convert),stringsAsFactors = F)
  
  for (atlas in list_atlas){
    id_roi_atlas<-which(df_src$atlas==atlas)

    for (id_roi in id_roi_atlas){
      if (is.na(df_src[id_roi,"group_1"])){
        df_src[id_roi,c("group_2","group_3")]<-df_convert[df_convert$atlas==atlas & is.na(df_convert$group_1),c("group_2","group_3")]
      }else{
        id_match<-which(df_convert$atlas==atlas & df_convert$group_1==df_src[id_roi,"group_1"])
        df_src[id_roi,c("group_2","group_3")]<-df_convert[id_match,c("group_2","group_3")]
      }
    }
  }
  write.csv(df_src,file.path(path_root,file_dst),row.names = F)
}