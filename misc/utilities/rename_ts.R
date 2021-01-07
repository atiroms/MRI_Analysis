#**************************************************
# Description =====================================
#**************************************************
# Rename Harvard-Oxford time-series data.

library(data.table)

#list_target_dir<-c("D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/400_ts_acompcor",
#                   "D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/410_ts_acompcor_gsr",
#                   "D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/420_ts_aroma",
#                   "D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/430_ts_aroma_gsr")
list_target_dir<-c("D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/440_ts_36p")

for (target_dir in list_target_dir){
  df_src<-as.data.frame(fread(file.path(target_dir,"output","atl-HarvardOxford_ts.csv")))
  colnames(df_src)<-c(colnames(df_src)[1:3],sprintf("ho112_%05d",seq(112)))
  write.csv(df_src,file.path(target_dir,"output","atl-ho112_ts.csv"),row.names=F)
}
                   