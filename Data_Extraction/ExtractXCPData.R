#****************************************
# Description ===========================
#****************************************

# R script to extract ROI average BOLD signal from XCP data
# execute extract.xcp() for data extraction.


#****************************************
# Parameters ============================
#****************************************
path.exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
dir.in   <- "02_ts_test_5sub"
dir.out  <- "03_ts_extract"

prefix.file.input<-"sub-"
suffix.file.input<-"_power264_ts.1D"

list.id.subj<-c(14,19,26,28,29)

atlas.roi<-"Power"
list.id.roi<-seq(264)


#****************************************
# Root path finder============
#****************************************
func.path<-function(list.path.root = c("D:/atiroms","C:/Users/atiro"),
                         path.exp.=path.exp,
                         dir.in.=dir.in,
                         dir.out.=dir.out){
  path.root<-NA
  for(p in list.path.root){
    if(file.exists(p)){
      path.root<-p
    }
  }
  if(is.na(path.root)){
    print("Error: root path could not be found.")
  }
  path.script <- file.path(path.root,"GitHub/MRI_Analysis")
  path.common <- file.path(path.root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
  path.in     <- file.path(path.root,path.exp.,dir.in.)
  path.out    <- file.path(path.root,path.exp.,dir.out.)
  output <- list("script"=path.script,"input"=path.in,"output"=path.out,"common"=path.common)
  return(output)
}

paths<-func.path()


#****************************************
# Function library ======================
#****************************************
source(file.path(paths$script,"Functionalities/Functions.R"))


#****************************************
# Data extraction =======================
#****************************************
extract.xcp<-function(paths.=paths,
                      dict.roi.=dict.roi,
                      prefix.file.input.=prefix.file.input,
                      suffix.file.input.=suffix.file.input,
                      list.id.subj.=list.id.subj,
                      atlas.roi.=atlas.roi,
                      list.id.roi.=list.id.roi
                      ){
  
  nullobj<-func.createdirs(paths)
  dict.roi<-func.dict.roi(paths)
  
  output<-data.frame(matrix(ncol=length(list.id.roi.)+2, nrow=0))
  for (id.subj in list.id.subj){
    file.input<-paste(prefix.file.input., sprintf("%05d", id.subj), suffix.file.input., sep="")
    ts.subj<-read.csv(file.path(paths.$input,"output",file.input),header=F,sep=" ")
    df.roi<-dict.roi.[which(dict.roi.$Atlas==atlas.roi.),]
    list.id.roi<-as.character(df.roi$ID_long)
    ts.subj<-cbind(id.subj, seq(dim(ts.subj)[1]),ts.subj)
    colnames(ts.subj)<- c("ID_pnTTC","timeframe",list.id.roi)
    output<-rbind(output, ts.subj)
  }
  
  write.csv(output, file.path(paths.$output,"output","timeseries.csv"),row.names=F)
  return(output)
}