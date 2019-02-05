#### Description ####

# R script for loading clinical data


#### Parameters ####

common_dir<-file.path(parent_dir,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
clinical_file <- "CSUB.csv"


#### Libraries ####


#### Clinical data loading ####

clinical_data <- read.csv(file.path(common_dir,clinical_file))
for (i in 1:ncol(subject_subset)){
  clinical_data<-clinical_data[which(clinical_data[,colnames(subject_subset)[i]]==subject_subset[1,i]),]
}
subject_id<-clinical_data[,"ID_pnTTC"]
id_data <- clinical_data[,1:(which(names(clinical_data)=="Clinical")-1)]
clinical_data <- clinical_data[,(-2):(-which(names(clinical_data)=="Clinical"))]
#heatmap(as.matrix(clinical_data[,-1]),col=rich.colors(100),scale="column",Colv = NA, Rowv = NA)
n_subject<-length(subject_id)
n_clinical_data<-ncol(clinical_data)-1
