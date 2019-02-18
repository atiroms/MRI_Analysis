



library(multcomp)


sem <- function(x) sd(x)/sqrt(length(x))

dir<-"C:/Users/atiro/Dropbox/MRI/Statistics/Structural_FS/Structural_data/20180828_143654_LI_Male_Volume_Subcortex"
#dir<-"C:/Users/atiro/Dropbox/MRI/Statistics/Structural_FS/Structural_data/20180828_143723_LI_Female_Volume_Subcortex"
input_file<-"LI_Structure.csv"



LI_data<-read.csv(file.path(dir,input_file))

measures<-LI_data[which(LI_data$ID_pnTTC==LI_data[1,"ID_pnTTC"]),
                          -c(which(colnames(LI_data)=="ID_pnTTC"),
                             which(colnames(LI_data)=="L_value"),
                             which(colnames(LI_data)=="R_value"),
                             which(colnames(LI_data)=="Laterality_Index"))]

measures$ROI<-as.character(measures$ROI)
measures$L_ROI_ID<-as.character(measures$L_ROI_ID)
measures$R_ROI_ID<-as.character(measures$R_ROI_ID)
subject_id<-unique(LI_data$ID_pnTTC)
subject_id<-subject_id[order(subject_id)]

parent_dir <- "C:/Users/atiro"
common_dir<-file.path(parent_dir,"DropBox/MRI/Statistics/CommonData")
clinical_file <- "CSUB.csv"

clinical_data <- read.csv(file.path(common_dir,clinical_file))

TS<-NULL
Age<-NULL
for (i in subject_id){
  TS<-c(TS,clinical_data[which(clinical_data$ID_pnTTC==i),"W1_Tanner_Stage"])
  Age<-c(Age,clinical_data[which(clinical_data$ID_pnTTC==i),"W1_Age_at_MRI"])
}

id_omit<-which(is.na(TS))

TS<-TS[-id_omit]
Age<-Age[-id_omit]

TS_ave<-mean(TS)
Age_ave<-mean(Age)
TS<-TS-TS_ave
Age<-Age-Age_ave


output<-data.frame(matrix(nrow=0,ncol=8))
colnames(output)<-c("ROI","L_ROI_ID","R_ROI_ID","exp_var","beta","sigma","t","p")

means<-data.frame(matrix(nrow=0,ncol=3))
colnames(means)<-c("ROI","mean","sem")

for (i in 1:nrow(measures)){
  LI<-LI_data[which(LI_data$ROI==measures[i,"ROI"]),"Laterality_Index"]
  means<-rbind(means,cbind(ROI=measures[i,"ROI"],mean=mean(LI),sem=sem(LI)))
  LI<-LI[-id_omit]
  glmfit<-lm(LI~TS*Age)
#  summary(glmfit)
  
  contrast <- matrix(c(0, 1, 0, 0), 1)
  ttest<-summary(glht(glmfit, linfct = contrast))$test
  stats<-cbind(ROI=measures[i,"ROI"],R_ROI_ID=measures[i,"R_ROI_ID"],
               L_ROI_ID=measures[i,"L_ROI_ID"],exp_var="TS",beta=ttest$coefficients[1],
               sigma=ttest$sigma[1],t=ttest$tstat[1],p=ttest$pvalues[1])
  output<-rbind(output,stats)
  
  contrast <- matrix(c(0, 0, 1, 0), 1)
  ttest<-summary(glht(glmfit, linfct = contrast))$test
  stats<-cbind(ROI=measures[i,"ROI"],R_ROI_ID=measures[i,"R_ROI_ID"],
               L_ROI_ID=measures[i,"L_ROI_ID"],exp_var="Age",beta=ttest$coefficients[1],
               sigma=ttest$sigma[1],t=ttest$tstat[1],p=ttest$pvalues[1])
  output<-rbind(output,stats)
  
  contrast <- matrix(c(0, 0, 0, 1), 1)
  ttest<-summary(glht(glmfit, linfct = contrast))$test
  stats<-cbind(ROI=measures[i,"ROI"],R_ROI_ID=measures[i,"R_ROI_ID"],
               L_ROI_ID=measures[i,"L_ROI_ID"],exp_var="TS_Age_interaction",beta=ttest$coefficients[1],
               sigma=ttest$sigma[1],t=ttest$tstat[1],p=ttest$pvalues[1])
  output<-rbind(output,stats)
}


write.csv(output,file.path(dir,"GLM_LI_interaction.csv"))
write.csv(means,file.path(dir,"LI_means.csv"))
