#parent_dir <- "D:/atiroms"
parent_dir <- "C:/Users/atiro"
#common_dir<-file.path(parent_dir,"DropBox/MRI/Statistics/CommonData")
common_dir<-file.path(parent_dir,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
clinical_file <- "CSUB.csv"


#### Libraries ####
library(mgcv)

#### Clinical data loading ####

clinical_data <- read.csv(file.path(common_dir,clinical_file))

# select subjects with QC'ed T1 image and TS data in at least one timepoint
clinical_data <- clinical_data[which((clinical_data$W1_T1QC==1 & !is.na(clinical_data$W1_Tanner_Max))| (clinical_data$W2_T1QC==1) & !is.na(clinical_data$W2_Tanner_Max)),]

clinical_data_W1 <- clinical_data[which(clinical_data$W1_T1QC==1 & !is.na(clinical_data$W1_Tanner_Max)),c("ID_pnTTC","W1_Age_at_MRI","W1_Tanner_Max")]

clinical_data_W2 <- clinical_data[which(clinical_data$W2_T1QC==1 & !is.na(clinical_data$W2_Tanner_Max)),c("ID_pnTTC","W2_Age_at_MRI","W2_Tanner_Max")]

longitudinal_data <- rbind(data.frame(ID_pnTTC=clinical_data_W1$ID_pnTTC,Wave=1,Age_at_MRI=clinical_data_W1$W1_Age_at_MRI,Tanner_Max=clinical_data_W1$W1_Tanner_Max),data.frame(ID_pnTTC=clinical_data_W2$ID_pnTTC,Wave=2,Age_at_MRI=clinical_data_W2$W2_Age_at_MRI,Tanner_Max=clinical_data_W2$W2_Tanner_Max))

fit_gamm <- gamm(formula = Tanner_Max ~ s(Age_at_MRI), data=longitudinal_data,random=list(ID_pnTTC=~1))

plot(fit_gamm$gam,pages=1)
plot(fit_gamm$lme)

summary(fit_gamm$lme)
summary(fit_gamm$gam)
anova(fit_gamm$gam) 

random_effects <- ranef(fit_gamm$lme)$ID_pnTTC[[1]]


