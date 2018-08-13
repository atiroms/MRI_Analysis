#### Description ####

# R script to analyze clinical data.


#### Parameters ####

#parent_dir <- "D:/atiroms"
parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/TS_Clinical")
output_dir <- file.path(input_dir,"Clinical_data")

subject_subset <- data.frame(W1_T1QC=1)

#### Libraries ####
library(ggplot2)
library(ggpubr)
library(plyr)

#barplot(as.numeric(clinical_data$Tanner_Stage))


#### Data loading ####

source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
clinical_data$Sex<-factor(clinical_data$Sex, levels=c(1,2))
clinical_data <-clinical_data[which(!is.na(clinical_data$Sex)),]

paste("Number of subjects with clnical data:", nrow(clinical_data), sep=" ")
paste("Number of male subjects with clinical data:", length(which(clinical_data$Sex==1)), sep=" ")
paste("Number of female subjects with clinical data:", length(which(clinical_data$Sex==2)), sep=" ")


#### TS sex-wise barplot ####

TS_data<-clinical_data[which(!is.na(clinical_data$W1_Tanner_Stage)),]
paste("Number of subjects with Tanner stage data:", nrow(TS_data), sep=" ")
paste("Number of male subjects with Tanner stage data:", length(which(TS_data$Sex==1)), sep=" ")
paste("Number of female subjects with Tanner stage data:", length(which(TS_data$Sex==2)), sep=" ")


#TS_data$W1_Tanner_Stage<-as.character(TS_data$W1_Tanner_Stage)

#ggplot(data=TS_data,aes(W1_Tanner_Stage,fill=Sex)) +
#  geom_bar(width=0.5,position=position_dodge(width=0.6)) +
#ggplot(as.data.frame(with(TS_data, table(W1_Tanner_Stage = factor(W1_Tanner_Stage), Sex = factor(Sex)))),aes(factor(W1_Tanner_Stage), y = Freq, fill = factor(Sex))) +
ggplot(as.data.frame(with(TS_data, table(W1_Tanner_Stage = factor(W1_Tanner_Stage), Sex))),aes(factor(W1_Tanner_Stage), y = Freq, fill = Sex)) +
  geom_col(width=0.7, position = position_dodge(width=0.8)) + 
  scale_fill_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral"),drop=FALSE) +
  ggtitle("Tanner Stage for pn-TTC 1st Wave") +
  xlab("Tanner Stage") +
  ylab("Number of Subjects") +
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1), legend.position=c(0.9,0.9),panel.grid.major.x=element_blank())



#### TS-Age Correlation ####

#TS_Age_corr <- cor.test(x=clinical_data$Age_at_MRI, y=clinical_data$Tanner_Stage, method = 'spearman', exact=F)

#ggscatter(clinical_data, x = "Age_at_MRI", y = "Tanner_Stage",
#          add = "reg.line", conf.int = TRUE, size=0.5,
#          size=0.5,
#          cor.coef = TRUE, cor.method = "pearson",
#          title="pn-TTC 1st Wave TS-Age",
#          xlim=c(3700,4900),
#          xlab = "Age (day)", ylab = "Tanner Stage")


#### TS-Age sex-wise Correlation ####

TS_data <-TS_data[which(!is.na(TS_data$W1_Age_at_MRI)),]
#TS_data<-TS_data[which(TS_data$T1exist==1),]
TS_data_male <-TS_data[which(TS_data$Sex==1),]
TS_data_female <-TS_data[which(TS_data$Sex==2),]
TS_Age_Pearson_male <- cor.test(x=TS_data_male$W1_Age_at_MRI, y=TS_data_male$W1_Tanner_Stage, method = 'pearson', exact=F)
TS_Age_Pearson_female <- cor.test(x=TS_data_female$W1_Age_at_MRI, y=TS_data_female$W1_Tanner_Stage, method = 'pearson', exact=F)
TS_Age_Spearman_male <- cor.test(x=TS_data_male$W1_Age_at_MRI, y=TS_data_male$W1_Tanner_Stage, method = 'spearman', exact=F)
TS_Age_Spearman_female <- cor.test(x=TS_data_female$W1_Age_at_MRI, y=TS_data_female$W1_Tanner_Stage, method = 'spearman', exact=F)


#### TS-Age Regression ####

TS_Age_LM_male<-list(lm(W1_Tanner_Stage ~ W1_Age_at_MRI, data=TS_data_male),lm(W1_Tanner_Stage ~ poly(W1_Age_at_MRI,2), data=TS_data_male),lm(W1_Tanner_Stage ~ poly(W1_Age_at_MRI,3), data=TS_data_male))
TS_Age_LM_female<-list(lm(W1_Tanner_Stage ~ W1_Age_at_MRI, data=TS_data_female),lm(W1_Tanner_Stage ~ poly(W1_Age_at_MRI,2), data=TS_data_female),lm(W1_Tanner_Stage ~ poly(W1_Age_at_MRI,3), data=TS_data_female))
male_AIC<-NULL
female_AIC<-NULL
for (i in 1:3){
  male_AIC<-c(male_AIC,paste("male TS-Age", i ,"th polynomial regression AIC: ", AIC(TS_Age_LM_male[[i]]), sep=" "))
  female_AIC<-c(female_AIC,paste("female TS-Age", i ,"th polynomial regression AIC: ", AIC(TS_Age_LM_female[[i]]), sep=" "))
}
male_AIC
female_AIC

#ggplot(TS_data,aes(x=W1_Age_at_MRI, y=W1_Tanner_Stage)) +
#  geom_point(aes(color = Sex, shape = Sex), show.legend=F) +
ggplot(TS_data,aes(x=W1_Age_at_MRI, y=W1_Tanner_Stage, colour=Sex,fill=Sex)) +
  geom_point() +
  scale_colour_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral")) +
  scale_fill_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral")) +
  geom_smooth(method = "lm") +
#  geom_rug(aes(color =Sex), show.legend=F) +
  ggtitle("Tanner Stage vs Age for pn-TTC 1st Wave") +
  xlab("Age (day)") +
  ylab("Tanner Stage") +
#  guides(fill = guide_legend(reverse=TRUE)) +
#  scale_fill_discrete(breaks=c("Male","Female")) +
# guides(fill=guide_legend(title=NULL)) +
#  scale_fill_discrete(name="Sex") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),legend.justification=c(0,1), legend.position=c(0.1,0.9),legend.direction="horizontal",panel.grid.minor=element_blank())


#### TS residual-Age Plot####

TS_data_male$W1_Tanner_Stage_residual<-TS_Age_LM_male[[1]]$residuals
TS_data_female$W1_Tanner_Stage_residual<-TS_Age_LM_female[[1]]$residuals
TS_data<-rbind(TS_data_male,TS_data_female)
TS_data<-TS_data[order(TS_data$ID_pnTTC),]

ggplot(TS_data,aes(x=W1_Age_at_MRI, y=W1_Tanner_Stage_residual, colour=Sex,fill=Sex)) +
  geom_point() +
  scale_colour_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral")) +
  scale_fill_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral")) +
  ggtitle("Tanner Stage Residual vs Age for pn-TTC 1st Wave") +
  xlab("Age (day)") +
  ylab("Tanner Stage Residual") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),legend.justification=c(0,0), legend.position=c(0.01,0),panel.grid.minor=element_blank())
