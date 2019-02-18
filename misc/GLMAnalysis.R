####parameters###
#global covariate to be regressed out
covar3_name="BrainSegVolNotVent"
#covar3_name="AEIQ"
#covariate to be assessed
covar4_name="AA83SDQ_pb"
#working_dir <- "C:/Users/atiro/Documents/MRI/Statistics"
working_dir <- "G:/MRI/Statistics"
input_file_name <- "Freesurfer_segmentation02.csv"
output_file_name <- "GLM_Data.csv"
non_ROI_column <-7
#################

library(multcomp)
setwd(working_dir)
data <-read.csv(input_file_name)
femaleratio <- sum(data$Sex == 2)/nrow(data)
#summary(data)
time <- strftime(Sys.time(),"%Y_%m_%d_%H_%M_%S")
dir.create(paste(working_dir, "/GLM_data", sep=""))
data_dir=paste(working_dir, "/GLM_data/", time, sep="")
dir.create(data_dir)
dir.create(paste(data_dir, "/2covar", sep=""))
dir.create(paste(data_dir, "/2covar/graph_age", sep=""))
dir.create(paste(data_dir, "/3covar", sep=""))
dir.create(paste(data_dir, "/3covar/graph_age", sep=""))
dir.create(paste(data_dir, "/3covar/graph_covar3", sep=""))
dir.create(paste(data_dir, "/4covar", sep=""))
dir.create(paste(data_dir, "/4covar/graph_age", sep=""))
dir.create(paste(data_dir, "/4covar/graph_covar4", sep=""))


sex<-data$Sex
age<-(data$Age_at_MRI)
centered_age <- age-mean(age)
covar3 <- data[[covar3_name]]
centered_covar3 <- covar3-mean(covar3)
covar4 <- data[[covar4_name]]
centered_covar4 <- covar4-mean(covar4)


############################
###sex, age as regressors###
############################

statsdata.2covar<-data.frame(ROI=NA,Intercept_Coefficient=NA, Intercept_Sigma=NA, Intercept_tstat=NA, Intercept_p=NA, Sex_Diff_of_Age_Slope_Coefficient=NA, Sex_Diff_of_Age_Slope_Sigma=NA, Sex_Diff_of_Age_Slope_tstat=NA, Sex_Diff_of_Age_Slope_p=NA)

for (i in (non_ROI_column+1):ncol(data)){
	ROI_Name <- colnames(data)[i]
	meas<-data[,i]
	maledata<-subset(data, Sex==1, c("Age_at_MRI", ROI_Name))
	femaledata<-subset(data, Sex==2, c("Age_at_MRI", ROI_Name))

	png(paste(data_dir, "/2covar/graph_age/", ROI_Name, ".png", sep=""), width=500, height=500)
	plot(maledata, type="p", col="blue", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), main=ROI_Name, xlab="Age(day)", ylab="Volume, Area or Thickness")
	par(new=T)
	plot(femaledata, type="p", col="red", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), xlab="", ylab="")

	fit <- lm(meas ~ factor(sex) * centered_age)
#	summary(fit)

	coef <-fit$coefficients

	abline(coef[1]-coef[3]*mean(age), coef[3], col="blue")
	abline(coef[1]+coef[2]-(coef[3]+coef[4])*mean(age), coef[3]+coef[4], col="red")
	abline(coef[1]+femaleratio*coef[2]-(coef[3]+femaleratio*coef[4])*mean(age), coef[3]+femaleratio*coef[4], col="black")
	dev.off()

	###select contrast###
	#null hypothesis: male intercept=0
	#contrast <- matrix(c(1, 0, 0, 0), 1)
	#null hypothesis: female intercept=0
	#contrast <- matrix(c(1, 1, 0, 0), 1)
	#null hypothesis: (female intercept-male intercept)=0
	#contrast <- matrix(c(0, 1, 0, 0), 1)
	#null hypothesis: male slope=0
	#contrast <- matrix(c(0, 0, 1, 0), 1)
	#null hypothesis: female slope=0
	#contrast <- matrix(c(0, 0, 1, 1), 1)
	#null hypothesis: (female slope-male slope)=0
	#contrast <- matrix(c(0, 0, 0, 1), 1)
	#####################

	contrast <- matrix(c(0, 1, 0, 0), 1)
	t.intercept <- glht(fit, linfct = contrast)
	summary.intercept <- summary(t.intercept)$test
	stats.intercept <-c(summary.intercept$coefficients[1],summary.intercept$sigma[1],summary.intercept$tstat[1],summary.intercept$pvalues[1])

	contrast <- matrix(c(0, 0, 0, 1), 1)
	t.slope <- glht(fit, linfct = contrast)
	summary.ageslope <- summary(t.slope)$test
	stats.ageslope <-c(summary.ageslope$coefficients[1],summary.ageslope$sigma[1],summary.ageslope$tstat[1],summary.ageslope$pvalues[1])
	stats <-c(ROI_Name, stats.intercept, stats.ageslope)
	statsdata.2covar<-rbind(statsdata.2covar,stats)
}

statsdata.2covar <- statsdata.2covar[-1,]
#statsdata.2covar
write.csv(statsdata.2covar, paste(data_dir, "/2covar/", output_file_name, sep=""))



########################################
###sex, age, covariate3 as regressors###
########################################

statsdata.3covar<-data.frame(ROI=NA,Intercept_Coefficient=NA, Intercept_Sigma=NA, Intercept_tstat=NA, Intercept_p=NA, Sex_Difference_of_Age_Slope_Coefficient=NA, Sex_Difference_of_Age_Slope_Sigma=NA, Sex_Difference_of_Age_Slope_tstat=NA, Sex_Difference_of_Age_Slope_p=NA, Covariate3_Slope_Coefficient=NA, Covariate3_Slope_Sigma=NA, Covariate3_Slope_tstat=NA, Covariate3_Slope_p=NA)

for (i in (non_ROI_column+1):ncol(data)){
	ROI_Name <- colnames(data)[i]
	meas<-data[,i]

	fit <- lm(meas ~ factor(sex) * (centered_age + centered_covar3))
#	summary(fit)
	coef <-fit$coefficients

	maledata.age<-subset(data, Sex==1, c("Age_at_MRI", ROI_Name))
	femaledata.age<-subset(data, Sex==2, c("Age_at_MRI", ROI_Name))
	png(paste(data_dir, "/3covar/graph_age/", ROI_Name, ".png", sep=""), width=500, height=500)
	plot(maledata.age, type="p", col="blue", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), main=ROI_Name, xlab="Age(day)", ylab="Volume, Area or Thickness")
	par(new=T)
	plot(femaledata.age, type="p", col="red", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), xlab="", ylab="")
	abline(coef[1]-coef[3]*mean(age), coef[3], col="blue")
	abline(coef[1]+coef[2]-(coef[3]+coef[5])*mean(age), coef[3]+coef[5], col="red")
	abline(coef[1]+femaleratio*coef[2]-(coef[3]+femaleratio*coef[5])*mean(age), coef[3]+femaleratio*coef[5], col="black")
	dev.off()

	maledata.covar3<-subset(data, Sex==1, c(covar3_name, ROI_Name))
	femaledata.covar3<-subset(data, Sex==2, c(covar3_name, ROI_Name))
	png(paste(data_dir, "/3covar/graph_covar3/", ROI_Name, ".png", sep=""), width=500, height=500)
	plot(maledata.covar3, type="p", col="blue", xlim=c(min(covar3), max(covar3)), ylim=c(min(meas), max(meas)), main=ROI_Name, xlab=covar3_name, ylab="Volume, Area or Thickness")
	par(new=T)
	plot(femaledata.covar3, type="p", col="red", xlim=c(min(covar3), max(covar3)), ylim=c(min(meas), max(meas)), xlab="", ylab="")
	abline(coef[1]-coef[4]*mean(covar3), coef[4], col="blue")
	abline(coef[1]+coef[2]-(coef[4]+coef[6])*mean(covar3), coef[4]+coef[6], col="red")
	abline(coef[1]+femaleratio*coef[2]-(coef[4]+femaleratio*coef[6])*mean(covar3), coef[4]+femaleratio*coef[6], col="black")
	dev.off()

	###select contrast###
	#null hypothesis: male intercept=0
	#contrast <- matrix(c(1, 0, 0, 0, 0, 0), 1)
	#null hypothesis: female intercept=0
	#contrast <- matrix(c(1, 1, 0, 0, 0, 0), 1)
	#null hypothesis: (female intercept-male intercept)=0
	#contrast <- matrix(c(0, 1, 0, 0, 0, 0), 1)
	#null hypothesis: age slope=0
	#contrast <- matrix(c(0, 0, 1, 0, (nrow(femaledata.age)/nrow(data)), 0), 1)
	#null hypothesis: age slope in male=0
	#contrast <- matrix(c(0, 0, 1, 0, 0, 0), 1)
	#null hypothesis: age slope in female=0
	#contrast <- matrix(c(0, 0, 1, 0, 1, 0), 1)
	#null hypothesis: (age slope in female-age slope in male)=0
	#contrast <- matrix(c(0, 0, 0, 0, 1, 0), 1)
	#null hypothesis: covariate slope=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, (nrow(femaledata.age)/nrow(data))), 1)
	#null hypothesis: covariate slope in male=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, 0), 1)
	#null hypothesis: covariate slope in female=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, 1), 1)
	#null hypothesis: (covariate slope in female-covariate slope in male)=0
	#contrast <- matrix(c(0, 0, 0, 0, 0, 1), 1)
	#####################

	#sex difference of intercept
	contrast <- matrix(c(0, 1, 0, 0, 0, 0), 1)
	t.intercept <- glht(fit, linfct = contrast)
	summary.intercept <- summary(t.intercept)$test
	stats.intercept <-c(summary.intercept$coefficients[1],summary.intercept$sigma[1],summary.intercept$tstat[1],summary.intercept$pvalues[1])

	#sex difference of age slope
	contrast <- matrix(c(0, 0, 0, 0, 1, 0), 1)
	t.ageslope <- glht(fit, linfct = contrast)
	summary.ageslope <- summary(t.ageslope)$test
	stats.ageslope <-c(summary.ageslope$coefficients[1],summary.ageslope$sigma[1],summary.ageslope$tstat[1],summary.ageslope$pvalues[1])

	#covariate slope in both sexes
	contrast <- matrix(c(0, 0, 0, 1, 0, (nrow(femaledata.age)/nrow(data))), 1)
	t.covar3slope <- glht(fit, linfct = contrast)
	summary.covar3slope <- summary(t.covar3slope)$test
	stats.covar3slope <-c(summary.covar3slope$coefficients[1],summary.covar3slope$sigma[1],summary.covar3slope$tstat[1],summary.covar3slope$pvalues[1])

	stats <-c(ROI_Name, stats.intercept, stats.ageslope, stats.covar3slope)
	statsdata.3covar<-rbind(statsdata.3covar,stats)
}

statsdata.3covar <- statsdata.3covar[-1,]
#statsdata.3covar
write.csv(statsdata.3covar, paste(data_dir,"/3covar/", output_file_name, sep=""))

###################################################
###sex, age and 2 other covariates as regressors###
###################################################

statsdata.4covar<-data.frame(ROI=NA,Intercept_Coefficient=NA, Intercept_Sigma=NA, Intercept_tstat=NA, Intercept_p=NA, Sex_Difference_of_Age_Slope_Coefficient=NA, Sex_Difference_of_Age_Slope_Sigma=NA, Sex_Difference_of_Age_Slope_tstat=NA, Sex_Difference_of_Age_Slope_p=NA, Covariate4_Slope_Coefficient=NA, Covariate4_Slope_Sigma=NA, Covariate4_Slope_tstat=NA, Covariate4_Slope_p=NA)

for (i in (non_ROI_column+1):ncol(data)){
	ROI_Name <- colnames(data)[i]
	meas<-data[,i]

	fit <- lm(meas ~ factor(sex) * (centered_age + centered_covar3 + centered_covar4))
#	summary(fit)
	coef <-fit$coefficients

	maledata.age<-subset(data, Sex==1, c("Age_at_MRI", ROI_Name))
	femaledata.age<-subset(data, Sex==2, c("Age_at_MRI", ROI_Name))
	png(paste(data_dir, "/4covar/graph_age/", ROI_Name, ".png", sep=""), width=500, height=500)
	plot(maledata.age, type="p", col="blue", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), main=ROI_Name, xlab="Age(day)", ylab="Volume, Area or Thickness")
	par(new=T)
	plot(femaledata.age, type="p", col="red", xlim=c(min(age), max(age)), ylim=c(min(meas), max(meas)), xlab="", ylab="")
	abline(coef[1]-coef[3]*mean(age), coef[3], col="blue")
	abline(coef[1]+coef[2]-(coef[3]+coef[6])*mean(age), coef[3]+coef[6], col="red")
	abline(coef[1]+femaleratio*coef[2]-(coef[3]+femaleratio*coef[6])*mean(age), coef[3]+femaleratio*coef[6], col="black")
	dev.off()

	maledata.covar4<-subset(data, Sex==1, c(covar4_name, ROI_Name))
	femaledata.covar4<-subset(data, Sex==2, c(covar4_name, ROI_Name))
	png(paste(data_dir, "/4covar/graph_covar4/", ROI_Name, ".png", sep=""), width=500, height=500)
	plot(maledata.covar4, type="p", col="blue", xlim=c(min(covar4), max(covar4)), ylim=c(min(meas), max(meas)), main=ROI_Name, xlab=covar4_name, ylab="Volume, Area or Thickness")
	par(new=T)
	plot(femaledata.covar4, type="p", col="red", xlim=c(min(covar4), max(covar4)), ylim=c(min(meas), max(meas)), xlab="", ylab="")
	abline(coef[1]-coef[5]*mean(covar4), coef[5], col="blue")
	abline(coef[1]+coef[2]-(coef[5]+coef[8])*mean(covar4), coef[5]+coef[8], col="red")
	abline(coef[1]+femaleratio*coef[2]-(coef[5]+femaleratio*coef[8])*mean(covar4), coef[5]+femaleratio*coef[8], col="black")
	dev.off()

	###select contrast###
	#null hypothesis: male intercept=0
	#contrast <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0), 1)
	#null hypothesis: female intercept=0
	#contrast <- matrix(c(1, 1, 0, 0, 0, 0, 0 ,0), 1)
	#null hypothesis: (female intercept-male intercept)=0
	#contrast <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0), 1)
	#null hypothesis: age slope=0
	#contrast <- matrix(c(0, 0, 1, 0, 0, (nrow(femaledata.age)/nrow(data)), 0, 0), 1)
	#null hypothesis: age slope in male=0
	#contrast <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0), 1)
	#null hypothesis: age slope in female=0
	#contrast <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0), 1)
	#null hypothesis: (age slope in female-age slope in male)=0
	#contrast <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0), 1)
	#null hypothesis: covariate slope=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, 0, (nrow(femaledata.age)/nrow(data)), 0), 1)
	#null hypothesis: covariate slope in male=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, 0, 0, 0), 1)
	#null hypothesis: covariate slope in female=0
	#contrast <- matrix(c(0, 0, 0, 1, 0, 0, 1, 0), 1)
	#null hypothesis: (covariate slope in female-covariate slope in male)=0
	#contrast <- matrix(c(0, 0, 0, 0, 0, 0, 1, 0), 1)
	#null hypothesis: 4covar slope=0
	#contrast <- matrix(c(0, 0, 0, 0, 1, 0, 0, (nrow(femaledata.age)/nrow(data))), 1)
	#null hypothesis: 4covar slope in male=0
	#contrast <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0), 1)
	#null hypothesis: 4covar slope in female=0
	#contrast <- matrix(c(0, 0, 0, 0, 1, 0, 0, 1), 1)
	#null hypothesis: (4covar slope in female-covariate slope in male)=0
	#contrast <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1), 1)
	#####################

	#sex difference of intercept
	contrast <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0), 1)
	t.intercept <- glht(fit, linfct = contrast)
	summary.intercept <- summary(t.intercept)$test
	stats.intercept <-c(summary.intercept$coefficients[1],summary.intercept$sigma[1],summary.intercept$tstat[1],summary.intercept$pvalues[1])

	#sex difference of age slope
	contrast <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0), 1)
	t.ageslope <- glht(fit, linfct = contrast)
	summary.ageslope <- summary(t.ageslope)$test
	stats.ageslope <-c(summary.ageslope$coefficients[1],summary.ageslope$sigma[1],summary.ageslope$tstat[1],summary.ageslope$pvalues[1])

	#covar4 slope in both sexes
	contrast <- matrix(c(0, 0, 0, 0, 1, 0, 0, (nrow(femaledata.age)/nrow(data))), 1)
	t.covar4slope <- glht(fit, linfct = contrast)
	summary.covar4slope <- summary(t.covar4slope)$test
	stats.covar4slope <-c(summary.covar4slope$coefficients[1],summary.covar4slope$sigma[1],summary.covar4slope$tstat[1],summary.covar4slope$pvalues[1])

	stats <-c(ROI_Name, stats.intercept, stats.ageslope, stats.covar4slope)
	statsdata.4covar<-rbind(statsdata.4covar,stats)
}

statsdata.4covar <- statsdata.4covar[-1,]
#statsdata.4covar
write.csv(statsdata.4covar, paste(data_dir,"/4covar/", output_file_name, sep=""))