####################
####################
# load this file after StructuralAnalysis.R

duplication_factor<-1:10

partial_data<-structure_data[,c(2,3)]
duplicated_data<-partial_data


variance <- function(x){
#  var(x)
  var(x)*(length(x)-1)/length(x)
}

standarddeviation<-function(x){
  sqrt(variance(x))
}

FlattenCorr <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

CalcCorr<-function(input){
  corr <-rcorr(as.matrix(input), type="pearson")
  corr_flat <- FlattenCorr(corr$r, corr$P)
  output<-list(corr, corr_flat)
  return(output)
}

n_subject<-nrow(partial_data)
output<-data.frame(matrix(ncol=n_subject+2, nrow=length(duplication_factor)))
colnames(output)<-c("duplication_factor", "overall", sprintf("JK_%03d",1:n_subject))
for (i in 1:length(duplication_factor)){
  output[i,1]<-duplication_factor[i]
  overall_corr<-CalcCorr(duplicated_data)[[2]]$cor
  output[i,2]<-overall_corr
  for (j in 1:n_subject){
    jk_corr<-CalcCorr(duplicated_data[-j,])[[2]]$cor
    output[i,j+2]<-jk_corr
  }
  duplicated_data<-rbind(duplicated_data, partial_data)
}

output2<-data.frame(matrix(ncol=10, nrow=length(duplication_factor)))
jkc<-data.frame(matrix(ncol=n_subject+1, nrow=length(duplication_factor)))
bias<-data.frame(matrix(ncol=n_subject+1, nrow=length(duplication_factor)))
pseudovalue<-data.frame(matrix(ncol=n_subject+1, nrow=length(duplication_factor)))
colnames(output2)<-c("duplication_factor","mean_jkc","sd_jkc","cv_jkc","mean_bias","sd_bias","cv_bias","mean_pseudovalue","sd_pseudovalue","cv_pseudovalue")
colnames(jkc)<-c("duplication_factor", sprintf("jkc_%03d",1:n_subject))
colnames(bias)<-c("duplication_factor", sprintf("bias_%03d",1:n_subject))
colnames(pseudovalue)<-c("duplication_factor", sprintf("pseudovalue_%03d",1:n_subject))
for (i in 1:length(duplication_factor)){
  output2[i,1]<-duplication_factor[i]
  jkc[i,1]<-duplication_factor[i]
  bias[i,1]<-duplication_factor[i]
  pseudovalue[i,1]<-duplication_factor[i]
  jkc_individual<-as.numeric(output[i,c(-1,-2)])
  mean_jkc<-mean(jkc_individual)
  sd_jkc<-standarddeviation(jkc_individual)
  cv_jkc<-sd_jkc/mean_jkc
  bias_individual<-as.numeric(output[i,2]-output[i,3:(n_subject+2)])
  #ordinary pseudovalue
  pseudovalue_individual<-as.numeric((n_subject*i)*output[i,2]-(i*n_subject-1)*output[i,3:(n_subject+2)])
  #sqrt-modified version of pseudovalue
#  pseudovalue_individual<-as.numeric((1+sqrt(i*n_subject))*output[i,2]-(sqrt(i*n_subject))*output[i,3:(n_subject+2)])
  
  mean_bias<-mean(bias_individual)
  sd_bias<-standarddeviation(bias_individual)
  cv_bias<-sd_bias/mean_bias
  mean_pseudovalue<-mean(pseudovalue_individual)
  sd_pseudovalue<-standarddeviation(pseudovalue_individual)
  cv_pseudovalue<-sd_pseudovalue/mean_pseudovalue
  
  output2[i,2:10]<-c(mean_jkc,sd_jkc,cv_jkc,mean_bias,sd_bias,cv_bias,mean_pseudovalue,sd_pseudovalue,cv_pseudovalue)
  jkc[i,2:(1+n_subject)]<-jkc_individual
  bias[i,2:(1+n_subject)]<-bias_individual
  pseudovalue[i,2:(1+n_subject)]<-pseudovalue_individual
  
}

for (i in 1:length(duplication_factor)){
  hist(as.numeric(jkc[i,-1]), breaks=(0.001*(10:50)))
}

for (i in 1:length(duplication_factor)){
  hist(as.numeric(bias[i,-1]),breaks=(0.001*(-20:15)))
}

for (i in 1:length(duplication_factor)){
  hist(as.numeric(pseudovalue[i,-1]),breaks=(0.1*(-50:50)))
#  hist(as.numeric(pseudovalue[i,-1]),breaks=(0.01*(-30:30)))
}