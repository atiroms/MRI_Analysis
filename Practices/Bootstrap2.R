####################
####################
# load this file after StructuralAnalysis.R


#library(mefa)

#duplication_factor<-1:10
n_sampling<-1000

partial_data<-structure_data[,c(2,3)]
#rep_data<-rep(partial_data,each=n_sampling)



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


output3<-data.frame(matrix(ncol=nrow(partial_data)+1, nrow=n_sampling))
output4<-data.frame(matrix(ncol=4,nrow=nrow(partial_data)))
colnames(output3)<-c("sampling", sprintf("%03d_samples_corr",1:nrow(partial_data)))
colnames(output4)<-c("sampling", "mean", "sd", "sd*sqrt(n)")
for (j in 5:nrow(partial_data)){
  for (i in 1:n_sampling){
    output3[i,1]<-i
    sampled_data<-partial_data[sample(nrow(partial_data), size=j,replace=T), ]
    bootstrap_corr<-CalcCorr(sampled_data)[[2]]$cor
    output3[i,j+1]<-bootstrap_corr
  }
#  hist(as.numeric(output3[,j]),breaks=(0.02*(-50:50)))
  output4[j,1]<-j
  output4[j,2]<-mean(as.numeric(output3[,j]))
  output4[j,3]<-standarddeviation(as.numeric(output3[,j]))
  output4[j,4]<-standarddeviation(as.numeric(output3[,j]))*sqrt(j)
#  hist(as.numeric(output3[,j]))
}