library(parallel)

n_core<-detectCores()

clust<-makeCluster(n_core)

func<-function(x){
  return(c(x^2,x^3,param*x))
}


list_in<-1:5
param<-4

clusterExport(clust,"param")
parLapply(clust,list_in,func)

stopCluster(clust)
