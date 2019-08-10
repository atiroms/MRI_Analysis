library(parallel)

n_core<-detectCores()

clust<-makeCluster(n_core)

func<-function(x){
  return(param*x)
}


list_in<-list(1,2,3,4,5)
param<-4

clusterExport(clust,"param")

parLapply(clust,list_in,func)

parSapply(clust,list_in,func)

stopCluster(clust)
