#### Description ####

# R script for common MRI analysis functionalities


#### Parameters ####

common_dir<-file.path(parent_dir,"DropBox/MRI/Statistics/CommonData")
roi_file <- "ROI.csv"


#### Libraries ####


#### Factor to Numeric Converter ####
as.numeric.factor <- function(x) {
  as.numeric(levels(x))[x]
}


#### ID Converter ####

roi_data<-read.csv(file.path(common_dir,roi_file))

ConvertID<-function(input,dict,from_type,to_type){
  from_vec<-as.character(dict[,which(names(dict)==from_type)])
  to_vec<-as.character(dict[,which(names(dict)==to_type)])
  output<-NULL
  for (i in 1:length(input)){
    output<-c(output,to_vec[which(from_vec==input[i])])
  }
  return(output)
}


##### Directory Organization ####

if (!file.exists(output_dir)){
  dir.create(file.path(output_dir))
}

ExpDir<-function(exptype){
  timestamp <- strftime(Sys.time(),"%Y%m%d_%H%M%S")
  output<-file.path(output_dir, paste(timestamp,exptype,sep="_"))
  dir.create(output)
#  setwd(output)
  return(output)
}


#### General Correlation Calculation ####

CalcCorr<-function(input, dirname, fileprefix, plot=T,save=T){
  corr <-rcorr(as.matrix(input), type="pearson")
  n_node<-ncol(input)
  corr_flat<-data.frame(matrix(nrow=n_node*(n_node-1)/2,ncol=6))
  colnames(corr_flat)<-c("from","from_label","to","to_label","r","p")
  k<-0
  for (i in 1:(n_node-1)){
    for (j in (i+1):n_node){
      k<-k+1
      corr_flat[k,1:4]<-c(rownames(corr$r)[i],
                       ConvertID(rownames(corr$r)[i],roi_data,"ID_long","label_proper"),
                       colnames(corr$r)[j],
                       ConvertID(colnames(corr$r)[j],roi_data,"ID_long","label_proper"))
      corr_flat[k,5:6]<-c(corr$r[i,j],
                          corr$P[i,j])
    }
  }
  if (plot==T){
    graph<-CorrMatPlot(corr$r,paste(fileprefix,"Correlation Matrix"))
  }else{
    graph<-NULL
  }
  if (save==T){
#    write.csv(corr$r, file.path(dirname,paste(fileprefix,"r.csv",sep="_")))
#    write.csv(corr$P, file.path(dirname,paste(fileprefix,"p.csv",sep="_")))
    write.csv(corr_flat, file.path(dirname,paste(fileprefix,"rp.csv",sep="_")),row.names=F)
  }
  return(list(corr, corr_flat,graph))
}


#### Measure - Clinical data Correlation ####

MeasClinicalCorr<-function(input, dirname){
  n_measures<-ncol(input)-1
  output<-data.frame(matrix(ncol=2*n_clinical_data, nrow=n_measures))
  colnames(output)[2*(1:n_clinical_data)-1]<-paste(colnames(clinical_data)[-1],"r",sep="_")
  colnames(output)[2*(1:n_clinical_data)]<-paste(colnames(clinical_data)[-1],"p",sep="_")
  graph<-NULL
  for (i in 1:n_measures){
    gg<-NULL
    for (j in 1:n_clinical_data){
      cortest<-cor.test(as.numeric(input[,i+1]), clinical_data[,j+1], method="pearson")
      output[i,2*j-1]<-cortest$estimate
      output[i,2*j]<-cortest$p.value
      g<-CorrelationPlot(as.numeric(input[,i+1]),clinical_data[,j+1],colnames(input)[i+1],colnames(clinical_data)[j+1])
      gg<-c(gg,list(g))
    }
    graph<-c(graph,list(gg))
  }
  
  output<-cbind(measure=colnames(input)[-1],output)
  write.csv(output, file.path(dirname,"MeasClinicalCorr.csv"),row.names = F)
  return(list(output,graph))
}


#### Transform Correlation Matrix into Nodes and Edges ####

Corr2Graph<-function(input){
  nodes<-data.frame(label=as.character(colnames(input[[1]]$r)))
  nodes$label<-as.character(nodes$label)
  #  nodes<-rbind(data.frame(label=nodes[(nrow(nodes)/2+1):nrow(nodes),]),data.frame(label=nodes[(nrow(nodes)/2):1,]))
  nodes<-rowid_to_column(nodes, "id")
  edges<-left_join(input[[2]], nodes, by = c("from" = "label")) 
  edges<-rename(edges, from2 = id)
  edges<-left_join(edges, nodes, by = c("to" = "label"))
  edges<-rename(edges, to2 = id)
  #  edges<-rename(edges, r=cor)
  bonferroni<-p_uncorrected/nrow(edges)
  edges$plot<- edges$p<bonferroni
  #  edges<-edges[which(edges$sign==TRUE),]
  edges<-edges[,c("from2","to2","r","p","plot")]
  colnames(edges)<-c("from","to","r","p","plot")
  return(list(nodes,edges))
}