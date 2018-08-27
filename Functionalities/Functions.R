#### Description ####

# R script for common MRI analysis functionalities


#### Parameters ####

common_dir<-file.path(parent_dir,"DropBox/MRI/Statistics/CommonData")
roi_file <- "ROI.csv"


#### Libraries ####
library(tidyverse)


#### Factor to Numeric Converter ####
as.numeric.factor <- function(x) {
  if (class(x)=="factor"){
    return(as.numeric(levels(x))[x])
  }else{
    return(x)
  }
}


#### ID Converter ####

roi_data<-read.csv(file.path(common_dir,roi_file))

ConvertID<-function(input,dict,from_type,to_type){
  from_vec<-as.character(dict[,which(names(dict)==from_type)])
  to_vec<-as.character(dict[,which(names(dict)==to_type)])
  output<-NULL
  for (i in 1:length(input)){
    if(is.na(input[i])){
      # when input=NA, output=NA
      output_add<-NA
    }else if(is.na(which(from_vec==input[i])[1])){
      # when input is not NA but does not exist in the dictionary, output=input
      output_add<-input[i]
    }else if(is.na(to_vec[which(from_vec==input[i])])){
      # when input is not NA and exist in dictionary but corresponding output is NA, output=input
      output_add<-input[i]
    }else{
      # when input is not NA and existin dictionary and corresponding output in not NA, output=corresponeing output
      output_add<-to_vec[which(from_vec==input[i])]
    }
    output<-c(output,output_add)
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
  nodes<-rowid_to_column(nodes, "id")
  input$from<-as.character(input$from)
  input$to<-as.character(input$to)
  edges<-left_join(input[[2]], nodes, by = c("from" = "label")) 
  edges<-edges[,-which(colnames(edges)=="from")]
  edges<-rename(edges, from = id)
  edges<-left_join(edges, nodes, by = c("to" = "label"))
  edges<-edges[,-which(colnames(edges)=="to")]
  edges<-rename(edges, to = id)
  edges<-rename(edges, weight=r)
  edges<-MultCompCorr(edges)
  collabel<-colnames(edges)
  collabel<-collabel[-c(which(collabel=="from"),which(collabel=="to"))]
  collabel<-c("from","to",collabel)
  edges<-edges[,collabel]
  output<-list(nodes,edges)
  names(output)<-c("nodes","edges")
  return(output)
}


#### Transform GLM of FC results into Nodes and Edges ####

GLM_FC2Graph<-function(input_glm,input_nodes){
  nodes<-data.frame(label=as.character(input_nodes))
  nodes$label<-as.character(nodes$label)
  nodes<-rowid_to_column(nodes, "id")
  input_glm$from<-as.character(input_glm$from)
  input_glm$to<-as.character(input_glm$to)
  edges<-left_join(input_glm, nodes, by = c("from" = "label")) 
  edges<-edges[,-which(colnames(edges)=="from")]
  edges<-rename(edges, from = id)
  edges<-left_join(edges, nodes, by = c("to" = "label"))
  edges<-edges[,-which(colnames(edges)=="to")]
  edges<-rename(edges, to = id)
  edges<-rename(edges, weight=t)
  collabel<-colnames(edges)
  collabel<-collabel[-c(which(collabel=="from"),which(collabel=="to"))]
  collabel<-c("from","to",collabel)
  edges<-edges[,collabel]
  output<-list(nodes,edges)
  names(output)<-c("nodes","edges")
  return(output)
}


#### Add columns with multiple comparison corrected p values from column named "p" ####

MultCompCorr<-function(input,n=NULL){
  output<-input
  output$p_Bonferroni<-p.adjust(output$p,method = "bonferroni")
  output$p_Holm_Bonferroni<-p.adjust(output$p,method = "holm")
  output$p_Hochberg<-p.adjust(output$p,method = "hochberg")
  output$p_Hommel<-p.adjust(output$p,method = "hommel")
  output$p_Benjamini_Hochberg<-p.adjust(output$p,method="BH")
  output$p_Benjamini_Yekutieli<-p.adjust(output$p,method="BY")
  if(!is.null(n)){
    output$p_Benjamini_Hochberg_n<-p_BH.adjust(output$p,n)
  }
  return(output)
}


#### BH multiple comparison correcton with n smaller than length(p) ####
p_BH.adjust<-function(p,n=NULL){
  calc<-data.frame(p_input=p)
  calc<-rowid_to_column(calc,"id")
  calc<-calc[order(calc$p_input),]
  calc$rank<-1:nrow(calc)
  if (is.null(n)){
    n_comparisons<-nrow(calc)
  }else{
    n_comparisons<-n
  }
  calc$p_adjusted<-calc$p_input*n_comparisons/calc$rank
  calc[nrow(calc),"p_BH"]<-calc[nrow(calc),"p_adjusted"]
  for (i in (nrow(calc)-1):1){
    calc[i,"p_BH"]<-min(calc[i+1,"p_BH"],calc[i,"p_adjusted"])
  }
  calc<-calc[order(calc$id),]
  output<-calc$p_BH
  return(output)
}