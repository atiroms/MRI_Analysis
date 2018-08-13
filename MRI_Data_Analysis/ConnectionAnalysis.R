#### Description ####

# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.
# 


#### Parameters ####

working_dir<-"G:/MRI/Statistics/Connection"
#working_dir<-"D:/MRI/Statistics/Connection"

#structural_inputfile <- "182_FS_DK_Thickness_SCA.csv"

#connection_file <- "HO_FC_r_TSexist_male.csv"
#connection_file <- "HO_FC_r_TSexist_female.csv"
#connection_file <- "Power_FC_r_TSexist_male.csv"
connection_file <- "Power_FC_r_TSexist_female.csv"

data_subdir<-"Connection_data"

commondata_dir<-"G:/MRI/Statistics/CommonData"
#commondata_dir<-"D:/MRI/Statistics/CommonData"

roiid_file<-"ROI_All.csv"
input_roiid_type<-"ID_long"
graph_roiid_type<-"label_proper"
clinical_file <- "CSUB_Clinical_Data.csv"

p_uncorrected<-0.001

#n_components<-10
#n_components<-30
n_components<-5
tsne_dims<-2
tsne_perplexity<-30
tsne_max_itr<-1000

cost<-seq(0.04,0.25,0.01)


#### Libraries ####
library(multcomp)
library(FactoMineR)
library(ica)
library(Rtsne)
library(ggpubr)
library(RColorBrewer)
library(igraph)
library(qgraph)


#### Data Loading ####
#structural_data<-read.csv(file.path(working_dir, structural_inputfile))
connection_data<-read.csv(file.path(working_dir, connection_file))
connections<-colnames(connection_data)[-1]
n_connections<-length(connections)
subject_id<-unique(connection_data$ID_pnTTC)
n_subject<-length(subject_id)

connections_matrix<-data.frame(matrix(nrow=n_connections,ncol=2))
for (i in 1:n_connections){
  split<-strsplit(connections[i], "_")
  connections_matrix[i,1]<-split[[1]][[1]]
  connections_matrix[i,2]<-split[[1]][[2]]
}
regions<-data.frame(label=c(connections_matrix[,1],connections_matrix[,2]))
regions<-as.character(unique(regions[,1]))
regions<-regions[order(regions)]

clinical_data <-read.csv(file.path(commondata_dir,clinical_file))
clinical_row<-NULL
for (i in 1:n_subject){
  clinical_row<-c(clinical_row,which(clinical_data$ID_pnTTC==subject_id[i]))
}
clinical_data<-clinical_data[clinical_row,]
rownames(clinical_data)<-NULL


#### Directory Organization ####

if (!file.exists(file.path(working_dir,data_subdir))){
  dir.create(file.path(working_dir, data_subdir))
}

ExpDir<-function(exptype){
  timestamp <- strftime(Sys.time(),"%Y%m%d_%H%M%S")
  data_dir<-file.path(working_dir, data_subdir, paste(timestamp,exptype,sep="_"))
  dir.create(data_dir)
  return(data_dir)
}


#### ID Converter ####

roiid_data<-read.csv(file.path(commondata_dir,roiid_file))
roiid_data$ID_long<-as.character(roiid_data$ID_long)
roiid_data$label<-as.character(roiid_data$label)
roiid_data$label_proper<-as.character(roiid_data$label_proper)
roiid_data$label_long<-as.character(roiid_data$label_long)
roiid_data$label_conn<-as.character(roiid_data$label_conn)
roiid_data$label_fs<-as.character(roiid_data$label_fs)

ConvertID<-function(input,iddata,inputcolumn,outputcolumn){
  idcolname<-data.frame(label=colnames(iddata))
  input_col<-iddata[,which(idcolname$label==inputcolumn)]
  output_col<-iddata[,which(idcolname$label==outputcolumn)]
  output<-data.frame(matrix(nrow=length(input),ncol=1))
  for (i in 1:length(input)){
    input_value<-input[[i]]
    output_row<-which(input_col==input_value)
    output[i,1]<-output_col[output_row]
  }
  output<-output[,1]
  return(output)
}


#### Component-Clinical Plotting ####

ComponentPlot<-function(input, title){
  graph<-NULL
  for (i in 2:ncol(clinical_data)){
    gg<-NULL
    for (j in 1:ncol(input)){
      g<-ggscatter(cbind(clinical_data,input), x = sprintf("Dim_%02d",j), y = colnames(clinical_data)[i],
                   add = "reg.line", conf.int = TRUE, size=0.5,
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = paste("Dimension",j,sep=" "), ylab = colnames(clinical_data)[i])
      #      g<-ggplot(cbind(clinical_data,input),aes(x=sprintf("Dim_%02d",j), y=colnames(clinical_data)[i])) +
      #        geom_point() +
      #        geom_smooth(method = "lm") +
      #        geom_rug() +
      #        ggtitle(paste(colnames(clinical_data)[i],title,"Dim",j,sep=" ")) +
      #        xlab(paste("Dimension",j,sep=" ")) +
      #        ylab(colnames(clinical_data)[i]) +
      #        theme_classic()
      gg<-c(gg,list(g))
    }
    graph<-c(graph,list(gg))
  }
  return(graph)
}


#### Two Component-Clinical Color Plotting ####

TwoComponentPlot<-function(input, title){
  for (j in 1:(ncol(input)-1)){
    dims<-c(j,j+1)
    plot(input[,dims],pch=20,main=title,xlab=paste("Dim",j,sep=" "), ylab=paste("Dim",(j+1),sep=" "))
    text(as.numeric(input[,dims[1]]), as.numeric(input[,dims[2]]),labels=subject_id, cex= 0.7, pos=3, offset=0.2)
  }
  for (i in 2:ncol(clinical_data)){
    colors <- brewer.pal(10, "Spectral")[as.integer(cut(as.vector(clinical_data[[i]]), 10))]
    for (j in 1:(ncol(input)-1)){
      dims<-c(j,j+1)
      plot(input[,dims], pch=20, main=paste(title,colnames(clinical_data)[i],sep=" "), col=colors,xlab=paste("Dim",j,sep=" "), ylab=paste("Dim",(j+1),sep=" "))
    }
  }
}


#### Component - Clinical data Correlation ####

ComponentClinicalCorr<-function(input, dirname){
  n_clinical_data<-ncol(clinical_data)-1
  output<-data.frame(matrix(ncol=2*n_clinical_data, nrow=n_components))
  colnames(output)[2*(1:n_clinical_data)-1]<-paste(colnames(clinical_data)[-1],"r",sep="_")
  colnames(output)[2*(1:n_clinical_data)]<-paste(colnames(clinical_data)[-1],"p",sep="_")
  rownames(output)<-sprintf("Dim_%02d",1:n_components)
  for (i in 1:n_components){
    for (j in 1:n_clinical_data){
      cortest<-cor.test(as.numeric(input[,i]), clinical_data[,j+1], method="pearson")
      output[i,2*j-1]<-cortest$estimate
      output[i,2*j]<-cortest$p.value
    }
  }
  write.csv(output, file.path(dirname,"ComponentClinical.csv"))
  return(output)
}


#### GLM Analysis  ####

DoGLM_FC<-function(){
  dirname<-ExpDir("GLMClinical")
  output<-data.frame(matrix(ncol=3, nrow=n_connections))
  output[,1]<-connections
  output[,2]<-ConvertID(connections_matrix[,1],roiid_data,input_roiid_type,"label_proper")
  output[,3]<-ConvertID(connections_matrix[,2],roiid_data,input_roiid_type,"label_proper")
  colnames(output)<-c("Connection_ID","row","column")
  
  TS<-clinical_data$Tanner_Stage
  age<-clinical_data$Age_at_MRI
  centered_TS<-TS-mean(TS)
  centered_age<-age-mean(age)

  output_add<-data.frame(matrix(ncol=6,nrow=n_connections))
  
  for (i in 1:n_connections){
    ObjV<-connection_data[,i+1]
    glmfit <- lm(ObjV ~ centered_TS + centered_age)
    contrast <- matrix(c(0, 1, 0), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("TS_Age_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  for (i in 1:n_connections){
    ObjV<-connection_data[,i+1]
    glmfit <- lm(ObjV ~ centered_TS)
    contrast <- matrix(c(0, 1), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("TS",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("TS_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  for (i in 1:n_connections){
    ObjV<-connection_data[,i+1]
    glmfit <- lm(ObjV ~ centered_age)
    contrast <- matrix(c(0, 1), 1)
    ttest <- summary(glht(glmfit, linfct = contrast))$test
    stats <-c(ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],AIC(glmfit),BIC(glmfit))
    output_add[i,]<-stats
    collabel<-paste("Age",c("beta","sigma","t","p"),sep="_")
    collabel<-c(collabel,"AIC","BIC")
    collabel<-paste("Age_model",collabel,sep="_")
    colnames(output_add)<-collabel
  }
  output<-cbind(output,output_add)
  
  
  best_model<-data.frame(matrix(ncol=1, nrow=n_connections))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (k in 1:n_connections){
      best_id<-which.min(xic[k,])
      if (best_id==1){
        best_model[k,1]<-"TS_Age"
      }else if (best_id==2){
        best_model[k,1]<-"TS"
      }else if (best_id==3){
        best_model[k,1]<-"Age"
      }
    }
    colnames(best_model)<-paste(i,"best_model",sep="_")
    output<-cbind(output,best_model)
  }
  write.csv(output, file.path(dirname,"GLMresult.csv"),row.names=F)
  return(output)
}


#### Principal Component Analysis ####

DoPCA_FC<-function(){
  dirname<-ExpDir("PCA_FC")
  data<-connection_data[-1]
  indexcolumn<-connection_data[1]
  data<-data.matrix(data)
  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
  varfactor<-data.frame(pca$var$coord)
  varfactor<-cbind(connections,ConvertID(connections_matrix[,1],roiid_data,input_roiid_type,"label_proper"),ConvertID(connections_matrix[,2],roiid_data,input_roiid_type,"label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(pca$ind$coord)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  pairs(indfactor)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-pca$eig
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  ComponentClinicalCorr(indfactor[-1],dirname)
  graph<-ComponentPlot(indfactor[-1],"PCA")
  return(list(varfactor,indfactor,varianceaccounted,graph))
}


#### Independent Component Analysis ####

DoICA_FC<-function(){
  dirname<-ExpDir("ICA_FC")
  data<-connection_data[-1]
  indexcolumn<-connection_data[1]
  data<-data.matrix(data)
  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
  varfactor<-data.frame(ica$M)
  varfactor<-cbind(connections,ConvertID(connections_matrix[,1],roiid_data,input_roiid_type,"label_proper"),ConvertID(connections_matrix[,2],roiid_data,input_roiid_type,"label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(ica$S)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  pairs(indfactor)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-ica$vafs
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  ComponentClinicalCorr(indfactor[-1],dirname)
  graph<-ComponentPlot(indfactor[-1],"ICA")
  
  return(list(varfactor,indfactor,varianceaccounted,graph))
}


#### t-SNE Analysis ####

DoTSNE_FC<-function(){
  dirname<-ExpDir("tSNE_FC")
  indexcolumn<-connection_data[1]
  data<-data.matrix(connection_data[-1])
  tsne_fc<-Rtsne(data, dims = tsne_dims, perplexity=tsne_perplexity, verbose=TRUE, max_iter = tsne_max_itr)
  
  indfactor<-data.frame(tsne_fc$Y)
  colnames(indfactor)<-sprintf("Dim_%02d",1:tsne_dims)
  indfactor<-cbind(indexcolumn,indfactor)
  

  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  graph<-TwoComponentPlot(indfactor[-1], "t-SNE_FC")
  return(list(tsne_fc,graph))
}


#### Graph object creation ####

Edges2Graph<-function(input){
  edgename<-colnames(input)[-1]
  n_inputedges<-length(edgename)
  edges<-data.frame(matrix(ncol=3, nrow=n_inputedges))
  for (i in 1:n_inputedges){
    split<-strsplit(edgename[[i]], "_")
    edges[i,1]<-split[[1]][[1]]
    edges[i,2]<-split[[1]][[2]]
  }
  edges[,3]<-as.numeric(input[1,-1])
  colnames(edges)<-c("from","to","weight")
  nodes<-c(edges[,1],edges[,2])
  nodes<-data.frame(id=unique(nodes))
  nodes$id<-nodes$id[order(nodes$id)]
  nodes$label_proper<-ConvertID(nodes$id,roiid_data,"ID_long","label_proper")
  
  output <- graph.data.frame(d = edges, vertices = nodes, directed = F)
  
  return(output)
}


#### Subset edges according to desired cost ####

SubsetEdges<-function(input_graph, input_cost){
  n_nodes<-length(V(input_graph))
  n_edges<-n_nodes*(n_nodes-1)/2
  n_edges4cost<-as.integer(n_nodes*(n_nodes-1)/2*input_cost)
  edges2delete<-head(order(E(input_graph)$weight),(n_edges-n_edges4cost))
  output<-delete.edges(input_graph,edges2delete)
  return(output)
}


#### Iterate over costs ####

ItrCost<-function(input_graph, input_cost=cost){
  n_nodes<-length(V(input_graph))
  output<-data.frame()
  for (i in input_cost){
    subgraph<-SubsetEdges(input_graph,i)
    E(subgraph)$weight<-1
    ## graph-level metrics
    metric<-data.frame(matrix(nrow=1,ncol=0))
    # characteristic path length
    metric<-cbind(metric,averagepathlength=average.path.length(subgraph))
    # global efficiency
    eff<-1/(shortest.paths(subgraph))
    eff[!is.finite(eff)]<-0
    metric<-cbind(metric,globalefficiency=mean(eff,na.rm=TRUE))
    # global clustering coefficient
    metric<-cbind(metric,globalclusteringcoefficient=transitivity(subgraph))
    # average clustering coefficient
    metric<-cbind(metric,averageclusteringcoefficient=transitivity(subgraph,type="average"))
    # local efficiency
    # modularity
    # small-worldness
    suppressWarnings(metric<-cbind(metric,smallworldindex=smallworldIndex(subgraph)$index))
    
    ## node-level metrics
    # degree centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-centr_degree(subgraph)$res
    colnames(nodemetric)<-paste("degreecentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    # betweenness centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-centr_betw(subgraph)$res
    colnames(nodemetric)<-paste("betweennesscentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    # eigenvector centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-eigen_centrality(subgraph)$vector
    colnames(nodemetric)<-paste("eigenvectorcentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    
    output<-rbind(output,metric)
  }
  output<-data.frame(cost=input_cost,output)
  
  average<-data.frame(matrix(ncol=ncol(output),nrow=1))
  colnames(average)<-colnames(output)
  for (i in 2:ncol(output)){
    average[1,i]<-mean(output[,i])
  }
  output<-rbind(average,output)
  return(output)
}


#### Graph Theoretical Analysis, subject-wise ####

DoGTA<-function(){
  dirname<-ExpDir("GTA")
  output<-data.frame()
  for (i in 1:n_subject){
    subject_graph<-Edges2Graph(connection_data[i,])
    subject_metric<-ItrCost(subject_graph)
    output<-rbind(output,subject_metric[1,-1])
  }
  output<-data.frame(ID_pnTTC=subject_id,output)
  write.csv(output, file.path(dirname,"GTA.csv"),row.names=F)
  return(output)
}


#### Circular Plotting ####

CircularPlot<-function(input){
  nodes<-input[[1]]
  nodes$angle <- 90 - 360 * nodes$id / nrow(nodes)
  nodes$hjust<-ifelse(nodes$angle < -90, 1, 0)
  nodes$angle<-ifelse(nodes$angle < -90, nodes$angle+180, nodes$angle)
  nodes$label_proper<-ConvertID(as.character(nodes$label),roiid,input_roiid_type,graph_roiid_type)
  graph_data <- graph_from_data_frame(d = input[[2]], vertices = nodes, directed = TRUE)
  output<-ggraph(graph_data, layout = "linear",circular = TRUE) + 
    geom_edge_arc(aes(color=r, alpha=r),width=2) +
    geom_node_text(aes(x = x*1.03, y=y*1.03, label=label_proper, angle = angle, hjust=hjust), size=3, alpha=1) +
    geom_node_point(aes(x=x, y=y),size=2, alpha=1) +
    scale_edge_color_continuous(low = "green", high = "red")+
    scale_edge_alpha_continuous(range = c(0.1, 1))+
    expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
    theme_void()
  return(output)
  #    theme_graph(background=NA)
  #    geom_edge_arc(aes(color=cor, width=2),alpha = 0.2) + 
  #    scale_edge_width(range = c(0.2, 2)) +
  #    geom_node_text(aes(label = label)) +
  #    labs(edge_width = "Letters") +
  #    scale_colour_manual(values= rep(brewer.pal(9,"Paired") , 30)) +
  #    scale_size_continuous(range = c(0.1,10) ) +
  #    scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  #    theme(legend.position="none", plot.margin=unit(c(0,0,0,0),"cm")) +
  #    theme_graph(background=NA, caption_size = NA)
  #    theme_void()
}

