#### Description ####

# R script to analyze Structural and Functional data
# 
# 

#### Parameters ####

working_dir<-"G:/MRI/Statistics/MultiModal"
#working_dir<-"G:/MRI/Statistics/Functional_CONN_HO"
#working_dir<-"G:/MRI/Statistics/Structural_FS"

#for SCA_FC_All
structural_inputfile <- "182_FS_DK_Thickness_SCA.csv"
#structural_inputfile <- "182_FS_DK_Volume_SCA.csv"
#structural_inputfile <- "182_FS_DK_Area_SCA.csv"
functional_inputfile <- "182_CONN_DK_FC_All.csv"

#for JKZ_FCZ
#structural_inputfile <- "182_FS_DK_Thickness_JKZ.csv"
#structural_inputfile <- "182_FS_DK_Volume_JKZ.csv"
#structural_inputfile <- "182_FS_DK_Area_JKZ.csv"
#functional_inputfile <- "182_CONN_DK_FC.csv"

#for JKPV_FC
#structural_inputfile <- "182_FS_DK_Thickness_JKPV.csv"
#structural_inputfile <- "182_FS_DK_Volume_JKPV.csv"
#structural_inputfile <- "182_FS_DK_Area_JKPV.csv"
#functional_inputfile <- "182_CONN_DK_FC.csv"

#others
#structural_inputfile <- "Freesurfer_CorticalArea.csv"
#functional_inputfile <- "CONN_ROI_Timeseries_Power.csv"
#functional_inputfile <- "CONN_ROI_Timeseries_HO.csv"

data_subdir<-"MultiModal_data"

#subjid_type<-"ID_T1QC_rsfMRIexist_CONNvoxelQC10"
subjid_type<-"ID_T1QC_rsfMRIexist_CONNvoxelQC20"
#subjid_type<-"ID_T1QC_rsfMRIexist"
#subjid_type<-"ID_T1QC"
#subjid_type<-"ID_pnTTC"


commondata_dir<-"G:/MRI/Statistics/CommonData"

subjid_file<-"CSUB_ID.csv"

roiid_file<-"ROI_All.csv"
input_roiid_type<-"label_fs"
graph_roiid_type<-"label_proper"
clinical_file <- "CSUB_Clinical_Data.csv"

p_uncorrected<-0.001


#### Libraries ####
library(ggpubr)
library(tidyverse)
library(ggraph)
library(igraph)


#### Data Loading ####
structural_data<-read.csv(file.path(working_dir, structural_inputfile))
functional_data<-read.csv(file.path(working_dir, functional_inputfile))

subjid_data<-read.csv(file.path(commondata_dir,subjid_file))

roiid<-read.csv(file.path(commondata_dir,roiid_file))
roiid$ID_long<-as.character(roiid$ID_long)
roiid$label<-as.character(roiid$label)
roiid$label_proper<-as.character(roiid$label_proper)
roiid$label_long<-as.character(roiid$label_long)
roiid$label_conn<-as.character(roiid$label_conn)
roiid$label_fs<-as.character(roiid$label_fs)



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


#### Subject ID Extraction ####

subjidcolname<-data.frame(label=colnames(subjid_data))
n_subject<-max(na.omit(subjid_data[,which(subjidcolname$label==subjid_type)]))
subject_id<-ConvertID(1:n_subject,subjid_data,subjid_type,"ID_pnTTC")

#### Data Subsetting According to Subject ID file ####

SubsetSubj<-function(input){
  dirname<-ExpDir("SubjectSubset")
  selectedrow<-numeric()
  for (i in 1:n_subject){
    selectedrow<-c(selectedrow,which(input$ID_pnTTC==subject_id[i]))
  }
  output<-input[selectedrow,]
  write.csv(output,file.path(dirname,"SubsetSubject.csv"),row.names=F)
  return(output)
}


#### Data stratification ####

TransposeFC<-function(input){
#  dirname<-ExpDir("SubjectSubset")
  output<-data.frame(matrix(nrow=2,ncol=nrow(input)))
  colnames(output)<-paste(as.character(input$row),as.character(input$column),sep="_")
  rownames(output)<-colnames(input)[3:4]
  output[1,]<-input[,3]
  output[2,]<-input[,4]
#  write.csv(output,file.path(dirname,output_file),row.names=F)
  return(output)  
}


#### Data Subsetting According to Multimodal ROI Intersection ####

IntersectROI<-function(structural_input=structural_data,functional_input=functional_data){
  structural_edge<-colnames(structural_input)[-1]
  functional_edge<-colnames(functional_input)[-1]
  intersect_edge<-intersect(structural_edge,functional_edge)
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


Corr2Graph<-function(input_node, input_edge, input_corr){
  nodes<-data.frame(label=input_node)
  nodes<-rowid_to_column(nodes, "id")
  edges<-left_join(input_edge, nodes, by = c("row" = "label")) 
  edges<-rename(edges, from = id)
  edges<-left_join(edges, nodes, by = c("column" = "label"))
  edges<-rename(edges, to = id)
  #  edges<-rename(edges, r=cor)
  bonferroni<-p_uncorrected/nrow(edges)
  edges$sign<- edges$p<bonferroni
  edges<-edges[which(edges$sign==TRUE),]
  edges<-select(edges, from, to, r, p)
  return(list(nodes,edges))
}


StrFuncGraphs<-function(input_edge, input_structure, input_function){
  intersect_edge<-data.frame(matrix(ncol=2, nrow=length(input_edge)))
  colnames(input_edge)<-c("row","column")
  for (i in 1:length(input_edge)){
    split<-strsplit(inputedge[[i]], "_")
    intersect_edge[i,1]<-split[[1]][[1]]
    intersect_edge[i,2]<-split[[1]][[2]]
  }
  intersect_node<-data.frame(label=c(intersect_edge[,1],intersect_edge[,2]))
  intersect_node<-as.character(unique(intersect_node[,1]))
  intersect_node<-intersect_node[order(intersect_node)]
  
  structural_subset_edge<-cbind(intersect_edge,t(input_structure))
  rownames(structural_subset_edge)<-NULL
  structural_graph<-CircularPlot(Corr2Graph(intersect_node, structural_subset_edge))
  functional_subset_edge<-cbind(intersect_edge,t(input_function))
  rownames(functional_subset_edge)<-NULL
  functional_graph<-CircularPlot(Corr2Graph(intersect_node, functional_subset_edge))
  
  return(list(structural_graph,functional_graph))
}


#### SCA and FC_All ####

SCA_FC_All<-function(){
  dirname<-ExpDir("SCA_FC_All")
  structural_input<-TransposeFC(structural_data)
  functional_input<-TransposeFC(functional_data)
  structural_edge<-colnames(structural_input)
  functional_edge<-colnames(functional_input)
  intersect_edge<-intersect(structural_edge,functional_edge)
  structural_edgeid<-functional_edgeid<-numeric()
  for (i in 1:length(intersect_edge)){
    structural_edgeid<-c(structural_edgeid,which(structural_edge==intersect_edge[i]))
    functional_edgeid<-c(functional_edgeid,which(functional_edge==intersect_edge[i]))
  }
  structural_subset<-structural_input[,structural_edgeid]
  structural_subset_r<-as.numeric(structural_subset[1,])
  functional_subset<-functional_input[,functional_edgeid]
  functional_subset_r<-as.numeric(functional_subset[1,])
  s_f_r<-data.frame(edge=intersect_edge,structural=structural_subset_r,functional=functional_subset_r)
  corr<-cor.test(structural_subset_r,functional_subset_r,method="pearson")
  output_corr<-data.frame(r=corr$estimate,t=corr$statistic,p=corr$p.value)
  write.csv(s_f_r,file.path(dirname,"SCA_FC_All_Edges.csv"),row.names=F)
  write.csv(output_corr,file.path(dirname,"SCA_FC_All_Statistic.csv"),row.names=F)
  
  strfuncgraphs<-StrFuncGraphs(intersect_edge,structural_subset,functional_subset)
  
  multimodal_graph<-ggscatter(s_f_r, x = "structural", y = "functional",
                   add = "reg.line", conf.int = TRUE, size=0.5,
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Structural Covariance", ylab = "Functional Correlation")
  
  return(list(s_f_r,corr,strfuncraphs[[1]],strfuncgraphs[[2]], multimodal_graph))
}


#### Subject-wise Correlations of z-score of Jackknife Parameter Estimate and of Functional Correlation ####

JKZ_FCZ_Subject<-function(){
  dirname<-ExpDir("JKZ_FCZ_Subject")
  structural_edge<-colnames(structural_data)[-1]
  functional_edge<-colnames(functional_data)[-1]
  intersect_edge<-intersect(structural_edge,functional_edge)
  structural_edgeid<-functional_edgeid<-numeric()
  for (i in 1:length(intersect_edge)){
    structural_edgeid<-c(structural_edgeid,which(structural_edge==intersect_edge[i]))
    functional_edgeid<-c(functional_edgeid,which(functional_edge==intersect_edge[i]))
  }
  structural_subset<-structural_data[,structural_edgeid]
  functional_subset<-functional_data[,functional_edgeid]
  functional_zscore<-functional_subset
  for (i in 1:length(intersect_edge)){
    group_mean<-mean(as.numeric(functional_subset[,i]))
    group_sd<-sd(as.numeric(functional_subset[,i]))
    functional_zscore[,i]<-(as.numeric(functional_zscore[,i])-group_mean)/group_sd
  }
  subject_wise_s_f_r<-data.frame(matrix(ncol=4,nrow=n_subject))
  colnames(subject_wise_s_f_r)<-c("ID_pnTTC","r","t","p")
  subject_wise_s_f_r[,1]<-subject_id
  for (j in 1:n_subject){
    s_f_r<-data.frame(structural=as.numeric(structural_subset[j,]),functional=as.numeric(functional_zscore[j,]))
    cortest<-cor.test(as.numeric(s_f_r$structural),as.numeric(s_f_r$functional),method="pearson")
    r<-as.numeric(cortest$estimate)
    t<-as.numeric(cortest$statistic)
    p<-as.numeric(cortest$p.value)
    subject_wise_s_f_r[j,2:4]<-c(r,t,p)
    ggscatter(s_f_r, x = "structural", y = "functional", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "z-score of Structural Jackknife Estimate", ylab = "z-score of Functional Correlation")
  }
  write.csv(subject_wise_s_f_r,file.path(dirname,"JKZ_FCZ_Subject.csv"),row.names=F)
  return(subject_wise_s_f_r)
}


#### Edge-wise Correlations of z-score of Jackknife Parameter Estimate and of Functional Correlation ####

JKZ_FCZ_Edge<-function(){
  dirname<-ExpDir("JKZ_FCZ_Edge")
  structural_edge<-colnames(structural_data)[-1]
  functional_edge<-colnames(functional_data)[-1]
  intersect_edge<-intersect(structural_edge,functional_edge)
  n_edge<-length(intersect_edge)
  structural_edgeid<-functional_edgeid<-numeric()
  for (i in 1:n_edge){
    structural_edgeid<-c(structural_edgeid,which(structural_edge==intersect_edge[i]))
    functional_edgeid<-c(functional_edgeid,which(functional_edge==intersect_edge[i]))
  }
  structural_subset<-structural_data[,structural_edgeid]
  functional_subset<-functional_data[,functional_edgeid]
  functional_zscore<-functional_subset
  for (i in 1:n_edge){
    group_mean<-mean(as.numeric(functional_subset[,i]))
    group_sd<-sd(as.numeric(functional_subset[,i]))
    functional_zscore[,i]<-(as.numeric(functional_zscore[,i])-group_mean)/group_sd
  }
  edge_wise_s_f_r<-data.frame(matrix(ncol=4,nrow=n_edge))
  colnames(edge_wise_s_f_r)<-c("edge","r","t","p")
  edge_wise_s_f_r[,1]<-intersect_edge
  for (j in 1:n_edge){
    s_f_r<-data.frame(structural=as.numeric(structural_subset[,j]),functional=as.numeric(functional_zscore[,j]))
    cortest<-cor.test(as.numeric(s_f_r$structural),as.numeric(s_f_r$functional),method="pearson")
    r<-as.numeric(cortest$estimate)
    t<-as.numeric(cortest$statistic)
    p<-as.numeric(cortest$p.value)
    edge_wise_s_f_r[j,2:4]<-c(r,t,p)
#    ggscatter(s_f_r, x = "structural", y = "functional", 
#              add = "reg.line", conf.int = TRUE, 
#              cor.coef = TRUE, cor.method = "pearson",
#              xlab = "z-score of Structural Jackknife Estimate", ylab = "z-score of Functional Correlation")
  }
  write.csv(edge_wise_s_f_r,file.path(dirname,"JKZ_FCZ_Edge.csv"),row.names=F)
  return(edge_wise_s_f_r)
}


#### Subject-wise Corrlations of Jackknife Pseudovalue and Functional Correlation ####

JKPV_FC_Subject<-function(){
  dirname<-ExpDir("JKPV_FC_Subject")
  structural_edge<-colnames(structural_data)[-1]
  functional_edge<-colnames(functional_data)[-1]
  intersect_edge<-intersect(structural_edge,functional_edge)
  structural_edgeid<-functional_edgeid<-numeric()
  for (i in 1:n_subject){
    structural_edgeid<-c(structural_edgeid,which(structural_edge==intersect_edge[i]))
    functional_edgeid<-c(functional_edgeid,which(functional_edge==intersect_edge[i]))
  }
  structural_subset<-structural_data[,structural_edgeid]
  functional_subset<-functional_data[,functional_edgeid]
  subject_wise_s_f_r<-data.frame(matrix(ncol=4,nrow=n_subject))
  colnames(subject_wise_s_f_r)<-c("ID_pnTTC","r","t","p")
  subject_wise_s_f_r[,1]<-subject_id
  for (j in 1:n_subject){
    s_f_r<-data.frame(structural=as.numeric(structural_subset[j,]),functional=as.numeric(functional_subset[j,]))
    cortest<-cor.test(as.numeric(s_f_r$structural),as.numeric(s_f_r$functional),method="pearson")
    r<-as.numeric(cortest$estimate)
    t<-as.numeric(cortest$statistic)
    p<-as.numeric(cortest$p.value)
    subject_wise_s_f_r[j,2:4]<-c(r,t,p)
#    ggscatter(s_f_r, x = "structural", y = "functional", 
#              add = "reg.line", conf.int = TRUE, 
#              cor.coef = TRUE, cor.method = "pearson",
#              xlab = "Structural Jackknife Pseudovalue", ylab = "Functional Correlation")
  }
  write.csv(subject_wise_s_f_r,file.path(dirname,"JKPV_FC_Subject.csv"),row.names=F)
  return(subject_wise_s_f_r)
}


#### Edge-wise Corrlations of Jackknife Pseudovalue and Functional Correlation ####

JKPV_FC_Edge<-function(){
  dirname<-ExpDir("JKPV_FC_Edge")
  structural_edge<-colnames(structural_data)[-1]
  functional_edge<-colnames(functional_data)[-1]
  intersect_edge<-intersect(structural_edge,functional_edge)
  n_edge<-length(intersect_edge)
  structural_edgeid<-functional_edgeid<-numeric()
  for (i in 1:n_edge){
    structural_edgeid<-c(structural_edgeid,which(structural_edge==intersect_edge[i]))
    functional_edgeid<-c(functional_edgeid,which(functional_edge==intersect_edge[i]))
  }
  structural_subset<-structural_data[,structural_edgeid]
  functional_subset<-functional_data[,functional_edgeid]
  edge_wise_s_f_r<-data.frame(matrix(ncol=4,nrow=n_edge))
  colnames(edge_wise_s_f_r)<-c("edge","r","t","p")
  edge_wise_s_f_r[,1]<-intersect_edge
  for (j in 1:n_edge){
    s_f_r<-data.frame(structural=as.numeric(structural_subset[,j]),functional=as.numeric(functional_subset[,j]))
    cortest<-cor.test(as.numeric(s_f_r$structural),as.numeric(s_f_r$functional),method="pearson")
    r<-as.numeric(cortest$estimate)
    t<-as.numeric(cortest$statistic)
    p<-as.numeric(cortest$p.value)
    edge_wise_s_f_r[j,2:4]<-c(r,t,p)
#    ggscatter(s_f_r, x = "structural", y = "functional", 
#              add = "reg.line", conf.int = TRUE, 
#              cor.coef = TRUE, cor.method = "pearson",
#              xlab = "Structural Jackknife Pseudovalue", ylab = "Functional Correlation")
  }
  write.csv(edge_wise_s_f_r,file.path(dirname,"JKPV_FC_Edge.csv"),row.names=F)
  return(edge_wise_s_f_r)
}