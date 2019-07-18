
#### Circular Plotting ####

CircularPlot<-function(input,pvalue_type,input_title){
  subset_edges<-input$edges
  subset_edges<-subset_edges[which(subset_edges[,pvalue_type]<p_corrected),]
  ordered_nodes<-input$nodes
  nodes_label<-ConvertID(ordered_nodes$label,roi_data,"ID_long","label_proper")
  r_nodes<-grep("^R ",nodes_label)
  l_nodes<-rev(grep("^L ",nodes_label))
  ordered_nodes<-rbind(ordered_nodes[r_nodes,],ordered_nodes[c(-r_nodes,-l_nodes),],ordered_nodes[l_nodes,])
  ordered_nodes$angle <- 90 - 360 * ((1:nrow(ordered_nodes))-0.5) / nrow(ordered_nodes)
  ordered_nodes$hjust<-ifelse(ordered_nodes$angle < -90, 1, 0)
  ordered_nodes$angle<-ifelse(ordered_nodes$angle < -90, ordered_nodes$angle+180, ordered_nodes$angle)
  graph_data <- graph_from_data_frame(d = subset_edges, vertices = ordered_nodes, directed = F)
  color_limits <- max(abs(max(subset_edges$weight)),abs(min(subset_edges$weight)))
  color_limits <- c(-color_limits,color_limits)
  fig<-ggraph(graph_data, layout = "linear",circular = T) +
    geom_node_text(aes(x = x*1.03, y=y*1.03,
                       label=ConvertID(label,roi_data,"ID_long","label_proper"),
                       angle = angle, hjust=hjust,vjust=0.2),
                   size=2.5, alpha=1) +
    geom_node_point(aes(x=x, y=y),size=1, alpha=1,colour="grey50") +
    scale_edge_color_gradientn(colors=matlab.like2(100),limits=color_limits,na.value="grey50")+
    expand_limits(x = c(-2, 2), y = c(-2, 2))+
    ggtitle(input_title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1), legend.position=c(1,1))
  if (nrow(subset_edges)>0){
    fig<-fig+
      geom_edge_arc(aes(color=weight),width=1,alpha=0.5)
  }
  return(fig)
}


#### Basic correlation plotting ####

CorrelationPlot<-function(data_x,data_y,label_x,label_y){
  data<-data.frame(x=data_x,y=data_y)
  fig<-ggplot(data,aes(x=x, y=y)) +
    geom_point() +
    geom_smooth(method = "lm", colour="black", fill="grey50") +
    ggtitle(paste(label_y,"-",label_x,"Correlation",sep=" ")) +
    xlab(label_x) +
    ylab(label_y) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor=element_blank())
  return(fig)
}


#### Heatmap Plotting ####

HeatmapPlot<-function(input,title,xlabel,xelements,scale_data){
  colnames(input)[-1]<-xelements
  if(scale_data){
    input[,-1]<-scale(input[,-1])
    legend_title<-"Z-score"
  }else{
    legend_title<-NULL
  }
  tidyinput<-gather(input,ROI,measure,2:ncol(input))
  fig<-ggplot(tidyinput, aes(ROI, factor(ID_pnTTC))) +
    geom_tile(aes(fill = measure)) +
    scale_fill_gradientn(colors = matlab.like2(100),name=legend_title) +
    scale_y_discrete(limits = factor(rev(unique(tidyinput$ID_pnTTC)))) +
    scale_x_discrete(limits = xelements,position="top") +
    ggtitle(title) +
    xlab(xlabel) +
    ylab("pnTTC ID") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=5,angle = 90,vjust=0,hjust=0),
          axis.text.y = element_text(size=2,vjust=0), 
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(fig)
}


#### Correlation Matrix Plotting ####

CorrMatPlot<-function(input,title){
  input<-data.frame(input)
  colnames(input)<-ConvertID(colnames(input),roi_data,"ID_long","label_proper")
  input<-rownames_to_column(input, "row")
  input$row<-ConvertID(input$row,roi_data,"ID_long","label_proper")
  tidyinput<-gather(input,column,r,2:ncol(input))
  fig<-ggplot(tidyinput, aes(column, row)) +
    geom_tile(aes(fill = r)) +
    scale_fill_gradientn(colors = matlab.like2(100),name="r") +
    scale_y_discrete(limits = rev(input$row)) +
    scale_x_discrete(limits = input$row, position="top") +
    ggtitle(title) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=5,angle = 90,vjust=0.3,hjust=0),
          axis.text.y = element_text(size=5),
          axis.title=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(fig)
}