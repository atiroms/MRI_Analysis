#**************************************************
# Description =====================================
#**************************************************

# R script for common visualization functionalities in MRI analysis


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(ggraph)
library(igraph)
library(colorRamps)
library(purrr)


#**************************************************
# Plot PCA/ICA result =============================
#**************************************************

plot_ca<-function(df_src,list_name_covar,n_dim){
  df_plot<-df_src
  df_plot$ses<-as.factor(df_plot$ses)
  df_plot$ID_pnTTC<-as.factor(df_plot$ID_pnTTC)
  list_plot<-list()
  for (i_dim in 1:(n_dim-1)){
    list_plot_dim<-list()
    for (name_covar in list_name_covar){
      df_plot_subset<-df_plot[,c("ses","ID_pnTTC",name_covar,sprintf("dim_%02d",i_dim),sprintf("dim_%02d",i_dim+1))]
      colnames(df_plot_subset)<-c("ses","ID_pnTTC","color","x","y")
      plot<-(ggplot(df_plot_subset)
             + aes(x=x,y=y,label=ID_pnTTC,shape=ses)
             + scale_shape_manual(name=NULL,labels=c("1st wave","2nd wave"),values=c(3,4))
             + geom_point(size=2,aes(color=color))
             + scale_color_gradientn(colors = matlab.like2(100),name=name_covar)
             #+ geom_text_repel(size=2)
             + geom_path(aes(group=ID_pnTTC),size=0.5,alpha=0.5)
             #+ ggtitle("PCA of FC")
             + xlab(sprintf("Dimension %02d",i_dim))
             + ylab(sprintf("Dimension %02d",i_dim+1))
             + theme_light()
             + theme(plot.title = element_text(hjust = 0.5))
      )
      list_plot_dim_covar<-list(plot)
      names(list_plot_dim_covar)<-name_covar
      list_plot_dim<-c(list_plot_dim,list_plot_dim_covar)
      #ggsave(paste("atl-",atlas,"_dim-",sprintf("%02d",i_dim),"-",sprintf("%02d",i_dim+1),"_cov-",name_covar,"_pca_fc.eps",sep=""),plot=plot,device=cairo_ps,
      #       path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    }
    list_plot_dim<-list(list_plot_dim)
    names(list_plot_dim)<-as.character(i_dim)
    list_plot<-c(list_plot,list_plot_dim)
  }
  return(list_plot)
}


#**************************************************
# GAMM plot =======================================
#**************************************************
# modified from voxel/plotGAM

plot_gamm<-function(mod_gamm,df_join_measure_roi,spec_graph){
  #df_src <- mod_gamm$model
  key_df_src<-c("value",names(mod_gamm$var.summary))
  if ("ID_pnTTC" %in% key_df_src){
    key_df_src<-key_df_src
  }else{
    key_df_src<-c(key_df_src, "ID_pnTTC")
  }
  df_src <- df_join_measure_roi[key_df_src]
  
  plot<-ggplot()
  # add prediction line + ribbon to plot
  if (!is.null(spec_graph[["smooth"]])){
    for (name_smooth in names(spec_graph[["smooth"]])){
      spec_smooth<-spec_graph[["smooth"]][[name_smooth]]
      df_smooth <- data.frame(x = seq(min(df_src[spec_graph[["x_axis"]]]),
                                      max(df_src[spec_graph[["x_axis"]]]),
                                      length.out=200))
      names(df_smooth) <- spec_graph[["x_axis"]]
      for (i in names(df_src)) {
        if (i != spec_graph[["x_axis"]]) {
          if (any(class(df_src[i][,1])[1] == c("numeric", "integer","boolean"))) {
            df_smooth[, dim(df_smooth)[2] + 1] <- mean(df_src[i][,1])
            names(df_smooth)[dim(df_smooth)[2]] <- i
          }
          else if (any(class(df_src[i][,1])[1] == c("character", "factor","ordered"))) {
            df_smooth[, dim(df_smooth)[2] + 1] <- df_src[i][1,1]
            names(df_smooth)[dim(df_smooth)[2]] <- i
          }
        }
      }
      if (!is.null(spec_smooth[["fix"]])){
        for (var in names(spec_smooth[["fix"]])){
          df_smooth[[var]]<-spec_smooth[["fix"]][[var]]
        }
      }
      #df_smooth <- cbind(df_smooth,
      #                   as.data.frame(predict.gam(mod_gamm, df_smooth, se.fit = TRUE)))
      df_smooth <- cbind(df_smooth,
                         as.data.frame(predict.gam(mod_gamm, df_smooth, exclude="s(ID_pnTTC)",se.fit = TRUE)))
      plot <- (plot
               + geom_line(aes(x=!!df_smooth[,1],y=!!df_smooth[["fit"]]),
                           color=spec_smooth[["color"]],size=0.5,alpha=spec_smooth[["alpha"]]))
      if (spec_smooth[["ribbon"]]){
        plot <- (plot
                 + geom_ribbon(aes(x=!!df_smooth[,1],
                                   ymax = !!df_smooth[["fit"]]+1.96*!!df_smooth[["se.fit"]],
                                   ymin = !!df_smooth[["fit"]]-1.96*!!df_smooth[["se.fit"]],
                                   linetype=NA),
                               fill=spec_smooth[["color"]],alpha = 0.2*spec_smooth[["alpha"]]))
      }
    }
  }
  
  # add point + path to plot
  if (!is.null(spec_graph[["point"]])){
    df_point<-data.frame()
    for (name_point in names(spec_graph[["point"]])){
      spec_point<-spec_graph[["point"]][[name_point]]
      df_point_series<-df_src
      for (var_subset in names(spec_point[["subset"]])){
        df_point_series<-df_point_series[df_point_series[[var_subset]]==spec_point[["subset"]][[var_subset]],]
      }
      df_point_series$name_series<-name_point
      df_point_series$color<-spec_point[["color"]]
      df_point_series$alpha<-spec_point[["alpha"]]
      df_point<-rbind(df_point,df_point_series)
    }
    plot <- (plot
             + geom_point(aes(x=df_point[[spec_graph[["x_axis"]]]],
                              y=df_point[,1]),
                          color=df_point[["color"]],shape=1,size=2,alpha=0.4*df_point[["alpha"]]))
    if (any(duplicated(df_point[["ID_pnTTC"]]))){
      plot <- (plot
               + geom_path(aes(x=df_point[[spec_graph[["x_axis"]]]],
                               y=df_point[,1],
                               group=df_point[["ID_pnTTC"]]),
                           color=df_point[["color"]],size=0.3,alpha=0.2*df_point[["alpha"]],linetype="dashed"))
    }
  }
  
  # add themes
  plot <- (plot
           + theme_light()
           + theme(plot.title = element_text(hjust = 0.5),
                   legend.position = c("top","left")))
  
  return(plot)
}


#**************************************************
# Plot correlation matrix in heatmap ==============
#**************************************************
plot_cor_heatmap<-function(input){
  input_tidy<-rownames_to_column(input,"row")
  input_tidy<-gather(input_tidy,key=column,value=r,2:ncol(input_tidy))
  plot<-(ggplot(input_tidy, aes(column, row))
         + geom_tile(aes(fill = r))
         + scale_fill_gradientn(colors = matlab.like2(100),name="r",limits=c(-1,1))
         + scale_y_discrete(limits = rev(rownames(input)))
         + scale_x_discrete(limits = colnames(input), position="top")
         + theme_light()
         + theme(axis.text.x = element_text(size=700/ncol(input),angle = 90,vjust=0,hjust=0),
                 axis.text.y = element_text(size=700/ncol(input)),
                 #axis.title=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank()))
    return(plot)
}


#**************************************************
# Helper functions for ggpairs() ==================
#**************************************************
# corr + heatmap plot for upper half of correlogram
custom_corr_heatmap <- function(data, mapping, ...) {
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ct <- cor.test(x,y)
  r <- unname(ct$estimate)
  r_text <- format(r, digits=3)
  sig_text <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "** ", "*  ", "   ")
  )
  text<-paste(r_text,sig_text,sep=" ")
  
  
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(r, seq(-1, 1, length=100))]
  
  ggally_text(
    label = text, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 5,
    color = I("black"),
    ...
  ) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

# small line width, xlim and ylim, loess
custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("black"),alpha=0.01,size=0.001,stroke = 0, shape = ".") + 
    geom_smooth(method = "lm", color = I("red"),size=0.1,linetype="dotted",fill=I("pink"), ...) +
    geom_smooth(method = "loess", color = I("blue"),size=0.1,linetype="dotted",fill=I("lightblue"), ...) +
    xlim(-1,1) +
    ylim(-1,1)
}

# xlim and ylim, theme change
custom_densityDiag <- function(data, mapping, ...){
  
  mapping <- mapping_color_to_fill(mapping)
  
  ggplot(data, mapping) +
    theme_minimal() +
    scale_y_continuous() +
    xlim(-1,1) +
    geom_density(linetype="blank", fill=I("grey"),...)
}


#**************************************************
# circular graph ==================================
#**************************************************
graph_circular<-function(input,type_pvalue,thr_pvalue){
  edge_plot<-input$edge
  edge_plot<-edge_plot[which(edge_plot[,type_pvalue]<thr_pvalue),]
  node_plot<-input$node
  r_node<-grep("^R ",node_plot$label)
  l_node<-rev(grep("^L ",node_plot$label))
  if (length(r_node)+length(l_node)>0){
    node_plot<-rbind(node_plot[r_node,],node_plot[c(-r_node,-l_node),],node_plot[l_node,])
  }
  node_plot$angle <- 90 - 360 * ((1:nrow(node_plot))-0.5) / nrow(node_plot)
  node_plot$hjust<-ifelse(node_plot$angle < -90, 1, 0)
  node_plot$angle<-ifelse(node_plot$angle < -90, node_plot$angle+180, node_plot$angle)
  data_igraph <- graph_from_data_frame(d = edge_plot, vertices = node_plot, directed = F)
  #limit_color <- max(abs(max(edge_plot$weight)),abs(min(edge_plot$weight)))
  #limit_color <- c(-limit_color,limit_color)
  limit_color <- c(-1,1)
  fig<-ggraph(data_igraph, layout = "linear",circular = T) +
    geom_node_text(aes(x = x*1.03, y=y*1.03,
                       label=label, angle = angle, hjust=hjust,vjust=0.2),
                   size=400/nrow(node_plot), alpha=1) +
    geom_node_point(aes(x=x, y=y),size=1, alpha=1,colour="grey50") +
    scale_edge_color_gradientn(colors=matlab.like2(100),limits=limit_color,na.value="grey50")+
    expand_limits(x = c(-2, 2), y = c(-2, 2))+
    #ggtitle(input_title) +
    theme_void() +
    #theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1), legend.position=c(1,1))
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  if (nrow(edge_plot)>0){
    fig<-fig+
      geom_edge_arc(aes(color=weight),width=1,alpha=0.5)
  }
  return(fig)
}