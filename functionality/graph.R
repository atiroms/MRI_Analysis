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
# GAMM plot =======================================
#**************************************************
# modified from voxel/plotGAM

plot_gamm<-function(mod_gamm,covar_x,color){
  df_src <- mod_gamm$model
  df_plot <- data.frame(x = seq(min(df_src[covar_x]),
                                max(df_src[covar_x]),
                                length.out=200))
  names(df_plot) <- covar_x
  for (i in names(df_src)[-1]) {
    if (i != covar_x) {
      if (any(class(df_src[i][,1])[1] == c("numeric", "integer","boolean"))) {
        df_plot[, dim(df_plot)[2] + 1] <- mean(df_src[i][,1])
        names(df_plot)[dim(df_plot)[2]] <- i
      }
      else if (any(class(df_src[i][,1])[1] == c("character", "factor","ordered"))) {
        df_plot[, dim(df_plot)[2] + 1] <- df_src[i][1,1]
        names(df_plot)[dim(df_plot)[2]] <- i
      }
    }
  }
  
  #for (i in 1:dim(df_plot)[2]) {
  #  if (class(df_plot[,i])[1] == "ordered" |  class(df_plot[,i])[1] == "factor") {
  #    warning("There are one or more factors in the model fit, please consider plotting by group since plot might be unprecise")
  #  }
  #}
  df_plot = cbind(df_plot, as.data.frame(predict.gam(mod_gamm, df_plot, se.fit = TRUE)))
  
  plot <- (ggplot(data=df_plot, aes(x=df_plot[,1]))
           + geom_line(aes(y=fit),
                       color=color,size=1)
           + geom_ribbon(data=df_plot, aes(ymax = fit+1.96*se.fit,
                                           ymin = fit-1.96*se.fit,
                                           linetype=NA),
                         fill=color,alpha = .3)
           + geom_point(data = df_src, aes(x=as_vector(df_src[covar_x]),
                                           y=df_src[,1]),
                        color=color, fill=color,size=3,alpha=.3)
           + geom_path(data = df_src, aes(x=as_vector(df_src[covar_x]),
                                          y=df_src[,1],
                                          group=as_vector(df_src["ID_pnTTC"])),
                       color=color,size=0.5,alpha=.2)
           #+ ggtitle("GAMM model")
           #+ ylab("Structural measure")
           #+ xlab("Tanner stage")
           + theme_light()
           + theme(plot.title = element_text(hjust = 0.5)))
  return(plot)
}


#**************************************************
# Plot correlation matrix in heatmap ==============
#**************************************************
cor_heatmap<-function(input){
  input_tidy<-gather(input,column,r,2:ncol(input))
  fig<-ggplot(input_tidy, aes(column, row)) +
    geom_tile(aes(fill = r)) +
    scale_fill_gradientn(colors = matlab.like2(100),name="r",limits=c(-1,1)) +
    scale_y_discrete(limits = rev(input$row)) +
    scale_x_discrete(limits = input$row, position="top") +
    #ggtitle(title) +
    theme_light() +
    theme(axis.text.x = element_text(size=700/ncol(input),angle = 90,vjust=0,hjust=0),
          axis.text.y = element_text(size=700/ncol(input)),
          axis.title=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
    return(fig)
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