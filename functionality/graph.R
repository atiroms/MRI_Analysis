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


#**************************************************
# Plot correlation matrix =========================
#**************************************************
plot_corrmat<-function(input,dict_roi,title){
  input<-data.frame(input)
  for(i in seq(ncol(input))){
    colnames(input)[i]<-as.character(dict_roi[which(dict_roi$ID_long==colnames(input)[i]),"label_proper"])
  }
  input<-rownames_to_column(input, "row")
  for(i in seq(nrow(input))){
    input$row[i]<-as.character(dict_roi[which(dict_roi$ID_long==input$row[i]),"label_proper"])
  }
  input_tidy<-gather(input,column,r,2:ncol(input))
  fig<-ggplot(input_tidy, aes(column, row)) +
    geom_tile(aes(fill = r)) +
    scale_fill_gradientn(colors = matlab.like2(100),name="r",limits=c(-1,1)) +
    scale_y_discrete(limits = rev(input$row)) +
    scale_x_discrete(limits = input$row, position="top") +
    ggtitle(title) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=700/ncol(input),angle = 90,vjust=0,hjust=0),
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