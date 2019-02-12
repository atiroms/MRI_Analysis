
library(ggplot2)
library(GGally)


data(flea)
ggpairs(flea, columns = 2:4)


fig<-ggpairs(flea, columns = 2:4,
        lower=list(continuous=wrap("points",alpha=0.1,size=0.001)),
        title="Title test.")


ggsave("test.eps",fig,device=cairo_ps,
       path="C:/Users/atiro/GitHub/MRI_Analysis/Practices",dpi=300,height=10,width=10,limitsize=F)