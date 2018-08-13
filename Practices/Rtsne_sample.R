## calling the installed package
train<- read.csv(file.choose()) ## Choose the train.csv file downloaded from the link above  
library(Rtsne)
library(RColorBrewer)
## Curating the database for analysis with both t-SNE and PCA
Labels<-train$label
train$label<-as.factor(train$label)
## for plotting
colors <- brewer.pal(10, "Spectral")[as.integer(cut(as.integer(train$label), 10))]
#colors = rainbow(length(unique(train$label)))
#names(colors) = unique(train$label)

## Executing the algorithm on curated data
tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
#tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 5000)
#exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))

## Plotting
#plot(tsne$Y, t='n', main="tsne")
#text(tsne$Y, labels=train$label, col=colors[train$label])
plot(tsne$Y, pch=20, main="t-SNE", col=colors)