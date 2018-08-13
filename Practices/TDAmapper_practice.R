library(TDAmapper)
library(ggplot2)
library(igraph)


First.Example.data = data.frame( x=2*cos(0.5*(1:100)), y=sin(1:100) )
qplot(First.Example.data$x,First.Example.data$y)


First.Example.dist = dist(First.Example.data)


First.Example.mapper <- mapper1D(distance_matrix = First.Example.dist,
                               filter_values = First.Example.data$x,
                               num_intervals = 6,
                               percent_overlap = 50,
                               num_bins_when_clustering = 10)

#First.Example.mapper


First.Example.graph <- graph.adjacency(First.Example.mapper$adjacency, mode="undirected")


plot(First.Example.graph, layout = layout.auto(First.Example.graph))


First.Example.mapper <- mapper2D(
  distance_matrix = First.Example.dist,
  filter_values = list(First.Example.data$x,First.Example.data$y),
  num_intervals = c(6,6),
  percent_overlap = 50,
  num_bins_when_clustering = 10)
