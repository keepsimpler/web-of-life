## read the web-of-life files of Bascompte et al.
weboflife = numeric()
for (i in 1:59) {
  networkname = paste('M_PL_', sprintf('%03d', i), sep='')
  weboflife[i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_PL_", sprintf("%03d", i), ".csv", sep=""), header=T, row.names = 1))
}
for (i in 1:30) {
  networkname = paste('M_SD_', sprintf('%03d', i), sep='')
  weboflife[59 + i] = networkname
  assign(networkname,
         read.csv(paste("~/Data/Bascompte/web-of-life/M_SD_", sprintf("%03d", i), ".csv", sep=""), header=T, row.names = 1))  
}

library(plyr)
library(bipartite)

res = ldply(weboflife, function(networkname) {
  network = get(networkname)
  network[network > 0] = 1
  edges = sum(network > 0)
  s1 = dim(network)[1]; s2 = dim(network)[2]
  ratio = s1 / s2
  nested.nodf = nested(network, method = 'NODF')
  #nested.nodf2 = nested(network, method = 'NODF2')
  network.adj = as.one.mode(network)
  degrees = rowSums(network.adj)
  network.adj.square = network.adj %*% network.adj
  degrees.square = rowSums(network.adj.square)
  c(heterogeneity = sum(degrees^2) / sum(degrees)^2, diversity = dim(network.adj)[1],
    disassortativity = sum(degrees.square / degrees), connectance = edges / (s1 * s2),
    k = edges / (s1 + s2) * ratio,
    nested.nodf = nested.nodf, ratio = ratio)
})
pairs(res)
